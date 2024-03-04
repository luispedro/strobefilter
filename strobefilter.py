from fastaq import fastq_iter, fasta_iter
import strobealign
from collections import namedtuple

FilterResults = namedtuple('FilterResults', ['nr_unigenes_kept', 'strategy'])

def extract_strobes(ifile, ip):
    import tempfile
    import numpy as np
    import os
    import gzip
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            import stly
            seen = stly.unordered_set_uint64_t()
            chunksize = 120_000_000
            seen.reserve(chunksize)
        except ImportError:
            seen = set()
            chunksize = 80_000_000
        n_fq = 0
        tmp_ix = 0
        def merge_tempfiles():
            print(f'Merging {len(seen)//10000/100.}m hashes')
            tempfile = f'{tmpdir}/tmp_{tmp_ix}.txt.gz'
            with gzip.open(tempfile, 'wt', compresslevel=1) as f:
                for s in seen:
                    f.write(f'{s}\n')
                if tmp_ix > 0:
                    print('Loading from previous chunk')
                    with gzip.open(f'{tmpdir}/tmp_{tmp_ix-1}.txt.gz', 'rt') as f2:
                        n_total = 0
                        n_written = 0
                        for line in f2:
                            ell = int(line.strip())
                            n_total += 1
                            if ell not in seen:
                                n_written += 1
                                f.write(line)
                    print(f'Copied {n_written//10000/100.}m hashes ({n_written/n_total:.2%}) from previous chunk, total {(n_written+len(seen))//10000/100.}m [{n_fq//1000/1000.}m reads; chunk {tmp_ix}]')
                    os.unlink(f'{tmpdir}/tmp_{tmp_ix-1}.txt.gz')
            seen.clear()
            return tempfile

        for _, seq,_ in fastq_iter(ifile):
            n_fq += 1
            rs = strobealign.randstrobes_query(seq, ip)
            for r in rs: seen.add(r.hash)
            if len(seen) > chunksize:
                merge_tempfiles()
                tmp_ix += 1
            if n_fq % 1_000_000 == 0 and n_fq < 10_000_000 or n_fq % 10_000_000 == 0:
                print(f'{n_fq//1000/1000.}m reads, {len(seen)//10000/100.}m hashes')
        if tmp_ix > 0:
            tempfile = merge_tempfiles()
            del seen
            r = np.loadtxt(tempfile, dtype=np.uint64)
        else:
            r = np.zeros(len(seen), dtype=np.uint64)
            for i, s in enumerate(seen):
                r[i] = s
        r.sort()
        return (r, n_fq)

def extract_strobes_to(dataset, sample, ofile):
    from preprocess import preprocess
    import numpy as np
    from os import makedirs, path
    makedirs(path.dirname(ofile), exist_ok=True)
    ip = strobealign.IndexParameters.from_read_length(100)
    with preprocess(dataset, sample) as ifile:
        seen, n_fq = extract_strobes(ifile, ip)
    np.save(ofile, seen)
    return ofile, n_fq, len(seen)

def read_strobed_fasta_chunk(fname):
    import numpy as np
    with open(fname, 'rb') as f:
        partial = None
        while ch := f.read(8*1024*1024):
            buf = np.frombuffer(ch, dtype=np.uint64)
            cur_p = 0
            if partial is not None:
                ch_s = int(partial[0])
                cur_p = ch_s - len(partial) + 1
                chunk = np.concatenate([partial[1:], buf[:cur_p]])
                yield chunk
                partial = None
            while True:
                ch_s = int(buf[cur_p])
                if cur_p + ch_s + 1 >= len(buf):
                    partial = buf[cur_p:]
                    break
                yield buf[cur_p+1:cur_p+ch_s+1]
                cur_p += ch_s + 1
                if cur_p >= len(buf):
                    break

def subsample_strobed_fasta(fn):
    ofile = fn + '.subsampled'
    import random
    import numpy as np
    random.seed(123)
    sz = np.zeros(1, dtype=np.uint64)
    n = 0
    w = 0
    with open(ofile, 'wb') as f:
        for ch in read_strobed_fasta_chunk(fn):
            n += 1
            if random.random() < .01:
                sz[0] = len(ch)
                f.write(sz.data)
                f.write(ch.data)
                w += 1
            if n % 1_000_000 == 0 and n < 10_000_000 or n % 10_000_000 == 0:
                print(f'{n//1000/1000.}m chunks, {w//1000/1000.}m written')
    return ofile


def strobefilter_count(rmers, preprocfa, strategy='strict'):
    import numpy as np
    import zstandard as zstd
    if strategy not in ['strict', 'packed']:
        raise ValueError('Only strict and packed strategies are implemented')
    ip = strobealign.IndexParameters.from_read_length(100)
    seen = np.load(rmers, mmap_mode='r')
    if strategy == 'packed':
        seen = seen.view(dtype=np.uint32).copy()
        seen.sort()
        packed = np.zeros(2**32, dtype=np.uint8)
        packed[seen] = 1
        del seen
    else:
        seen32 = seen & 0xFFFFFFFF
        seen32.sort()
        pref = np.zeros(2**32, dtype=bool)
        pref[seen32] = True
        del seen32

    n = 0
    s = 0
    s1 = 0
    for hs in read_strobed_fasta_chunk(preprocfa):
        n += 1
        if strategy == 'packed':
            common = np.sum(
                        packed[
                            hs.view(dtype=np.uint32).reshape((-1, 2))
                            ].sum(1)
                        ==2)
        else:
            hs = hs[pref[hs.view(dtype=np.uint32)[::2]]]
            ps = np.searchsorted(seen, hs)
            ps %= len(seen) # when the hash is larger than the largest in seen, this would be out of bounds
            common = np.sum(seen[ps] == hs)
        s  += common > 0
        s1 += common > 1
        if n % 1_000_000 == 0 and n < 10_000_000 or n % 10_000_000 == 0:
            print(f'{n//1000/1000.}m unigenes, {s/n:.5%} selected, {s1/n:.5%} with > 1 hit [{strategy}]')
    return [FilterResults(s, 'min1')
            ,FilterResults(s1, 'min2')]


def extract_fa_strobes(fafile):
    from os import makedirs, path
    import numpy as np
    ofile = 'preproc-data/preproc-fasta/' + path.basename(fafile) + '.bin'
    makedirs(path.dirname(ofile), exist_ok=True)
    ip = strobealign.IndexParameters.from_read_length(100)
    size = np.zeros(1, dtype=np.uint64)
    with open(ofile, "wb") as f:
        n = 0
        for h, seq in fasta_iter(fafile):
            rs = strobealign.randstrobes_query(seq, ip)
            hs = set(rs.hash for rs in rs)
            hs = np.array(list(hs), dtype=np.uint64)
            hs.sort()
            size[0] = len(hs)
            f.write(size.data)
            f.write(hs.data)
            n += 1
            if n % 1_000_000 == 0 and n < 10_000_000 or n % 10_000_000 == 0:
                print(f'{n//1000/1000.}m preprocessed ({fafile})')
    return ofile
