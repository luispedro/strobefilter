from fastaq import fastq_iter, fasta_iter
import strobealign
from collections import namedtuple

FilterResults = namedtuple('FilterResults', ['nr_unigenes_kept', 'strategy'])

def merge_sorted(tempfiles, merge_chunksize=256*1024*1024):
    import numpy as np
    chunks = [np.load(f, mmap_mode='r') for f in tempfiles]

    rs = []
    while chunks:
        min_val = min(ch[min(merge_chunksize, len(ch)-1)] for ch in chunks)
        next_chunks = []
        to_merge = []
        for ch in chunks:
            p = 0
            while p + 1024 < len(ch) and ch[p+1024] < min_val:
                p += 1024
            while p < len(ch) and ch[p] <= min_val:
                p += 1
            to_merge.append(ch[:p])
            if p < len(ch):
                next_chunks.append(ch[p:])
        to_merge = np.concatenate(to_merge)
        to_merge.sort()
        rs.append(np.unique(to_merge))
        chunks = next_chunks
    return np.concatenate(rs)

def extract_strobes(ifile, ip):
    import tempfile
    import numpy as np
    tempfiles = []
    next_ix = 0
    chunksize = 512_000_000
    chunk = np.zeros(chunksize, dtype=np.uint64)
    with tempfile.TemporaryDirectory() as tmpdir:
        n_fq = 0
        for _, seq,_ in fastq_iter(ifile):
            n_fq += 1
            rs = strobealign.randstrobes_query(seq, ip)
            for r in rs:
                chunk[next_ix] = r.hash
                next_ix += 1
                if next_ix == chunksize:
                    chunk.sort()
                    tname = f'{tmpdir}/chunk{len(tempfiles)}.npy'
                    np.save(tname, np.unique(chunk))
                    tempfiles.append(tname)
                    print(f'Wrote chunk {len(tempfiles)} to disk {tempfiles[-1]}')
                    next_ix = 0
            if n_fq % 1_000_000 == 0 and n_fq < 10_000_000 or n_fq % 10_000_000 == 0:
                print(f'{n_fq//1000/1000.}m reads')
        if next_ix > 0:
            chunk = np.unique(chunk[:next_ix])
            tname = f'{tmpdir}/chunk{len(tempfiles)}.npy'
            np.save(tname, chunk)
            tempfiles.append(tname)
            print(f'Wrote chunk {len(tempfiles)} to disk {tempfiles[-1]}')
        print(f'Merging {len(tempfiles)} chunks')
        r = merge_sorted(tempfiles)
        return (r, n_fq)

def extract_strobes_to(dataset, sample, ofile, preproc='25q/45ell'):
    from preprocess import preprocess
    import numpy as np
    from os import makedirs, path
    makedirs(path.dirname(ofile), exist_ok=True)
    ip = strobealign.IndexParameters.from_read_length(100)
    with preprocess(dataset, sample, preproc) as ifile:
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
    import random
    import numpy as np
    random.seed(123)
    sz = np.zeros(1, dtype=np.uint64)
    n = 0
    w = 0
    ofile = fn + '.subsampled'
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
    if strategy not in ['strict', 'packed']:
        raise ValueError('Only strict and packed strategies are implemented')
    seen = np.load(rmers)
    if strategy == 'packed':
        seen = seen.view(dtype=np.uint32).copy()
        seen.sort()
        packed = np.zeros(2**32, dtype=np.uint8)
        packed[seen] = 1
        del seen
    else:
        seen32 = seen.view(dtype=np.uint32)[::2]
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
