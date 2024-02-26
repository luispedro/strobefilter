from fastaq import fastq_iter, fasta_iter
import strobealign
from collections import namedtuple

FilterResults = namedtuple('FilterResults', ['nr_unigenes_kept', 'strategy'])

def extract_strobes(fqs, ip):
    import tempfile
    import subprocess
    import numpy as np
    import os
    import gzip
    with tempfile.TemporaryDirectory() as tmpdir:
        seen = set()
        CHUNK = 60_000_000
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

        for ifile in fqs:
            print(f'extracting strobes from {ifile}')
            for _, seq,_ in fastq_iter(ifile):
                n_fq += 1
                rs = strobealign.randstrobes_query(seq, ip)
                for r in rs: seen.add(r.hash)
                if len(seen) > CHUNK:
                    merge_tempfiles()
                    tmp_ix += 1
                if n_fq % 1_000_000 == 0 and n_fq < 10_000_000 or n_fq % 10_000_000 == 0:
                    print(f'{n_fq//1000/1000.}m reads, {len(seen)//10000/100.}m hashes')
        tempfile = merge_tempfiles()
        r = np.loadtxt(tempfile, dtype=np.uint64)
        r.sort()
        return r

def extract_strobes_to(fqs, ofile):
    import numpy as np
    from os import makedirs, path
    makedirs(path.dirname(ofile), exist_ok=True)
    ip = strobealign.IndexParameters.from_read_length(100)
    seen = extract_strobes(fqs, ip)
    np.save(ofile, seen)
    return ofile

def strobefilter_count(rmers, ffile, strategy='strict'):
    import numpy as np
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
        import stly
        print(f'Building unordered_set with {len(seen):,} elements')
        seens = stly.unordered_set_uint64_t()
        seens.reserve(len(seen))
        print(f'Allocated unordered_set with {seens.capacity():,} elements')
        for s in seen:
            seens.add(s)
        seen = seens
        print(f'Built unordered_set with {len(seen):,} elements')

    n = 0
    s = 0
    s1 = 0
    for h, seq in fasta_iter(ffile):
        n += 1
        rs = strobealign.randstrobes_query(seq, ip)
        hs = set(r.hash for r in rs)
        if strategy == 'packed':
            hs = np.array(list(hs), dtype=np.uint64)
            common = np.sum(
                        packed[
                            hs.view(dtype=np.uint32).reshape((-1, 2))
                            ].sum(1)
                        ==2)
        else:
            common = sum((h in seen) for h in hs)
        s  += common > 0
        s1 += common > 1
        if n % 1_000_000 == 0 and n < 10_000_000 or n % 10_000_000 == 0:
            print(f'{n//1000/1000.}m unigenes, {s/n:.5%} selected, {s1/n:.5%} with > 1 hit')
    return [FilterResults(s, 'min1')
            ,FilterResults(s1, 'min2')]

