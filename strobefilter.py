from fastaq import fastq_iter, fasta_iter
import strobealign
from collections import namedtuple

FilterResults = namedtuple('FilterResults', ['nr_unigenes_kept', 'strategy'])

def extract_strobes(fqs, ip):
    import tempfile
    import subprocess
    import numpy as np
    with tempfile.TemporaryDirectory() as tmpdir:
        tempfile = f'{tmpdir}/tmp.txt'
        seen = set()
        n_fq = 0
        p_n = 0
        with open(tempfile, 'w') as f:
            for ifile in fqs:
                print(f'extracting strobes from {ifile}')
                for _, seq,_ in fastq_iter(ifile):
                    n_fq += 1
                    rs = strobealign.randstrobes_query(seq, ip)
                    for r in rs: seen.add(r.hash)
                    if len(seen) > 15_000_000:
                        print(f'writing {len(seen)} hashes (from {n_fq - p_n} reads) to {tempfile}')
                        p_n = n_fq
                        for s in seen:
                            f.write(f'{s}\n')
                        seen.clear()
                    if n_fq % 1_000_000 == 0 and n_fq < 10_000_000 or n_fq % 10_000_000 == 0:
                        print(f'{n_fq//1000/1000.}m reads, {len(seen)//10000/100.}m hashes')
            for s in seen:
                f.write(f'{s}\n')
            seen.clear()
        print(f'sorting and removing duplicates from {tempfile}')
        subprocess.check_call(['sort', '-u', '-o', tempfile+'.u', tempfile])
        r = np.loadtxt(tempfile+'.u', dtype=np.uint64)
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
        seen = set(seen)

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
    return [FilterResults(s, 'strict')
            ,FilterResults(s1, 'min2')]

