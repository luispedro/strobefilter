from fastaq import fastq_iter, fasta_iter
import strobealign
from collections import namedtuple

FilterResults = namedtuple('FilterResults', ['nr_reads', 'nr_strobes', 'nr_unigenes', 'strategy'])

def extract_strobes(fqs, ip):
    seen = set()
    n_fq = 0
    for ifile in fqs:
        print(f'Extracting strobes from {ifile}')
        for _, seq,_ in fastq_iter(ifile):
            n_fq += 1
            rs = strobealign.randstrobes_query(seq, ip)
            for r in rs: seen.add(r.hash)
            if n_fq % 1_000_000 == 0 and n_fq < 10_000_000 or n_fq % 10_000_000 == 0:
                print(f'{n_fq//1000/1000.}m reads, {len(seen)//10000/100.}m hashes')
    return seen

def extract_strobes_to(fqs, ofile):
    import numpy as np
    from os import makedirs, path
    makedirs(path.dirname(ofile), exist_ok=True)
    ip = strobealign.IndexParameters.from_read_length(100)
    seen = extract_strobes(fqs, ip)
    seena = np.array(list(seen), dtype=np.uint64)
    np.save(ofile, seena)
    return ofile

def strobefilter_count(rmers, ffile, strategy='strict'):
    import numpy as np
    if strategy != 'strict':
        raise ValueError('Only strict strategy is implemented')
    ip = strobealign.IndexParameters.from_read_length(100)
    seen = np.load(rmers, mmap_mode='r')
    seen = set(seen)

    n = 0
    s = 0
    s1 = 0
    for h, seq in fasta_iter(ffile):
        rs = strobealign.randstrobes_query(seq, ip)
        hs = set(rs.hash for rs in rs)
        n += 1
        common = sum((h in seen) for h in hs)
        s  += common > 0
        s1 += common > 1
        if n % 1_000_000 == 0 and n < 10_000_000 or n % 10_000_000 == 0:
            print(f'{n//1000/1000.}m unigenes, {s/n:.5%} selected, {s1/n:.5%} with > 1 hit')
    return [FilterResults(n_fq, len(seen), s, 'strict')
            ,FilterResults(n_fq, len(seen), s1, 'min2')]

