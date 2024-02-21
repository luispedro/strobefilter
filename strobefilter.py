from fastaq import fastq_iter, fasta_iter
import strobealign
from collections import namedtuple

FilterResults = namedtuple('FilterResults', ['nr_reads', 'nr_strobes', 'nr_unigenes'])

def strobefilter_count(fqs, strategy='strict'):
    if strategy != 'strict':
        raise ValueError('Only strict strategy is implemented')
    ip = strobealign.IndexParameters.from_read_length(100)

    seen = set()
    n_fq = 0
    for ifile in fqs:
        for _, seq,_ in fastq_iter(ifile):
            n_fq += 1
            rs = strobealign.randstrobes_query(seq, ip)
            seen.update([rs.hash for rs in rs])

            if n_fq % 100_000 == 0:
                print(f'{n_fq//1000/1000.}m reads, {len(seen)//10000/100.}m hashes')

    n = 0
    s = 0
    for h, seq in fasta_iter(ffile):
        rs = strobealign.randstrobes_query(seq, ip)
        hs = set(rs.hash for rs in rs)
        n += 1
        has_any = bool(hs & seen)
        s += has_any
        if n % 100_000 == 0:
            print(f'{n//1000/1000.}m unigenes, {s/n:.5%} selected')
    return FilterResults(n_fq, s, len(seen))

