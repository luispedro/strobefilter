import os
from os import path
from jug import TaskGenerator

import samples
from strobefilter import strobefilter_count, extract_strobes_to, extract_fa_strobes, subsample_strobed_fasta
from print_report import print_report
from download import download_if_needed

strobefilter_count = TaskGenerator(strobefilter_count)
extract_strobes_to = TaskGenerator(extract_strobes_to)
extract_fa_strobes = TaskGenerator(extract_fa_strobes)
subsample_strobed_fasta = TaskGenerator(subsample_strobed_fasta)
print_report = TaskGenerator(print_report)
download_if_needed = TaskGenerator(download_if_needed)

GMGCV1_HASH = 'fde21071406072134befbbf6aacca6c9e27604a0d3e4954c2dbc3ee4bfe9dbb5'
GMGCV1_PATH = 'data/GMGC10.95nr.fna'
GMGCV1_URL = 'http://gmgc.embl.de/downloads/v1.0/GMGC10.95nr.fna'

SUBCATALOGS = {
        'dog-gut':
            (
                'http://gmgc.embl.de/downloads/v1.0/subcatalogs/GMGC10.dog-gut.95nr.fna.gz',
                'data/GMGC10.dog-gut.95nr.fna.gz',
                '2b4de550829f8b8ddbb77c68de6e313964da413a6603b9ec3f470f7365ccb36b'
            ),
        'human-gut':
            (
                'http://gmgc.embl.de/downloads/v1.0/subcatalogs/GMGC10.human-gut.95nr.fna.gz',
                'data/GMGC10.human-gut.95nr.fna.gz',
                'ef48089b9b7e6c0bfa3d01108c471d8575bcfed50dc392a8d40cdcc9e48a8b2f',
            ),
        'marine':
            (
                'http://gmgc.embl.de/downloads/v1.0/subcatalogs/GMGC10.marine.95nr.fna.gz',
                'data/GMGC10.marine.95nr.fna.gz',
                'a9ba24184631c5190012dc26f1b0557c0806254900e2e1824a9dc81b30e8336a',
            ),
            }


@TaskGenerator
def nr_elements_in_strobed_fasta(f):
    from strobefilter import read_strobed_fasta_chunk
    n = 0
    for _ in read_strobed_fasta_chunk(f):
        n += 1
    return n


@TaskGenerator
def get_gmgcv1():
    return download_if_needed.f(GMGCV1_URL, GMGCV1_PATH, GMGCV1_HASH)

@TaskGenerator
def size_fastq(fastq, strobef : str):
    print(f'size_fastq({fastq}, {strobef})')
    from fastaq import fastq_iter
    import numpy as np
    n = 0
    for _ in fastq_iter(fastq):
        n += 1
    rmers = np.load(strobef, mmap_mode='r')
    return n, rmers.shape[0]

@TaskGenerator
def merge_partials(partials, oname):
    from strobefilter import merge_sorted
    import numpy as np
    r = merge_sorted(partials)
    np.save(oname, r)
    return oname

gmgc_v1 = get_gmgcv1()
fastrobes = extract_fa_strobes(gmgc_v1)
fastrobes_sub = subsample_strobed_fasta(fastrobes)

subcatalog = {}
nr_elements = {}
for key, (url, opath, hv) in SUBCATALOGS.items():
    fafile = download_if_needed(url, opath, hv)
    cat_fastrobes = extract_fa_strobes(fafile)
    subcatalog[key] = subsample_strobed_fasta(cat_fastrobes)
    nr_elements[key] = nr_elements_in_strobed_fasta(subcatalog[key])

size_fasta = nr_elements_in_strobed_fasta(fastrobes)
size_fasta_sub = nr_elements_in_strobed_fasta(fastrobes_sub)

results = {}

@TaskGenerator
def reorganize_results(results):
    import pandas as pd
    reordered = []
    for (dataset,preproc,sample,strat),(vs, n_db) in results.items():
        for (n,crit) in vs:
            reordered.append((dataset+'/'+sample, preproc, strat, crit, n, n_db))
    reordered = pd.DataFrame(reordered, columns=['sample', 'preproc', 'strategy', 'criterion', 'nr_unigenes_kept', 'initial_db'])
    reordered.eval('frac_kept = nr_unigenes_kept/initial_db', inplace=True)
    reordered = pd.pivot(reordered, index='sample', values='frac_kept', columns=['strategy', 'preproc', 'criterion'])
    return reordered

for dataset,ss,habitat in [
        (samples.ZELLER_STUDY, samples.ZELLER_SAMPLES, 'human-gut'),
        (samples.DOG_STUDY, samples.DOG_SAMPLES, 'dog-gut'),
        (samples.TARA_STUDY, samples.TARA_SAMPLES, 'marine'),
        ]:
    for preproc in ['25q-45ell', 'passthru']:
        ps = []
        for sample in ss:
            rmers = extract_strobes_to(dataset,
                                    sample,
                                    path.join('preproc-data',
                                                'strobes',
                                                dataset,
                                                f'{sample}.{preproc}.npy'),
                                    preproc=preproc)
            ps.append(rmers[0])
            for st in ['strict']:
                results[dataset, preproc, sample, st] = (strobefilter_count(rmers[0], fastrobes_sub, strategy=st),
                                                size_fasta_sub)
                results[f'{dataset}:{habitat}', preproc, sample, st] = (strobefilter_count(rmers[0], subcatalog[habitat], strategy=st),
                                                          nr_elements[habitat])
        combined = merge_partials(ps,
                       path.join('preproc-data',
                                 'strobes',
                                 dataset,
                                 f'{dataset}-combined.{preproc}.npy'))
        results[f'{dataset}-combined' , preproc, '0', st] = (strobefilter_count(combined, fastrobes_sub, strategy=st),
                                        size_fasta_sub)
        results[f'{dataset}-combined:{habitat}', preproc, sample, st] = (strobefilter_count(combined, subcatalog[habitat], strategy=st),
                                                  nr_elements[habitat])

    #for st in ['strict']:
    #    results[dataset, st] = strobefilter_count(ps, fastrobes, strategy=st)

final = reorganize_results(results)
print_report(final, 'output_report.md')
