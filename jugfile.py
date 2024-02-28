import os
from os import path
from jug import TaskGenerator

import samples
from strobefilter import strobefilter_count, extract_strobes_to, extract_fa_strobes, subsample_strobed_fasta
from print_report import print_report

strobefilter_count = TaskGenerator(strobefilter_count)
extract_strobes_to = TaskGenerator(extract_strobes_to)
extract_fa_strobes = TaskGenerator(extract_fa_strobes)
subsample_strobed_fasta = TaskGenerator(subsample_strobed_fasta)
print_report = TaskGenerator(print_report)

GMGCV1_HASH = 'fde21071406072134befbbf6aacca6c9e27604a0d3e4954c2dbc3ee4bfe9dbb5'
GMGCV1_PATH = 'data/GMGC10.95nr.fna'
GMGCV1_URL = 'http://gmgc.embl.de/downloads/v1.0/GMGC10.95nr.fna'

@TaskGenerator
def nr_elements_in_strobed_fasta(f):
    from strobefilter import read_strobed_fasta_chunk
    n = 0
    for _ in read_strobed_fasta_chunk(f):
        n += 1
    return n

def hash_file(file_path):
    import hashlib
    h = hashlib.sha256()
    with open(file_path, 'rb') as file:
        while True:
            chunk = file.read(h.block_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()

@TaskGenerator
def get_gmgcv1():
    import requests
    target = GMGCV1_PATH
    if path.exists(target):
        print(f'Target {target} exists. Double checking hash.')
        hd = hash_file(target)
        if hd == GMGCV1_HASH:
            print(f'Hash matches. Skipping download.')
            return target
        else:
            print(f'Hash mismatch. Cowardly refusing to overwrite. Expected {GMGCV1_HASH}, got {hd}.')
            raise Exception('Hash mismatch')
    print(f'Downloading {target}.')

    r = requests.get(GMGC_URL, allow_redirects=True, stream=True)

    with open(target, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    print(f'Finished downloading {target}.')
    if hash_file(target) != GMGCV1_HASH:
        os.unlink(target)
        raise Exception('Hash mismatch')
    return target

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

gmgc_v1 = get_gmgcv1()

fastrobes = extract_fa_strobes(gmgc_v1)
fastrobes_sub = subsample_strobed_fasta(fastrobes)

size_fasta = nr_elements_in_strobed_fasta(fastrobes)
size_fasta_sub = nr_elements_in_strobed_fasta(fastrobes_sub)

results = {}

@TaskGenerator
def reorganize_results(results):
    import pandas as pd
    reordered = []
    for (dataset,sample,strat),vs in results.items():
        for (n,crit) in vs:
            reordered.append((dataset+'/'+sample, strat, crit, n))
    reordered = pd.DataFrame(reordered, columns=['sample', 'strategy', 'criterion', 'nr_unigenes_kept'])
    reordered = pd.pivot(reordered, index='sample', values='nr_unigenes_kept', columns=['strategy', 'criterion'])
    return reordered

for dataset,ss in [
        (samples.DOG_STUDY, samples.DOG_SAMPLES),
        (samples.ZELLER_STUDY, samples.ZELLER_SAMPLES)]:
    ps = []
    for sample in ss:
        rmers = extract_strobes_to(dataset,
                                    sample,
                                    path.join('preproc-data',
                                                'strobes',
                                                dataset,
                                                f'{sample}.preproc.npy'))
        ps.append(rmers[0])
        for st in ['strict', 'packed']:
            if len(ps) < 5:
                results[dataset, sample, st+'-full'] = strobefilter_count(rmers[0], fastrobes, strategy=st)
            results[dataset, sample, st] = strobefilter_count(rmers[0], fastrobes_sub, strategy=st)
    #for st in ['strict']:
    #    results[dataset, st] = strobefilter_count(ps, fastrobes, strategy=st)

final = reorganize_results(results)
print_report(final, size_fasta_sub, 'output_report.md')
