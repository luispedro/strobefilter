import os
from os import path
from jug import TaskGenerator, Tasklet

import samples
from preprocess import preprocess
from strobefilter import strobefilter_count, extract_strobes_to

preprocess = TaskGenerator(preprocess)
strobefilter_count = TaskGenerator(strobefilter_count)
extract_strobes_to = TaskGenerator(extract_strobes_to)

GMGCV1_HASH = 'fde21071406072134befbbf6aacca6c9e27604a0d3e4954c2dbc3ee4bfe9dbb5'
GMGCV1_PATH = 'data/GMGC10.95nr.fna'
GMGCV1_URL = 'http://gmgc.embl.de/downloads/v1.0/GMGC10.95nr.fna'


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

results = {}
sizes = {}

# Jug 2.4 will support lambdas in Tasklets, but this is not yet released.
# So we need to declare this as a separate function.
def pp_strobe(p):
    return 'preproc-data/strobes/' + p + '.npy'

for study,ss in [
        (samples.DOG_STUDY, samples.DOG_SAMPLES),
        (samples.ZELLER_STUDY, samples.ZELLER_SAMPLES)]:
    ps = []
    for p in ss:
        pp = preprocess(study, p)
        ps.append(pp)
        rmers = extract_strobes_to([pp], Tasklet(pp, pp_strobe))
        sizes[p] = size_fastq(pp, rmers)
        for st in ['strict', 'packed']:
            strobe = strobefilter_count(rmers, gmgc_v1, strategy=st)
            results[(p, st)] = strobe
    #for st in ['strict']:
    #    results[study, st] = strobefilter_count(ps, gmgc_v1, strategy=st)



