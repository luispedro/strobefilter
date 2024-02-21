import os
from os import path
from jug import TaskGenerator

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

gmgc_v1 = get_gmgcv1()

