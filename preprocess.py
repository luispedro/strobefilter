import subprocess
from os import makedirs, path

NGLESS_THREADS = 4

def preprocess(dataset, sample, compression: str = 'zst', preprocess: str = '25q/45ell') -> str:
    if preprocess != '25q/45ell':
        raise NotImplementedError('Only 25q/45ell is supported')

    oname = f'preproc-data/{dataset}/{sample}.preproc.fq.{compression}'
    makedirs(path.dirname(oname), exist_ok=True)
    subprocess.check_call(
            ['ngless',
             f'--threads={NGLESS_THREADS}',
             'preprocess.ngl',
             f'data/{dataset}/{sample}',
             oname])
    return oname

