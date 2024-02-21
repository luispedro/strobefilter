import subprocess
from os import makedirs

NGLESS_THREADS = 4

def preprocess(dataset, sample, compression: str = 'zst') -> str:
    makedirs(f'preproc/{dataset}', exist_ok=True)

    oname = f'preproc/{dataset}/{sample}.preproc.fq.{compression}'
    subprocess.check_call(
            ['ngless',
             f'--threads={NGLESS_THREADS}',
             'preprocess.ngl',
             f'{dataset}/{sample}',
             oname])
    return oname

