import subprocess
from contextlib import contextmanager
NGLESS_THREADS = 4

@contextmanager
def preprocess(dataset, sample, preprocess: str = '25q/45ell') -> str:
    if preprocess != '25q/45ell':
        raise NotImplementedError('Only 25q/45ell is supported')

    proc = subprocess.Popen(
            ['ngless',
             f'--threads={NGLESS_THREADS}',
             'preprocess.ngl',
             f'data/{dataset}/{sample}'],
            stdout=subprocess.PIPE,
            text=True,
            )
    yield proc.stdout
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f'Preprocessing failed for {dataset}/{sample}')

