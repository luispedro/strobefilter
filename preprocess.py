import subprocess
from contextlib import contextmanager
import threading
NGLESS_THREADS = 4

def stream_gzips_to(files, fd):
    for f in files:
        assert f.endswith('.gz')
        print(f'Unzipping {f}')
        subprocess.check_call(['zcat', f], stdout=fd)

@contextmanager
def preprocess(dataset, sample, preprocess: str = '25q/45ell') -> str:
    if preprocess not in ['25q/45ell', 'passthru']:
        raise NotImplementedError(f'Only 25q/45ell or passthru supported (got {preprocess})')

    if preprocess == '25q/45ell':
        script = 'preprocess.ngl'
        proc = subprocess.Popen(
                ['ngless',
                 f'--threads={NGLESS_THREADS}',
                 script,
                 f'data/{dataset}/{sample}'],
                stdout=subprocess.PIPE,
                text=True,
                )
        yield proc.stdout
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f'Preprocessing failed for {dataset}/{sample}')
    else:
        from glob import glob
        import os
        fs = glob(f'data/{dataset}/{sample}/*')
        for f in fs:
            assert f.endswith('.gz')
        read_fd, write_df = os.pipe()
        read_f = os.fdopen(read_fd, 'r')
        write_f = os.fdopen(write_df, 'w')
        t = threading.Thread(target=stream_gzips_to, args=(fs, write_f))
        t.start()
        yield read_f
        write_f.close()
        read_f.close()
        t.join()

