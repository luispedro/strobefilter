def open_for(fname):
    if fname.endswith('.gz'):
        import gzip
        return gzip.open
    elif fname.endswith('.bz2'):
        import bz2
        return bz2.open
    elif fname.endswith('.xz'):
        import lzma
        return lzma.open
    elif fname.endswith('.zst') or fname.endswith('.zstd'):
        import zstandard as zstd
        return zstd.open
    else:
        return open

def fasta_iter(fname, full_header=False):
    '''Iterate over a (possibly gzipped) FASTA file

    Parameters
    ----------
    fname : str
        Filename.
            If it ends with .gz, gzip format is assumed
            If .bz2 then bzip2 format is assumed
            if .xz, then lzma format is assumerd
    full_header : boolean (optional)
        If True, yields the full header. Otherwise (the default), only the
        first word

    Yields
    ------
    (h,seq): tuple of (str, str)
    '''
    header = None
    chunks = []
    op = open_for(fname)
    with op(fname, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if header is not None:
                    yield header,''.join(chunks)
                line = line[1:].strip()
                if not line:
                    header = ''
                elif full_header:
                    header = line.strip()
                else:
                    header = line.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)

def fastq_iter(ifile):
    op = open_for(ifile)
    with op(ifile, 'rt') as f:
        while True:
            name = f.readline().strip()
            if not name:
                break
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            yield name, seq, qual

