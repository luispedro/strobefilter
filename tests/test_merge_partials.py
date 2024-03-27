from strobefilter import merge_chunks
def merge_chunks_in_mem(tempfiles):
    import numpy as np
    chunks = []
    for f in tempfiles:
        with open(f, 'rb') as g:
            chunks.append(np.frombuffer(g.read(), dtype=np.uint64))
    return np.unique(np.concatenate(chunks))


def test_merge_chunks(tmpdir):
    import numpy as np
    fs = [np.random.randint(0, 1024, 128) for i in range(16)]

    tempfiles = []
    for i, f in enumerate(fs):
        f.sort()
        f = np.unique(f)
        oname = f'{tmpdir}/chunk{i}.raw'
        with open(oname, 'wb') as g:
            g.write(f.data)
        tempfiles.append(oname)
    merge_chunks(tempfiles, f'{tmpdir}/merged.raw', merge_chunksize=16)
    np.testing.assert_array_equal(merge_chunks_in_mem(tempfiles), np.fromfile(f'{tmpdir}/merged.raw', dtype=np.uint64))


