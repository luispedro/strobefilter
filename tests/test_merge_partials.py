from strobefilter import merge_sorted
def merge_chunks_in_mem(tempfiles):
    import numpy as np
    chunks = [np.load(f) for f in tempfiles]
    return np.unique(np.concatenate(chunks))


def test_merge_chunks(tmpdir):
    import numpy as np
    fs = [np.random.randint(0, 1024, 256) for i in range(16)]

    tempfiles = []
    for i, f in enumerate(fs):
        f.sort()
        f = np.unique(f)
        oname = f'{tmpdir}/chunk{i}.npy'
        np.save(oname, f)
        tempfiles.append(oname)
    r = merge_sorted(tempfiles, merge_chunksize=80)
    np.testing.assert_array_equal(merge_chunks_in_mem(tempfiles), r)


