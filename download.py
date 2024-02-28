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

def download_if_needed(url, target, expected_hash):
    from os import path
    import requests
    if path.exists(target):
        print(f'Target {target} exists. Double checking hash.')
        hd = hash_file(target)
        if hd == expected_hash:
            print(f'Hash matches. Skipping download.')
            return target
        else:
            print(f'Hash mismatch. Cowardly refusing to overwrite. Expected {expected_hash}, got {hd}.')
            raise Exception('Hash mismatch')
    print(f'Downloading {target}.')

    r = requests.get(url, allow_redirects=True, stream=True)

    with open(target, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    print(f'Finished downloading {target}.')
    if hash_file(target) != expected_hash:
        os.unlink(target)
        raise Exception('Hash mismatch')
    return target

