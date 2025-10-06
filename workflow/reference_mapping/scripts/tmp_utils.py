def get_model_path(model_path, tmpdir):
    from urllib.parse import urlparse
    import requests
    import shutil

    model_path = Path(model_path)

    # check if already an existing file path
    if model_path.exists():
        suffix = model_path.suffix
        assert suffix == '.pt', f'Incorrect file suffix {suffix}'
        return model_path

    # check if valid URL
    if not urlparse(str(model_path)).scheme in ('http', 'https', 'ftp'):
        raise ValueError(f'Invalid model path or URL: {model_path}')
    
    tmpdir = Path(tmpdir.name)
    
    # Download the model file to a temporary location
    model_path = tmpdir / model_path.name
    with requests.get(model_path, stream=True) as r:
        r.raise_for_status()
        with open(model_path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    if model_path.endswith('.zip'):
        with zipfile.ZipFile(model_path, 'r') as zip_ref:
            zip_ref.extractall(tmpdir.name)
        model_path = tmpdir / 'model.pt'
        assert model_path.exists(), f'Model file {model_path} not found in {tmpdir.name}'
    
    return model_path


with tempfile.TemporaryDirectory(dir=snakemake.tmpdir) as tmpdir:
    model_path = get_model_path(model_path, Path(tmpdir))
    # TODO: check for scPoli models, since it still uses an old API
