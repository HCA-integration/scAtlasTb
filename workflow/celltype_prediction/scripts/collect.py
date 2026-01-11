from utils.io import read_anndata, write_zarr_linked

input_file = snakemake.input[0]
input_celltypist = snakemake.input.celltypist
input_predict_sex = snakemake.input.predict_sex
output_file = snakemake.output[0]

files_to_keep = ['obs']

print(f'Read file: {input_file}...', flush=True)
adata = read_anndata(input_file, obs='obs')

for file in input_celltypist:
    print(f'Read file: {file}...', flush=True)
    obs = read_anndata(file, obs='obs').obs
    adata.obs[obs.columns] = obs

for file in input_predict_sex:
    print(f'Read file: {file}...', flush=True)
    _ad = read_anndata(file, obs='obs', uns='uns')
    adata.obs[_ad.obs.columns] = _ad.obs
    adata.uns['predict_sex'] = _ad.uns['predict_sex']
    files_to_keep.append('uns')

print(f'Write file: {output_file}...', flush=True)
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
)