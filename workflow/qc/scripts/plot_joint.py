from pathlib import Path
from pandas.api.types import is_numeric_dtype
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
# from matplotlib.axes import Axes
import seaborn as sns
from tqdm import tqdm
import traceback
from joblib import Parallel, delayed
import scanpy as sc
from pprint import pformat
import logging

logging.basicConfig(level=logging.INFO)
sns.set_theme(style='white')
sc.set_figure_params(frameon=False, fontsize=10, dpi_save=200, vector_friendly=True)

from utils.io import read_anndata
from qc_utils import parse_parameters, get_thresholds, plot_qc_joint


input_zarr = snakemake.input.zarr
output_joint = snakemake.output.joint

output_joint = Path(output_joint)
output_joint.mkdir(parents=True, exist_ok=True)

logging.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, obs='obs', uns='uns')
print(adata, flush=True)

# get parameters
file_id = snakemake.wildcards.file_id
threads = snakemake.threads
dataset, hues = parse_parameters(adata, snakemake.params, filter_hues=True)
scautoqc_metrics = snakemake.params.get('scautoqc_metrics')
hues = list(set(hues+scautoqc_metrics+['qc_status']))


# if no cells filtered out, save empty plots
if adata.obs.shape[0] == 0:
    logging.info('No data, skip plotting...')
    exit()


def create_figure(df, png_file, density_png, joint_title, **kwargs):
    g = plot_qc_joint(df, **kwargs)
    
    # remove legend from joint plot image
    if g.ax_joint.legend_ is not None:
        g.ax_joint.get_legend().remove()
    
    # save plot temporarily
    plt.tight_layout()
    plt.savefig(png_file, bbox_inches='tight')
    plt.close('all')
    
    # assemble figure and update plot
    f, axes = plt.subplots(1, 2, figsize=(20, 10))
    axes[0].imshow(mpimg.imread(png_file))
    axes[1].imshow(mpimg.imread(density_png))

    # move legend from joint plot to right of figure
    handles, labels = g.ax_joint.get_legend_handles_labels()
    markerscale = (80 / kwargs.get('s', 20)) ** 0.5
    g.ax_joint.legend(markerscale=markerscale)
    if handles and labels:
        axes[0].legend(
            handles=handles,
            labels=labels,
            markerscale=markerscale,
            loc='center right',
            bbox_to_anchor=(0, 0.5),
            borderaxespad=0.5,
        )
    
    for ax in axes.ravel():
        ax.set_axis_off()
    plt.suptitle(joint_title, fontsize=16)
    
    # save final plot
    plt.tight_layout(w_pad=0.05)
    plt.savefig(png_file, bbox_inches='tight')
    plt.close('all')


def call_plot(df, x, y, log_x, log_y, hue, scatter_plot_kwargs, density_png, density_log_png):
    # logging.info(f'Joint QC plots for hue={hue}...') 
    joint_title = f'Joint QC for\n{dataset}\nmargin hue: {hue}'
       
    plot_path = output_joint / f'hue={hue}'
    plot_path.mkdir(exist_ok=True)

    # determine plotting parameters
    if is_numeric_dtype(df[hue]):
        palette = 'plasma'
        legend = 'brief'
    else:
        palette = None # if df[hue].nunique() < 50 else 'plasma'
        legend = df[hue].nunique() <= 30
    
    scatter_plot_kwargs |= dict(
        palette=palette,
        legend=legend,
        marginal_kwargs=dict(
            palette=palette,
            legend=False,
            stat='density',
        ),
    )
    
    # plot joint QC on regular scale
    create_figure(
        df,
        png_file=plot_path / f'{x}_vs_{y}.png',
        density_png=density_png,
        joint_title=joint_title,
        x=x,
        y=y,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title='',
        **scatter_plot_kwargs,
    )

    # plot in log scale 
    log_x_prefix = f'log_{log_x}_' if log_x > 1 else ''
    log_y_prefix = f'log_{log_y}_' if log_y > 1 else ''
    
    create_figure(
        df,
        png_file=plot_path / f'{log_x_prefix}{x}_vs_{log_y_prefix}{y}.png',
        density_png=density_log_png,
        joint_title=joint_title,
        x=x,
        y=y,
        log_x=log_x,
        log_y=log_y,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title='',
        **scatter_plot_kwargs,
    )


thresholds = get_thresholds(
    threshold_keys=scautoqc_metrics,
    autoqc_thresholds=adata.uns['scautoqc_ranges'],
    user_thresholds=snakemake.params.get('thresholds'),
)
logging.info(f'\n{pformat(thresholds)}')

scatter_plot_kwargs = dict(
    s=8,
    alpha=.8,
    linewidth=0,
)

kde_plot_kwargs = dict(
    fill=True,
    cmap='plasma',
    alpha=.8,
)
if adata.n_obs > 5e4:
    kde_plot_kwargs |= dict(
        bw_adjust=4,      # smoother for large data
        gridsize=75,      # moderate grid
        thresh=0.2,       # ignore more low-density
        levels=5,         # fewer contours
    )

coordinates = [
    ('n_counts', 'n_genes', 10, 10),
    ('n_genes', 'percent_mito', 2, 1),
    ('n_genes', 'percent_ribo', 2, 1),
    ('n_genes', 'percent_hb', 2, 1),
]
# filter to configured metrics
coordinates = [
    c for c in coordinates if
    all(x in scautoqc_metrics for x in c[:2])
]

# # subset to max of 300k cells due to high computational cost
density_data = adata.obs.sample(n=int(min(300_000, adata.n_obs)), random_state=42)
# density_data = adata.obs

for x, y, log_x, log_y in coordinates:
    logging.info(f'Joint QC plots per {x} vs {y}...')
    
    # temporary files
    density_png = output_joint / f'{x}_vs_{y}_density.png'
    density_log_png = output_joint / f'log_{x}_vs_{y}_density.png'
    
    try:
        logging.info('Plotting density...')
        plot_qc_joint(
            density_data,
            x=x,
            y=y,
            # main_plot_function=Axes.hexbin,
            main_plot_function=sns.kdeplot,
            x_threshold=thresholds[x],
            y_threshold=thresholds[y],
            title='',
            **kde_plot_kwargs,
        )
    except ValueError as e:
        logging.error(f'Error in plotting density: {e}')
        traceback.print_exc()
        
        logging.info('Retry with adjusted parameters...')
        kde_plot_kwargs.pop('bw_adjust', None)
        kde_plot_kwargs.pop('gridsize', None)
        
        plot_qc_joint(
            density_data,
            x=x,
            y=y,
            # main_plot_function=Axes.hexbin,
            main_plot_function=sns.kdeplot,
            x_threshold=thresholds[x],
            y_threshold=thresholds[y],
            title='',
            **kde_plot_kwargs,
        )
    plt.tight_layout()
    plt.savefig(density_png, bbox_inches='tight')   
    
    try:
        logging.info('Plotting density for log scale...')
        plot_qc_joint(
            density_data,
            x=x,
            y=y,
            main_plot_function=sns.kdeplot,
            log_x=log_x,
            log_y=log_y,
            x_threshold=thresholds[x],
            y_threshold=thresholds[y],
            title='',
            **kde_plot_kwargs,
        )
    except ValueError as e:
        logging.error(f'Error in plotting density for log scale: {e}')
        traceback.print_exc()
        
        logging.info('Retry with adjusted parameters...')
        kde_plot_kwargs.pop('bw_adjust', None)
        kde_plot_kwargs.pop('gridsize', None)
        
        plot_qc_joint(
            density_data,
            x=x,
            y=y,
            main_plot_function=sns.kdeplot,
            log_x=log_x,
            log_y=log_y,
            x_threshold=thresholds[x],
            y_threshold=thresholds[y],
            title='',
            **kde_plot_kwargs,
        )
    plt.tight_layout()
    plt.savefig(density_log_png, bbox_inches='tight')

    
    def safe_call_plot(*args, **kwargs):
        try:
            return call_plot(*args, **kwargs)
        except Exception as e:
            logging.error(f"Error in plot job: {e}")
            traceback.print_exc()
            return None
    
    Parallel(n_jobs=threads)(
        delayed(safe_call_plot)(
            df=adata.obs,
            x=x,
            y=y, 
            log_x=log_x,
            log_y=log_y,
            hue=hue,
            scatter_plot_kwargs=scatter_plot_kwargs,
            density_png=density_png,
            density_log_png=density_log_png
        ) for hue in tqdm(hues)
    )
    
    # remove redundant plots
    density_png.unlink()
    density_log_png.unlink()
