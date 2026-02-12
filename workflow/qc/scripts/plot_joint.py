from pathlib import Path
import warnings
from pandas.api.types import is_numeric_dtype
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import seaborn as sns
import numpy as np
from tqdm import tqdm
import traceback
from joblib import Parallel, delayed
import scanpy as sc
from pprint import pformat
import logging

logging.basicConfig(level=logging.INFO)
warnings.filterwarnings(
    'ignore',
    message=r'.*pkg_resources is deprecated as an API.*',
    category=UserWarning,
)
sns.set_theme(style='white')
sc.set_figure_params(frameon=False, fontsize=10, dpi_save=200, vector_friendly=True)

from utils.io import read_anndata
from qc_utils import parse_parameters, get_thresholds, plot_qc_joint, QC_FLAGS


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
scautoqc_metrics = snakemake.params.get('scautoqc_metrics', [])
hues = list(dict.fromkeys(hues + QC_FLAGS + ['qc_status']))

logging.info(f'{hues=}')


# if no cells filtered out, save empty plots
if adata.obs.shape[0] == 0:
    logging.info('No data, skip plotting...')
    exit()


def _render_figure_to_image(fig):
    canvas = FigureCanvas(fig)
    canvas.draw()
    image = np.asarray(canvas.buffer_rgba())
    return image


def create_figure(df, png_file, density_img, joint_title, **kwargs):
    g = plot_qc_joint(df, **kwargs)

    handles, labels = g.ax_joint.get_legend_handles_labels()
    if g.ax_joint.legend_ is not None:
        g.ax_joint.get_legend().remove()

    joint_img = _render_figure_to_image(g.fig)
    plt.close(g.fig)

    f, axes = plt.subplots(1, 2, figsize=(20, 10))
    axes[0].imshow(joint_img)
    axes[1].imshow(density_img)

    markerscale = (80 / kwargs.get('s', 20)) ** 0.5
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
    f.suptitle(joint_title, fontsize=16)
    f.subplots_adjust(top=0.9, wspace=0.05)
    f.savefig(png_file, bbox_inches='tight')
    plt.close(f)


def call_plot(x, y, log_x, log_y, hue, scatter_plot_kwargs, density_img, density_log_img):
    # logging.info(f'Joint QC plots for hue={hue}...') 
    joint_title = f'Joint QC for\n{dataset}\nmargin hue: {hue}'
       
    plot_path = output_joint / f'hue={hue}'
    plot_path.mkdir(exist_ok=True)

    # determine plotting parameters
    if is_numeric_dtype(adata.obs[hue]):
        palette = 'plasma'
        legend = 'brief'
    else:
        n_unique = adata.obs[hue].nunique()
        n_max_cat = 102
        if 10 < n_unique <= n_max_cat:
            palette = sc.pl.palettes.default_102
            palette = {cat: palette[i] for i, cat in enumerate(adata.obs[hue].unique())}
            legend = 'full'
        elif n_unique > n_max_cat:
            palette = 'turbo'
            legend = False
        else:
            palette = None
            legend = 'auto'
    
    plot_kwargs = dict(scatter_plot_kwargs)
    plot_kwargs.update(
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
        adata.obs,
        png_file=plot_path / f'{x}_vs_{y}.png',
        density_img=density_img,
        joint_title=joint_title,
        x=x,
        y=y,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds.get(x),
        y_threshold=thresholds.get(y),
        title='',
        **plot_kwargs,
    )

    # plot in log scale 
    log_x_prefix = f'log_{log_x}_' if log_x > 1 else ''
    log_y_prefix = f'log_{log_y}_' if log_y > 1 else ''
    
    create_figure(
        adata.obs,
        png_file=plot_path / f'{log_x_prefix}{x}_vs_{log_y_prefix}{y}.png',
        density_img=density_log_img,
        joint_title=joint_title,
        x=x,
        y=y,
        log_x=log_x,
        log_y=log_y,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds.get(x),
        y_threshold=thresholds.get(y),
        title='',
        **plot_kwargs,
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
    ('n_genes', 'scrublet_score', 2, 1),
    ('n_counts', 'scrublet_score', 10, 1),
]
# filter to configured metrics
coordinates = [
    c for c in coordinates if
    all(x in adata.obs.columns for x in c[:2])
]

# reduce obs to required columns only
required_columns = sorted({*hues, *[c for coords in coordinates for c in coords[:2]]})
adata.obs = adata.obs[required_columns].copy()

# subset to max of 300k cells due to high computational cost
density_data = adata.obs.sample(n=int(min(300_000, adata.n_obs)), random_state=42)
# density_data = adata.obs

def _plot_density_with_retry(df, log_x=None, log_y=None):
    density_kwargs = dict(kde_plot_kwargs)
    try:
        plot_qc_joint(
            df,
            x=x,
            y=y,
            main_plot_function=sns.kdeplot,
            log_x=log_x,
            log_y=log_y,
            x_threshold=thresholds.get(x),
            y_threshold=thresholds.get(y),
            title='',
            **density_kwargs,
        )
    except ValueError as e:
        logging.error(f'Error in plotting density: {e}')
        traceback.print_exc()

        logging.info('Retry with adjusted parameters...')
        density_kwargs.pop('bw_adjust', None)
        density_kwargs.pop('gridsize', None)
        plot_qc_joint(
            df,
            x=x,
            y=y,
            main_plot_function=sns.kdeplot,
            log_x=log_x,
            log_y=log_y,
            x_threshold=thresholds.get(x),
            y_threshold=thresholds.get(y),
            title='',
            **density_kwargs,
        )


for x, y, log_x, log_y in coordinates:
    logging.info(f'Joint QC plots for "{x}" vs "{y}"...')
    
    # temporary files
    density_png = output_joint / f'{x}_vs_{y}_density.png'
    density_log_png = output_joint / f'log_{x}_vs_{y}_density.png'
    
    logging.info('Plotting density...')
    _plot_density_with_retry(density_data)
    plt.tight_layout()
    plt.savefig(density_png, bbox_inches='tight')
    plt.close('all')

    logging.info('Plotting density for log scale...')
    _plot_density_with_retry(density_data, log_x=log_x, log_y=log_y)
    plt.tight_layout()
    plt.savefig(density_log_png, bbox_inches='tight')
    plt.close('all')

    density_img = mpimg.imread(density_png)
    density_log_img = mpimg.imread(density_log_png)

    
    def safe_call_plot(*args, **kwargs):
        try:
            return call_plot(*args, **kwargs)
        except Exception as e:
            logging.error(f"Error in plot job: {e}")
            traceback.print_exc()
            return None
    
    batch_size = max(1, len(hues) // threads)
    results = Parallel(n_jobs=threads, batch_size=batch_size, return_as='generator')(
        delayed(safe_call_plot)(
            x=x,
            y=y, 
            log_x=log_x,
            log_y=log_y,
            hue=hue,
            scatter_plot_kwargs=scatter_plot_kwargs,
            density_img=density_img,
            density_log_img=density_log_img
        ) for hue in hues
    )
    # Consume generator to trigger execution and display progress
    for _ in tqdm(results, total=len(hues), desc=f'Plotting {x} vs {y} using {threads=}, {batch_size=}'):
        pass
    
    logging.info('Removing temporary files...')
    density_png.unlink()
    density_log_png.unlink()
