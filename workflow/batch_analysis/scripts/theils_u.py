import logging
logging.basicConfig(level=logging.INFO)
from pprint import pformat
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, List, Callable
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from pandas.api.types import is_numeric_dtype
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

# ============================================================================
# Theil's U Computation
# ============================================================================

class TheilsUAnalyzer:
    """Compute and preprocess Theil's U correlation matrix for categorical variables"""
    
    def __init__(
        self,
        obs: pd.DataFrame,
        sample_key: str,
        covariates: List[str],
        max_unique_continuous: int = 10
    ):
        """
        Initialize Theil's U analyzer
        
        Args:
            obs: Observations dataframe
            sample_key: Column name(s) to use as sample identifier
            covariates: List of covariate column names to analyze
            max_unique_continuous: Skip numeric covariates with more unique values
        """
        self.obs = obs
        self.sample_key = sample_key
        self.covariates = covariates
        self.max_unique_continuous = max_unique_continuous
        
        self._valid_covariates: List[str] = []
        self._aggregated_df: Optional[pd.DataFrame] = None
        self._theils_u_matrix: Optional[pd.DataFrame] = None
    
    def prepare_sample_key(self) -> str:
        """Parse and prepare sample key column"""
        if not self.sample_key or self.sample_key == 'None':
            logging.info('Using index as sample key...')
            sample_key = 'index'
            self.obs[sample_key] = self.obs.index
            return sample_key
        
        sample_keys = [x.strip() for x in self.sample_key.split(',')]
        if len(sample_keys) > 1:
            sample_key = '--'.join(sample_keys)
            self.obs[sample_key] = (
                self.obs[sample_keys]
                .astype(str)
                .agg('-'.join, axis=1)
                .astype('category')
            )
            return sample_key
        
        return self.sample_key
    
    def filter_covariates(self) -> List[str]:
        """Filter covariates based on validity criteria"""
        logging.info(f'Filter covariates:\n{pformat(self.covariates)}')
        valid_covariates = []
        
        for covariate in self.covariates:
            if covariate not in self.obs.columns:
                logging.info(f'Skipping covariate not in obs: {covariate}')
                continue
            
            if self.obs[covariate].dropna().nunique() < 2:
                logging.info(f'Skipping covariate with fewer than 2 non-NA unique values: {covariate}')
                continue
            
            if (is_numeric_dtype(self.obs[covariate]) and 
                self.obs[covariate].nunique() > self.max_unique_continuous):
                logging.info(
                    f'Skipping numeric covariate with > {self.max_unique_continuous} '
                    f'unique values: {covariate}'
                )
                continue
            
            valid_covariates.append(covariate)
        
        self._valid_covariates = valid_covariates
        logging.info(f'Valid covariates:\n{pformat(valid_covariates)}')
        return valid_covariates
    
    def aggregate_by_sample(self, sample_key: str) -> pd.DataFrame:
        """Aggregate covariates by sample key using mode"""
        logging.info(f'Aggregate covariates by {sample_key}...')
        
        # Remove sample_key from covariates if present
        covariates = [c for c in self._valid_covariates if c != sample_key]
        assert len(covariates) > 0, 'No valid covariates found in obs'
        
        self._aggregated_df = self.obs[[sample_key] + covariates].groupby(
            sample_key, observed=True
        ).agg(
            lambda x: x.mode().iloc[0] if not x.mode().empty else np.nan
        )
        
        return self._aggregated_df
    
    def compute_theils_u(self) -> pd.DataFrame:
        """Compute pairwise Theil's U matrix"""
        logging.info(f'Compute Theil\'s U matrix...')
        
        df = self._aggregated_df
        cols = df.columns
        n = len(cols)
        
        # Precompute entropies
        entropies = {col: self._entropy(df[col]) for col in cols}
        u = np.zeros((n, n), dtype=float)
        
        for i, x in enumerate(cols):
            h_x = entropies[x]
            for j, y in enumerate(cols):
                if i == j or h_x == 0:
                    u[i, j] = 1.0
                else:
                    u[i, j] = (h_x - self._conditional_entropy(df[x], df[y])) / h_x
        
        self._theils_u_matrix = pd.DataFrame(u, index=cols, columns=cols)
        return self._theils_u_matrix
    
    @staticmethod
    def _entropy(x: pd.Series) -> float:
        """Shannon entropy H(X)"""
        p_x = x.value_counts(normalize=True)
        p_x = p_x[p_x > 0]  # filter zeros
        return -(p_x * np.log2(p_x)).sum()
    
    @staticmethod
    def _conditional_entropy(x: pd.Series, y: pd.Series) -> float:
        """Conditional entropy H(X | Y)"""
        p_xy = pd.crosstab(x, y, normalize=True)
        p_x_given_y = p_xy.div(p_xy.sum(axis=0), axis=1)
        log_p = np.where(p_x_given_y.values > 0, np.log2(p_x_given_y.values), 0.0)
        return -np.nansum(p_xy.values * log_p)
    
    def run(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Run complete analysis pipeline"""
        sample_key = self.prepare_sample_key()
        self.filter_covariates()
        self.aggregate_by_sample(sample_key)
        theils_u_matrix = self.compute_theils_u()
        
        return theils_u_matrix, self._aggregated_df


# ============================================================================
# Generalized Clustermap Plotter
# ============================================================================

@dataclass
class ClustermapConfig:
    """Configuration for clustermap visualization"""
    min_fig_size: float = 10.0
    size_per_item: float = 0.6
    min_fontsize: int = 6
    max_fontsize: int = 10
    fontsize_multiplier: float = 12.0
    bar_width: float = 0.1
    bar_gap: float = 0.02
    bar_label_threshold: float = 0.1
    dpi: int = 300
    cmap: str = "coolwarm"
    vmin: float = 0
    vmax: float = 1
    dendrogram_ratio: float = 0.1
    show_dendrogram: bool = True
    show_barplot: bool = True


class ClustermapPlotter:
    """Generalized clustermap plotter with optional dendrogram and barplot"""
    
    def __init__(
        self,
        data: pd.DataFrame,
        metadata_df: Optional[pd.DataFrame] = None,
        metadata_agg_func: Optional[Callable] = None,
        config: Optional[ClustermapConfig] = None
    ):
        """
        Initialize clustermap plotter
        
        Args:
            data: Data matrix to plot (rows and columns will be clustered)
            metadata_df: Optional dataframe for computing barplot metadata
            metadata_agg_func: Function to aggregate metadata (e.g., lambda col: col.nunique())
            config: Plot configuration
        """
        self.data = data
        self.metadata_df = metadata_df if metadata_df is not None else data
        self.metadata_agg_func = metadata_agg_func or (lambda col: col.nunique())
        self.config = config or ClustermapConfig()
        
        self._fig: Optional[Figure] = None
        self._heatmap_ax: Optional[Axes] = None
        self._clustermap = None  # Store clustermap object for barplot access
        self._reordered_labels: Optional[pd.Index] = None
        self._fig_size: float = 0
        self._fontsize: int = 0
    
    def plot(
        self,
        title: str,
        output_path: str,
        cbar_label: str = "Value",
        bar_xlabel: str = "# unique"
    ) -> None:
        """
        Create and save the complete plot
        
        Args:
            title: Plot title
            output_path: Path to save figure
            cbar_label: Label for colorbar
            bar_xlabel: Label for barplot x-axis
        """
        self._calculate_plot_parameters()
        
        if self.config.show_dendrogram and self.data.shape[0] > 2:
            self._create_clustermap(cbar_label)
        else:
            self._create_simple_heatmap(cbar_label)
        
        if self.config.show_barplot:
            self._add_barplot(bar_xlabel)
        
        self._configure_heatmap_axis()
        
        # Use tight_layout to automatically adjust spacing, reserving top space for title
        self._fig.suptitle(title, y=0.98, fontsize=12)
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.config.dpi)
        plt.close(self._fig)
    
    def _calculate_plot_parameters(self) -> None:
        """Calculate figure size and fontsize"""
        n_items = len(self.data)
        self._fig_size = max(self.config.min_fig_size, n_items * self.config.size_per_item)
        self._fontsize = max(
            self.config.min_fontsize,
            min(self.config.max_fontsize, int((self._fig_size / n_items) * self.config.fontsize_multiplier))
        )
    
    def _create_clustermap(self, cbar_label: str) -> None:
        """Create clustered heatmap with dendrogram"""
        # Convert similarity to distance (diagonal is already 0 since self-similarity = 1)
        distances = squareform(1 - self.data.values, checks=False)
        col_linkage = linkage(distances, method='average')

        self._clustermap = sns.clustermap(
            self.data,
            col_linkage=col_linkage,
            row_linkage=col_linkage,  # Use same linkage for rows
            col_cluster=True,
            row_cluster=True,  # Enable row clustering for axis
            annot=True,
            fmt=".2f",
            cmap=self.config.cmap,
            vmin=self.config.vmin,
            vmax=self.config.vmax,
            figsize=(self._fig_size, self._fig_size),
            cbar_kws={"label": cbar_label},
            cbar_pos=(0.02, 0.8, 0.03, 0.15),  # Position colorbar top-left: (left, bottom, width, height)
            dendrogram_ratio=self.config.dendrogram_ratio,
            xticklabels=True,
            yticklabels=False,
            annot_kws={"size": self._fontsize},
        )
        
        self._fig = self._clustermap.fig
        self._heatmap_ax = self._clustermap.ax_heatmap
        # Get reordered labels from row dendrogram for consistent ordering
        self._reordered_labels = self.data.index[self._clustermap.dendrogram_row.reordered_ind]
    
    def _create_simple_heatmap(self, cbar_label: str) -> None:
        """Create simple heatmap without clustering"""
        self._fig, self._heatmap_ax = plt.subplots(figsize=(self._fig_size, self._fig_size))
        
        sns.heatmap(
            self.data,
            ax=self._heatmap_ax,
            annot=True,
            fmt=".2f",
            cmap=self.config.cmap,
            vmin=self.config.vmin,
            vmax=self.config.vmax,
            cbar_kws={"label": cbar_label},
            xticklabels=True,
            yticklabels=False,
            annot_kws={"size": self._fontsize},
        )
        
        self._reordered_labels = self.data.index
    
    def _add_barplot(self, bar_xlabel: str) -> None:
        """Add barplot in place of row dendrogram"""
        # Clear the row dendrogram axis and replace with barplot
        bar_ax = self._clustermap.ax_row_dendrogram
        bar_ax.clear()
        
        # Compute metadata values in reordered labels order
        metadata_values = [
            self.metadata_agg_func(self.metadata_df[col]) 
            for col in self._reordered_labels
        ]
        
        # Create horizontal barplot
        bars = bar_ax.barh(range(len(metadata_values)), metadata_values, color='black', alpha=0.7)
        
        # Match heatmap alignment
        bar_ax.set_ylim(self._heatmap_ax.get_ylim())
        bar_ax.invert_xaxis()  # Match dendrogram orientation
        bar_ax.set_xlabel(bar_xlabel, fontsize=self._fontsize)
        
        # Add value labels
        threshold = max(metadata_values) * self.config.bar_label_threshold
        for i, (bar, v) in enumerate(zip(bars, metadata_values)):
            label = f'n={v}'
            
            if v > threshold:
                # Inside bar (white text)
                bar_ax.text(
                    v * 0.5, i, label,
                    va='center', ha='center',
                    fontsize=self._fontsize,
                    color='white'
                )
            else:
                # Outside bar (black text)
                bar_ax.text(
                    max(metadata_values) * -0.05, i, label,
                    va='center', ha='right',
                    fontsize=self._fontsize,
                    color='black'
                )
        
        sns.despine(ax=bar_ax, left=True, bottom=False)
    
    def _configure_heatmap_axis(self) -> None:
        """Configure heatmap axis labels"""
        self._heatmap_ax.set_ylabel('')
        self._heatmap_ax.set_xticklabels(
            self._heatmap_ax.get_xticklabels(),
            rotation=90,
            ha='right',
            fontsize=self._fontsize
        )


# ============================================================================
# Snakemake Script
# ============================================================================
from utils.io import read_anndata


# Extract snakemake parameters
input_file = snakemake.input[0]
output_png = snakemake.output.plot

sample_key = snakemake.params.get('sample_key')
covariates = snakemake.params.get('covariates', [])
max_unique_continuous = snakemake.params.get('max_unique_continuous', 10)
na_strings = snakemake.params.get('na_strings', ['NA', 'NaN', 'nan', ''])
wildcards = snakemake.wildcards

title = f"Theil's U Heatmap\ndataset: {wildcards['dataset']}, file_id: {wildcards['file_id']}"

# Read data
logging.info(f'Read {input_file}...')
obs = read_anndata(input_file, obs='obs', verbose=False).obs

# Run Theil's U analysis
analyzer = TheilsUAnalyzer(
    obs=obs,
    sample_key=sample_key,
    covariates=covariates,
    max_unique_continuous=max_unique_continuous
)
theils_u_matrix, aggregated_df = analyzer.run()

# Plot
# logging.info(f'Plot Theil\'s U heatmap to {output_png}...')
# plotter = ClustermapPlotter(
#     data=theils_u_matrix,
#     metadata_df=aggregated_df,
#     metadata_agg_func=lambda col: col.nunique()
# )
# plotter.plot(
#     title=title,
#     output_path=output_png,
#     cbar_label="Theil's U",
#     bar_xlabel="# unique"
# )
# logging.info('Done.')



### old code
logging.info(f'Plot Theil\'s U heatmap to {output_png} ...')
width = max(10, len(theils_u_matrix) * 0.6)
height = max(10, len(theils_u_matrix) * 0.6)
n_cells = max(1, theils_u_matrix.shape[0])
cell_h_inches = height / n_cells
# scale annotation font size by cell height (inches → points)
annot_fontsize = max(6, min(10, cell_h_inches * 12))

# add unique value counts to labels
covariates = theils_u_matrix.columns
n_unique = {col: aggregated_df[col].nunique() for col in covariates}
labels_with_counts = [f"{col} (n={n_unique[col]})" for col in covariates]

# rename index and columns with counts
theils_u_matrix.index = labels_with_counts
theils_u_matrix.columns = labels_with_counts

heatmap_kwargs = dict(
    annot=True,
    fmt='.2f',
    cmap="coolwarm",
    vmin=0,
    vmax=1,
    cbar_kws={"label": "Theil's U"},
    xticklabels=True,
    yticklabels=True,
    annot_kws={"size": annot_fontsize},
)

title = f"Theil's U Heatmap of covariates\ndataset: {wildcards['dataset']}\nfile_id: {wildcards['file_id']}"

if theils_u_matrix.shape[0] > 2:
    import scipy
    
    dist = scipy.spatial.distance.squareform(1 - theils_u_matrix.values, checks=False)
    link = scipy.cluster.hierarchy.linkage(dist, method="average")

    g = sns.clustermap(
        theils_u_matrix,
        row_linkage=link,
        col_linkage=link,   # same clustering → same ordering → diagonal preserved
        **heatmap_kwargs,
        dendrogram_ratio=0.1,
        figsize=(width, height),
    )
    g.fig.suptitle(title, y=0.98)
    plt.tight_layout()
else:
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        theils_u_matrix,
        ax=ax,
        **heatmap_kwargs,
    )
    ax.set_title(title, pad=20)

plt.tight_layout()
plt.savefig(output_png, bbox_inches='tight')
