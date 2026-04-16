import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd

from dataclasses import dataclass
from typing import Literal
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import seaborn as sns


# ============================================================================
# Helper
# ============================================================================

def round_values(x, prefix='', n_digits=3):
    """Format numeric values for plot annotations with fixed or scientific notation."""
    if x < 10 ** (-n_digits):
        return f'{prefix}{x:.2e}'
    elif pd.notna(x):
        return f'{prefix}{x:.{n_digits}f}'
    return ''


# ============================================================================
# PCR Plotter
# ============================================================================

PlotKind = Literal['barplot', 'violin']


@dataclass
class PCRPlotConfig:
    """Configuration for PCR score plots."""
    # Shared sizing for all plot kinds: (height_base, height_per_covariate, width_scale)
    plot_sizing: tuple[float, float, float] = (4.0, 0.45, 1.2)
    min_width: float = 4.0
    panel_fraction: float = 0.14
    min_gap_inches: float = 0.7
    label_char_width_inches: float = 0.08
    label_gap_padding_inches: float = 0.2
    main_xlabel: str = 'Principal component regression $R^2$'
    left_margin: float = 0.04
    right_margin: float = 0.02
    bottom_margin: float = 0.11
    top_margin: float = 0.93
    dpi: int = 300


class PCRPlotter:
    """Plots PCR scores as a barplot and violin plot with a left-side n_covariates panel."""

    def __init__(self, df: pd.DataFrame, config: PCRPlotConfig | None = None):
        """Prepare sorted plotting data and precomputed label strings."""
        self.config = config or PCRPlotConfig()

        df = df.sort_values(['pcr', 'n_covariates', 'covariate'], ascending=[False, True, False])
        self._df = df
        self._covariate_order = list(df['covariate'].drop_duplicates())

        stats = df.groupby('covariate', sort=False).first().reindex(self._covariate_order).copy()
        stats['pcr_string'] = stats['pcr'].apply(round_values, prefix='pcr=')
        stats['signif'] = stats['p-val'].apply(
            lambda x: '**' if x <= 0.01 else '*' if x <= 0.05 else ''
        )
        stats['z_score'] = stats['z_score'].apply(round_values, prefix='z=', n_digits=2)
        stats['p-val'] = stats['p-val'].apply(round_values, prefix='p-val=', n_digits=3)
        self._covariate_stats = stats

        self._covariate_bar_labels = stats[
            ['pcr_string', 'z_score', 'p-val', 'signif']
        ].astype(str).agg(lambda x: ', '.join([s for s in x if s]), axis=1)

        self._n_covariates_values = stats.loc[self._covariate_order, 'n_covariates'].values
        self._covariate_is_categorical = self._resolve_categorical_flags(stats)
        self._num_covariates = len(self._covariate_order)


    def _resolve_categorical_flags(self, stats: pd.DataFrame) -> pd.Series:
        """Resolve per-covariate categorical status from optional metadata columns."""
        if 'covariate_type' in stats.columns:
            return stats['covariate_type'].astype(str).str.strip().str.lower().eq('categorical')

        # Backward-compatible default: treat all covariates as categorical.
        return pd.Series(True, index=stats.index, dtype=bool)

    def _style_main_axis(self, ax_main: Axes, title: str) -> None:
        """Apply shared styling for the main plotting axis."""
        ax_main.set(title=title, ylabel='', xlabel=self.config.main_xlabel)
        ax_main.tick_params(axis='x', labelrotation=90)
        sns.despine(ax=ax_main, top=True, right=True, left=False, bottom=False)

    def _annotate_permuted_std(self, ax_main: Axes) -> None:
        """Annotate permuted bars with precomputed standard deviation labels."""
        if len(ax_main.containers) < 2:
            return
        perm_std_values = pd.to_numeric(self._covariate_stats['perm_std'], errors='coerce').fillna(0.0)
        perm_labels = perm_std_values.apply(round_values, prefix='std=')
        x_min, x_max = ax_main.get_xlim()
        x_pad = (x_max - x_min) * 0.02
        for rect, label, std in zip(ax_main.containers[1], perm_labels, perm_std_values):
            x_tip = rect.get_x() + rect.get_width()
            y_mid = rect.get_y() + rect.get_height() / 2
            ax_main.text(
                x_tip + max(float(std), 0.0) + x_pad,
                y_mid,
                label,
                va='center',
                ha='left',
                fontsize=10,
                clip_on=False,
            )

    def _add_right_summary_twin_axis(self, ax_main: Axes) -> None:
        """Attach right-side summary labels aligned to y-ticks."""
        ax_right = ax_main.twinx()
        ax_right.set_ylim(ax_main.get_ylim())
        ax_right.set_yticks(ax_main.get_yticks())
        ax_right.set_yticklabels([self._covariate_bar_labels.get(c, '') for c in self._covariate_order])
        ax_right.tick_params(left=False, right=False, length=0)
        ax_right.spines['top'].set_visible(False)
        ax_right.spines['right'].set_visible(False)

    def _create_barplot(self, title: str, ax_left: Axes, ax_main: Axes) -> None:
        """Draw the stratified PCR barplot on the main axis."""
        sns.barplot(
            data=self._df,
            x='pcr', y='covariate', hue='permuted',
            order=self._covariate_order,
            errorbar='sd', dodge=True,
            err_kws={'linewidth': 1}, capsize=0.1,
            ax=ax_main,
        )
        _, x_max = ax_main.get_xlim()
        ax_main.set_xlim(0, x_max * 1.25)
        logging.info(self._covariate_bar_labels)
        ax_main.bar_label(ax_main.containers[0], labels=self._covariate_bar_labels, padding=10)
        self._annotate_permuted_std(ax_main)

        self._draw_left_panel(ax_left, ax_main.get_ylim())
        self._style_main_axis(ax_main, title)

    def _create_violin(self, title: str, ax_left: Axes, ax_main: Axes) -> None:
        """Draw violin and strip overlays for permuted and observed PCR values."""
        sns.boxenplot(
            data=self._df[self._df['permuted']],
            x='pcr', y='covariate',
            order=self._covariate_order,
            fill=False, dodge=False,
            ax=ax_main,
        )
        sns.stripplot(
            data=self._df[~self._df['permuted']],
            x='pcr', y='covariate',
            order=self._covariate_order,
            color='red', size=8, marker='o', dodge=False,
            ax=ax_main,
        )

        self._add_right_summary_twin_axis(ax_main)

        self._draw_left_panel(ax_left, ax_main.get_ylim())
        self._style_main_axis(ax_main, title)

    def _draw_left_panel(self, ax: Axes, y_lim: tuple, xlabel: str = '# covariates') -> None:
        """Draw the left marginal panel with inverted horizontal bar labels."""
        values = pd.Series(self._n_covariates_values, dtype=float)
        is_categorical = self._covariate_is_categorical.reindex(self._covariate_order).fillna(True).to_numpy()
        max_value = max(float(values.max()), 1.0)
        y_positions = list(range(len(values)))

        categorical_y = [y for y, is_cat in zip(y_positions, is_categorical) if is_cat]
        categorical_values = [v for v, is_cat in zip(values.values, is_categorical) if is_cat]
        if categorical_values:
            max_value = max(categorical_values)
            ax.barh(categorical_y, categorical_values, color='black', alpha=0.75, height=0.8)

        ax.set_xlim(max_value * 1.2, 0)
        ax.set_ylim(y_lim)
        ax.set_xlabel(xlabel)
        ax.set_yticks([])

        for y, v, is_cat in zip(y_positions, values.values, is_categorical):
            label = str(int(v)) if float(v).is_integer() else f'{v:.2f}'
            if is_cat:
                x_pos = v + max_value * 0.03
                ha = 'right'
            else:
                # Non-categorical covariates are shown as text-only values in the margin.
                x_pos = max_value * 0.06
                ha = 'right'
            ax.text(x_pos, y, label, va='center', ha=ha, color='black', clip_on=False)

        sns.despine(ax=ax, top=True, right=False, left=True, bottom=False)

    def _estimate_tick_label_gap_inches(self) -> float:
        """Estimate required gap from longest covariate label length."""
        cfg = self.config
        longest_len = max((len(str(c)) for c in self._covariate_order), default=0)
        return max(
            cfg.min_gap_inches,
            longest_len * cfg.label_char_width_inches + cfg.label_gap_padding_inches,
        )

    def plot(self, title: str, output_path: str, kind: PlotKind) -> None:
        """Create and save a PCR plot of the requested kind."""
        cfg = self.config
        height_base, height_per_cov, width_scale = cfg.plot_sizing
        height = height_base + self._num_covariates * height_per_cov
        width = max(cfg.min_width, height * width_scale)

        # Keep spacing logic simple and predictable with a label-length estimate.
        gap_inches = self._estimate_tick_label_gap_inches()
        gap_fraction = gap_inches / width

        fig = plt.figure(figsize=(width, height), dpi=cfg.dpi)

        usable_width = 1.0 - cfg.left_margin - cfg.right_margin
        usable_height = cfg.top_margin - cfg.bottom_margin

        left_panel_width = usable_width * cfg.panel_fraction
        gap_width = min(gap_fraction, usable_width * 0.5)
        main_width = usable_width - left_panel_width - gap_width
        if main_width <= 0:
            # Fallback to a safe minimum if labels are extremely long.
            gap_width = max(0.01, usable_width * 0.2)
            main_width = usable_width - left_panel_width - gap_width

        ax_left = fig.add_axes([cfg.left_margin, cfg.bottom_margin, left_panel_width, usable_height])
        ax_main = fig.add_axes([
            cfg.left_margin + left_panel_width + gap_width,
            cfg.bottom_margin,
            main_width,
            usable_height,
        ])

        logging.info(f'{kind}...')
        if kind == 'barplot':
            self._create_barplot(title, ax_left, ax_main)
        else:
            self._create_violin(title, ax_left, ax_main)

        logging.info(f'Save to {output_path}...')
        fig.savefig(output_path, bbox_inches='tight', dpi=cfg.dpi)
        plt.close(fig)


# ============================================================================
# Snakemake script
# ============================================================================

input_file = snakemake.input.tsv
output_bar = snakemake.output.barplot
output_violin = snakemake.output.violinplot
dataset = snakemake.wildcards.dataset
file_id = snakemake.wildcards.file_id
n_permute = snakemake.params.n_permute

logging.info('Read TSV...')
df = pd.read_table(input_file)

if df.shape[0] == 0:
    logging.info('Empty TSV, skip plotting')
    for output in [output_bar, output_violin]:
        fig, ax = plt.subplots()
        ax.axis('off')
        fig.savefig(output, bbox_inches='tight', dpi=300)
        plt.close(fig)
    exit(0)

title = 'Linear variance contribution of covariates'
title += f'\ndataset={dataset}\nfile_id={file_id}\n{n_permute} permutations'

PCRPlotter(df).plot(title, output_bar, kind='barplot')
PCRPlotter(df).plot(title, output_violin, kind='violin')