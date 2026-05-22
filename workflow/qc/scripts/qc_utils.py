import ast
from types import SimpleNamespace

import numpy as np
import pandas as pd
import anndata as ad


QC_FLAGS = [
    'n_counts',
    'n_genes',
    'percent_mito',
    'n_counts_mito',
    'percent_ribo',
    'n_counts_ribo',
    'percent_hb',
    'n_counts_hb',
]


def parse_parameters(adata: ad.AnnData, params: dict, filter_hues: bool = False):
    dataset = params.get('dataset', 'None')
    hues = params.get('hue', [])
    
    split_datasets = dataset.split('--')
    if len(split_datasets) > 1:
        dataset = ' '.join([split_datasets[0], split_datasets[-1]])

    if isinstance(hues, str):
        hues = [hues]
    hues = [hue for hue in hues if hue in adata.obs.columns]
    if filter_hues:
        max_groups = params.get('max_groups', 100)
        hues = {
            hue for hue in hues
            if (
                # Keep numeric columns regardless of unique value count (for continuous colormaps)
                pd.api.types.is_numeric_dtype(adata.obs[hue])
                # For categorical/non-numeric, apply max_groups filter
                or (1 < adata.obs[hue].nunique() < max_groups)
            )
        }
        hues = list(hues)
    if len(hues) == 0:
        hues = [None]

    return dataset, hues


def parse_autoqc(autoqc_thresholds: pd.DataFrame, thresholds: dict = None):
    if thresholds is None:
        thresholds = {}
    if autoqc_thresholds is None or autoqc_thresholds.empty:
        return thresholds

    has_side = 'side' in autoqc_thresholds.columns
    for key in autoqc_thresholds.index:
        low = autoqc_thresholds.loc[key, 'low']
        high = autoqc_thresholds.loc[key, 'high']
        side = autoqc_thresholds.loc[key, 'side'] if has_side else None

        low = None if pd.isna(low) else low
        high = None if pd.isna(high) else high

        if isinstance(side, str):
            side = side.lower()
            if side == 'max_only':
                low = None
            elif side == 'min_only':
                high = None

        thresholds[f'{key}_min'] = low
        thresholds[f'{key}_max'] = high

    return thresholds


def get_thresholds(
    threshold_keys: list = None,
    autoqc_thresholds: pd.DataFrame = None,
    user_thresholds: [str, dict] = None,
    transform: bool = True,
    init_nan: bool = False,
):
    """
    :param threshold_keys: keys in mappings that define different QC paramters
    :param autoqc_thresholds: autoQC thresholds as provided by sctk under adata.uns['scautoqc_ranges']
    :param user_thresholds: user defined thresholds
    :return: thresholds: dict of key: thresholds tuple
    """
    import ast
    
    # set default thresholds
    if threshold_keys is None:
        threshold_keys = ['n_counts', 'n_genes', 'percent_mito']
    
    # initialise thresholds
    if init_nan:
        thresholds = {f'{key}_min': np.nan for key in threshold_keys}
        thresholds |= {f'{key}_max': np.nan for key in threshold_keys}
    else:
        thresholds = {f'{key}_min': 0 for key in threshold_keys}
        thresholds |= {f'{key}_max': np.inf for key in threshold_keys}

    if autoqc_thresholds is not None:
        thresholds = parse_autoqc(autoqc_thresholds, thresholds)

    # update to user thresholds
    if user_thresholds is None:
        user_thresholds = {}
    elif isinstance(user_thresholds, str):
        user_thresholds = ast.literal_eval(user_thresholds)
    elif isinstance(user_thresholds, dict):
        pass
    else:
        ValueError('thresholds must be a dict or string')
    thresholds |= user_thresholds
    
    if transform:
        # transform to shape expected by plot_qc_joint
        return {
            key: (thresholds[f'{key}_min'], thresholds[f'{key}_max'])
            for key in threshold_keys
        }
    return {
        key: value
        for key, value in thresholds.items()
        if any([key.startswith(x) for x in threshold_keys])
    }


def apply_thresholds(
    adata: ad.AnnData,
    thresholds: dict,
    threshold_keys: list,
    column_name='passed_qc'
):
    """
    :param adata: AnnData object
    :param thresholds: dict of key: thresholds tuple as returned by get_thresholds
    """
    adata.obs[column_name] = True
    if adata.n_obs == 0:
        return
    for key in threshold_keys:
        lower, upper = thresholds[key]
        if lower is None and upper is None:
            continue
        if lower is None:
            passed = adata.obs[key] <= upper
        elif upper is None:
            passed = adata.obs[key] >= lower
        else:
            passed = adata.obs[key].between(lower, upper)
        adata.obs[column_name] = adata.obs[column_name] & passed.fillna(False)


def plot_qc_joint(
    df: pd.DataFrame,
    x: str,
    y: str,
    log_x: int = 1,
    log_y: int = 1,
    hue: str = None,
    main_plot_function=None,
    marginal_hue=None,
    x_threshold=None,
    y_threshold=None,
    x_threshold2=None,
    y_threshold2=None,
    threshold_color='black',
    threshold_color2='black',
    threshold_linestyle='-',
    threshold_linestyle2='--',
    title='',
    return_df=False,
    marginal_kwargs: dict = None,
    fig=None,
    subplot_spec=None,
    sharey=None,
    **kwargs,
):
    """
    Plot scatter plot with marginal histograms from df columns.

    :param df: observation dataframe
    :param x: df column for x axis
    :param y: df column for y axis
    :param log_x: log base for transforming x values. Default 1, no transformation
    :param log_y: log base for transforming y values. Default 1, no transformation
    :param hue: df column with annotations for color coding scatter plot points
    :param marginal_hue: df column with annotations for color coding marginal plot distributions
    :param x_threshold: tuple of (min, max) filter thresholds for x axis
    :param y_threshold: tuple of (min, max) filter thresholds for y axis
    :param fig: existing matplotlib Figure used when `subplot_spec` is provided
    :param subplot_spec: optional matplotlib SubplotSpec to draw into an existing grid layout.
        When provided, creates joint/marginal axes inside this slot instead of creating a new JointGrid figure.
    :param sharey: existing Axes to share y axis with (subplot_spec mode only)
    :param title: title text for plot
    :return: seaborn plot (and df dataframe with updated values, if `return_df=True`)
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    df = df.copy()  # prevent in-place mutation of caller's dataframe

    if main_plot_function is None:
        main_plot_function = sns.scatterplot

    def _normalize_threshold_pair(threshold):
        # None/NaN bounds mean open interval on that side.
        if threshold is None:
            return (0, np.inf)
        lower, upper = threshold
        lower = 0 if lower is None or pd.isna(lower) else lower
        upper = np.inf if upper is None or pd.isna(upper) else upper
        return (lower, upper)

    x_threshold = _normalize_threshold_pair(x_threshold)
    y_threshold = _normalize_threshold_pair(y_threshold)
    x_threshold2 = _normalize_threshold_pair(x_threshold2)
    y_threshold2 = _normalize_threshold_pair(y_threshold2)

    def log1p_base(_x, base):
        return np.log1p(_x) / np.log(base)

    if log_x > 1:
        x_log = f'log{log_x} {x}'
        df[x_log] = log1p_base(df[x], log_x)
        x_threshold  = log1p_base(x_threshold,  log_x)
        x_threshold2 = log1p_base(x_threshold2, log_x)
        x = x_log

    if log_y > 1:
        y_log = f'log{log_y} {y}'
        df[y_log] = log1p_base(df[y], log_y)
        y_threshold  = log1p_base(y_threshold,  log_y)
        y_threshold2 = log1p_base(y_threshold2, log_y)
        y = y_log

    def thresholds_equal(a, b):
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        return a.shape == b.shape and np.allclose(a, b, equal_nan=True)

    # Avoid plotting duplicate autoQC lines when they match user thresholds
    if thresholds_equal(x_threshold, x_threshold2):
        x_threshold2 = (0, np.inf)
    if thresholds_equal(y_threshold, y_threshold2):
        y_threshold2 = (0, np.inf)

    if marginal_kwargs is None:
        marginal_kwargs = dict(legend=False)

    if marginal_hue in df.columns:
        marginal_hue = None if df[marginal_hue].nunique() > 100 else marginal_hue
    use_marg_hue = marginal_hue is not None

    if not use_marg_hue:
        marginal_kwargs.pop('palette', None)

    if hue in df.columns and pd.api.types.is_categorical_dtype(df[hue]):
        # sort so smaller groups are plotted on top of larger groups
        hue_order = df[hue].value_counts(ascending=False).index
    else:
        hue_order = None

    def _draw_thresholds(ax_joint, ax_marg_x, ax_marg_y):
        # Plot autoQC (2) first, then user/updated (1) on top
        for t, t_def in zip(x_threshold2, (0, np.inf)):
            if t != t_def:
                ax_joint.axvline(x=t, color=threshold_color2, linestyle=threshold_linestyle2)
                ax_marg_x.axvline(x=t, color=threshold_color2, linestyle=threshold_linestyle2)
        for t, t_def in zip(x_threshold, (0, np.inf)):
            if t != t_def:
                ax_joint.axvline(x=t, color=threshold_color, linestyle=threshold_linestyle)
                ax_marg_x.axvline(x=t, color=threshold_color, linestyle=threshold_linestyle)
        for t, t_def in zip(y_threshold2, (0, np.inf)):
            if t != t_def:
                ax_joint.axhline(y=t, color=threshold_color2, linestyle=threshold_linestyle2)
                ax_marg_y.axhline(y=t, color=threshold_color2, linestyle=threshold_linestyle2)
        for t, t_def in zip(y_threshold, (0, np.inf)):
            if t != t_def:
                ax_joint.axhline(y=t, color=threshold_color, linestyle=threshold_linestyle)
                ax_marg_y.axhline(y=t, color=threshold_color, linestyle=threshold_linestyle)

    def _plot_on_axes(ax_joint, ax_marg_x, ax_marg_y):
        main_plot_function(
            data=df, x=x, y=y, hue=hue, hue_order=hue_order, ax=ax_joint, **kwargs,
        )

        if hue is not None and kwargs.get('legend', True):
            handles, labels = ax_joint.get_legend_handles_labels()
            if handles:
                ax_joint.legend(
                    handles=handles,
                    labels=labels,
                    markerscale=(60 / kwargs.get('s', 20)) ** 0.5,
                    fontsize=ax_joint.xaxis.label.get_size(),
                )

        hist_kwargs = dict(
            data=df,
            hue=marginal_hue,
            hue_order=hue_order if use_marg_hue else None,
            element='step' if use_marg_hue else 'bars',
            fill=False,
            bins=100,
            **marginal_kwargs,
        )
        sns.histplot(x=x, ax=ax_marg_x, **hist_kwargs)
        sns.histplot(y=y, ax=ax_marg_y, **hist_kwargs)
        _draw_thresholds(ax_joint, ax_marg_x, ax_marg_y)

        ax_marg_x.tick_params(axis='x', bottom=False, labelbottom=False)
        ax_marg_y.tick_params(axis='y', left=False, labelleft=False)
        ax_marg_x.spines['top'].set_visible(False)
        ax_marg_x.spines['right'].set_visible(False)
        ax_marg_y.spines['top'].set_visible(False)
        ax_marg_y.spines['right'].set_visible(False)
        for ax in (ax_marg_x, ax_marg_y):
            ax.set_xlabel('')
            ax.set_ylabel('')

        x_max = float(np.nanmax(df[x]))
        y_max = float(np.nanmax(df[y]))
        if x_max >= 1e4:
            ax_joint.ticklabel_format(axis='x', style='sci', scilimits=(4, 4))
            ax_marg_x.ticklabel_format(axis='x', style='sci', scilimits=(4, 4))
        if y_max >= 1e4:
            ax_joint.ticklabel_format(axis='y', style='sci', scilimits=(4, 4))
            ax_marg_y.ticklabel_format(axis='y', style='sci', scilimits=(4, 4))

        ax_joint.set_xlim(0, x_max)
        ax_joint.set_ylim(0, y_max)
        ax_joint.set_xlabel(x)
        ax_joint.set_ylabel(y)
        ax_joint.spines['top'].set_visible(False)
        ax_joint.spines['right'].set_visible(False)

        return SimpleNamespace(
            fig=ax_joint.figure,
            ax_joint=ax_joint,
            ax_marg_x=ax_marg_x,
            ax_marg_y=ax_marg_y,
            _figsize=ax_joint.figure.get_size_inches(),
        )

    if subplot_spec is not None:
        if fig is None:
            fig = getattr(subplot_spec, 'figure', None)
            if fig is None:
                fig = subplot_spec.get_gridspec().figure
        inner = subplot_spec.subgridspec(
            2, 2,
            height_ratios=[1, 4],
            width_ratios=[4, 1],
            hspace=0,
            wspace=0,
        )
        ax_joint  = fig.add_subplot(inner[1, 0], sharey=sharey)
        ax_marg_x = fig.add_subplot(inner[0, 0], sharex=ax_joint)
        ax_marg_y = fig.add_subplot(inner[1, 1], sharey=ax_joint)
        fig.add_subplot(inner[0, 1]).set_axis_off()

        # Hide tick labels at the shared boundary
        plt.setp(ax_marg_x.get_xticklabels(), visible=False)
        plt.setp(ax_marg_y.get_yticklabels(), visible=False)
        ax_marg_x.tick_params(axis='x', bottom=False)
        ax_marg_y.tick_params(axis='y', left=False)

        # When sharing y with an external axis, hide the y tick labels to avoid duplication
        if sharey is not None:
            plt.setp(ax_joint.get_yticklabels(), visible=False)
            ax_joint.set_ylabel('')

        # Remove facing spines so marginal and joint look connected
        ax_marg_x.spines['bottom'].set_visible(False)
        ax_marg_y.spines['left'].set_visible(False)

        g = _plot_on_axes(ax_joint, ax_marg_x, ax_marg_y)
        g.fig.suptitle(title, fontsize=14)
        if return_df:
            return g, df
        return g

    # JointGrid handles sharex/sharey internally
    g = sns.JointGrid(data=df, x=x, y=y, xlim=(0, df[x].max()), ylim=(0, df[y].max()))
    g = _plot_on_axes(g.ax_joint, g.ax_marg_x, g.ax_marg_y)
    g.fig.suptitle(title, fontsize=14)
    g._figsize = g.fig.get_size_inches()
    if return_df:
        return g, df
    return g
