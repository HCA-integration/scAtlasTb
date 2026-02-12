import ast
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
            if 1 < adata.obs[hue].nunique() < max_groups
        }
        hues = list(hues)
    assert len(hues) > 0, 'No valid hue columns provided for plotting'

    return dataset, hues


def parse_autoqc(autoqc_thresholds: pd.DataFrame, thresholds: dict = None):
    if thresholds is None:
        thresholds = {}
    thresholds |= {
        f'{key}_min': autoqc_thresholds.loc[key, 'low']
        for key in autoqc_thresholds.index
    }
    thresholds |={
        f'{key}_max': autoqc_thresholds.loc[key, 'high']
        for key in autoqc_thresholds.index
    }
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
        adata.obs[column_name] = adata.obs[column_name] \
            & adata.obs[key].between(*thresholds[key])


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
    title='',
    return_df=False,
    marginal_kwargs: dict=None,
    **kwargs,
):
    """
    Plot scatter plot with marginal histograms from df columns.

    :param df: observation dataframe
    :param x: df column for x axis
    :param y: df column for y axis
    :param log: log base for transforming values. Default 1, no transformation
    :param hue: df column with annotations for color coding scatter plot points
    :param marginal_hue: df column with annotations for color coding marginal plot distributions
    :param x_threshold: tuple of upper and lower filter thresholds for x axis
    :param y_threshold: tuple of upper and lower filter thresholds for y axis
    :param title: Title text for plot
    :return:
        seaborn plot (and df dataframe with updated values, if `return_df=True`)
    """
    import seaborn as sns
    
    if main_plot_function is None:
        main_plot_function = sns.scatterplot
    if not x_threshold:
        x_threshold=(0, np.inf)
    if not y_threshold:
        y_threshold=(0, np.inf)

    log_x = 1 if log_x is None else log_x
    log_y = 1 if log_y is None else log_y

    def log1p_base(_x, base):
        return np.log1p(_x) / np.log(base)

    if log_x > 1:
        x_log = f'log{log_x} {x}'
        df[x_log] = log1p_base(df[x], log_x)
        x_threshold = log1p_base(x_threshold, log_x)
        x = x_log
    
    if log_y > 1:
        y_log = f'log{log_y} {y}'
        df[y_log] = log1p_base(df[y], log_y)
        y_threshold = log1p_base(y_threshold, log_y)
        y = y_log
        
    if marginal_kwargs is None:
        marginal_kwargs = dict(legend=False)
    
    if marginal_hue in df.columns:
        marginal_hue = None if df[marginal_hue].nunique() > 100 else marginal_hue
    use_marg_hue = marginal_hue is not None
    
    if not use_marg_hue:
         marginal_kwargs.pop('palette', None)
    
    g = sns.JointGrid(
        data=df,
        x=x,
        y=y,
        xlim=(0, df[x].max()),
        ylim=(0, df[y].max()),
    )
    
    # main plot
    if hue is not None and hue in df.columns and not pd.api.types.is_numeric_dtype(df[hue]):
        hue_counts = df[hue].value_counts(dropna=False)
        # common -> rare (rarest last)
        df = df.assign(_hue_freq=df[hue].map(hue_counts)).sort_values('_hue_freq', ascending=False)
    
    g.plot_joint(main_plot_function, data=df, hue=hue, **kwargs)
    
    # marginal hist plot
    g.plot_marginals(
        sns.histplot,
        data=df,
        hue=marginal_hue,
        element='step' if use_marg_hue else 'bars',
        fill=False,
        bins=100,
        **marginal_kwargs,
    )

    g.fig.suptitle(title, fontsize=12)
    # workaround for patchworklib
    g._figsize = g.fig.get_size_inches()

    if hue is not None:
        # handles, labels = g.ax_joint.get_legend_handles_labels()
        markerscale = (80 / kwargs.get('s', 20)) ** 0.5
        g.ax_joint.legend(markerscale=markerscale)

    # x threshold
    for t, t_def in zip(x_threshold, (0, np.inf)):
        if t != t_def:
            g.ax_joint.axvline(x=t, color='red')
            g.ax_marg_x.axvline(x=t, color='red')

    # y threshold
    for t, t_def in zip(y_threshold, (0, np.inf)):
        if t != t_def:
            g.ax_joint.axhline(y=t, color='red')
            g.ax_marg_y.axhline(y=t, color='red')

    if return_df:
        return g, df
    return g
