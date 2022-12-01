import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import anndata as ad
import pandas as pd
import seaborn as sns
import scipy.sparse

import warnings
import collections.abc as cabc
from pathlib import Path
from types import MappingProxyType
from typing import Optional, Union, List, Sequence, Mapping, Any, Tuple

import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from pandas.api.types import is_categorical_dtype
from matplotlib import pyplot as pl, rcParams, ticker
from matplotlib import patheffects
from matplotlib.axes import Axes
from matplotlib.colors import is_color_like, Colormap
from scipy.sparse import issparse
from sklearn.utils import check_random_state

from scanpy.plotting import _utils
from scanpy.plotting._utils import matrix, _IGraphLayout, _FontWeight, _FontSize
from scanpy import _utils as _sc_utils, logging as logg
from scanpy._settings import settings
from scanpy._compat import Literal

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' 
        from: https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
        creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def get_act(adata, imputed=False, fill_value=np.nan):
    if imputed:
        key='activities_imputed'
    else:
        key='activities'
    adata = ad.AnnData(X=adata.obsm[key], obs=adata.obs, uns=adata.uns, obsm=adata.obsm, var=pd.DataFrame(adata.uns['gene_names'], index=adata.uns['gene_names']))
    adata[adata.obs['atac_match']==0].X = fill_value
    return adata

def get_chromvar(adata, fill_value=0):
    adata_chromvar = ad.AnnData(X=adata.obsm['chromvar'], obs=adata.obs, uns=adata.uns, obsm=adata.obsm, var=pd.DataFrame(adata.uns['motif_names'], index=adata.uns['motif_names']))

    #map motifs
    from pyjaspar import jaspardb
    jdb_obj = jaspardb(release='JASPAR2020')
    motifs = jdb_obj.fetch_motifs(
        collection = 'CORE',
        tax_group = ['vertebrates']
        )

    motif_names = np.array([np.array([motif.name, motif.matrix_id]) for motif in motifs])

    motif_map = dict(zip(motif_names[:,1], motif_names[:,0]))

    adata_chromvar.var.index = adata_chromvar.var.index.map(motif_map)
    adata_chromvar[ adata_chromvar.obs['atac_match']==0].X = fill_value
    
    return adata_chromvar


def get_row_colors(data, color='red'):
    #create a color palette with the same number of colors as unique values in the Source column
    network_pal = sns.color_palette(color, as_cmap=False, n_colors=len(data.unique()))

    #Create a dictionary where the key is the category and the values are the
    #colors from the palette we just created
    network_lut = dict(zip(data.unique(), network_pal))

    #get the series of all of the categories
    networks = data

    #map the colors to the series. Now we have a list of colors the same
    #length as our dataframe, where unique values are mapped to the same color
    network_colors = pd.Series(networks).map(network_lut)
    return network_colors

def plot_heatmap(adata, groupby, markers, layer=None, **kwargs):
    test = adata[adata.obs.sort_values(groupby).index]
    
    if layer:
        X = test[:, markers].layers[layer].A.T
    else:
        if scipy.sparse.issparse(test[:, markers].X):
            X = test[:, markers].X.A.T
        else:
            X = test[:, markers].X.T
        
    sns_df = pd.DataFrame(X,  index=markers, columns=test.obs_names)

    row_colors = get_row_colors(test.obs[groupby[0]], "crest").to_frame().join( get_row_colors(test.obs[groupby[1]], "flare"))
    
    row_colors.columns = [f'Absorption probability {groupby[0]}', "Pseudotime"]
    g = sns.clustermap(sns_df, 
                   col_cluster=False, 
                   standard_scale=0, 
                   figsize=(17,11), 
                   cmap='Purples', 
                   yticklabels=True,
                   dendrogram_ratio=0.1,
                    col_colors=row_colors,
                    **kwargs
                  )
    g.ax_heatmap.set(xticklabels=[], xticks=[]);
    g.ax_row_dendrogram.set_visible(False)
    return g


def paga_path(
    adata: AnnData,
    nodes: Sequence[Union[str, int]],
    keys: Sequence[str],
    use_raw: bool = True,
    annotations: Sequence[str] = ('dpt_pseudotime',),
    color_map: Union[str, Colormap, None] = None,
    color_maps_annotations: Mapping[str, Union[str, Colormap]] = MappingProxyType(
        dict(dpt_pseudotime='Greys')
    ),
    palette_groups: Optional[Sequence[str]] = None,
    n_avg: int = 1,
    groups_key: Optional[str] = None,
    xlim: Tuple[Optional[int], Optional[int]] = (None, None),
    title: Optional[str] = None,
    left_margin=None,
    ytick_fontsize: Optional[int] = None,
    title_fontsize: Optional[int] = None,
    show_node_names: bool = True,
    show_yticks: bool = True,
    show_colorbar: bool = True,
    legend_fontsize: Union[int, float, _FontSize, None] = None,
    legend_fontweight: Union[int, _FontWeight, None] = None,
    normalize_to_zero_one: bool = False,
    as_heatmap: bool = True,
    return_data: bool = False,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    ax: Optional[Axes] = None,
) -> Optional[Axes]:
    """\
    Gene expression and annotation changes along paths in the abstracted graph.
    Parameters
    ----------
    adata
        An annotated data matrix.
    nodes
        A path through nodes of the abstracted graph, that is, names or indices
        (within `.categories`) of groups that have been used to run PAGA.
    keys
        Either variables in `adata.var_names` or annotations in
        `adata.obs`. They are plotted using `color_map`.
    use_raw
        Use `adata.raw` for retrieving gene expressions if it has been set.
    annotations
        Plot these keys with `color_maps_annotations`. Need to be keys for
        `adata.obs`.
    color_map
        Matplotlib colormap.
    color_maps_annotations
        Color maps for plotting the annotations. Keys of the dictionary must
        appear in `annotations`.
    palette_groups
        Ususally, use the same `sc.pl.palettes...` as used for coloring the
        abstracted graph.
    n_avg
        Number of data points to include in computation of running average.
    groups_key
        Key of the grouping used to run PAGA. If `None`, defaults to
        `adata.uns['paga']['groups']`.
    as_heatmap
        Plot the timeseries as heatmap. If not plotting as heatmap,
        `annotations` have no effect.
    show_node_names
        Plot the node names on the nodes bar.
    show_colorbar
        Show the colorbar.
    show_yticks
        Show the y ticks.
    normalize_to_zero_one
        Shift and scale the running average to [0, 1] per gene.
    return_data
        Return the timeseries data in addition to the axes if `True`.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on \\{`'.pdf'`, `'.png'`, `'.svg'`\\}.
    ax
         A matplotlib axes object.
    Returns
    -------
    A :class:`~matplotlib.axes.Axes` object, if `ax` is `None`, else `None`.
    If `return_data`, return the timeseries data in addition to an axes.
    """
    ax_was_none = ax is None

    if groups_key is None:
        if 'groups' not in adata.uns['paga']:
            raise KeyError(
                'Pass the key of the grouping with which you ran PAGA, '
                'using the parameter `groups_key`.'
            )
        groups_key = adata.uns['paga']['groups']
    groups_names = adata.obs[groups_key].cat.categories

    if 'dpt_pseudotime' not in adata.obs.keys():
        raise ValueError(
            '`pl.paga_path` requires computation of a pseudotime `tl.dpt` '
            'for ordering at single-cell resolution'
        )

    if palette_groups is None:
        _utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
        palette_groups = adata.uns[f'{groups_key}_colors']

    def moving_average(a):
        return _sc_utils.moving_average(a, n_avg)

    ax = pl.gca() if ax is None else ax

    X = []
    x_tick_locs = [0]
    x_tick_labels = []
    groups = []
    anno_dict = {anno: [] for anno in annotations}
    if isinstance(nodes[0], str):
        nodes_ints = []
        groups_names_set = set(groups_names)
        for node in nodes:
            if node not in groups_names_set:
                raise ValueError(
                    f'Each node/group needs to be in {groups_names.tolist()} '
                    f'(`groups_key`={groups_key!r}) not {node!r}.'
                )
            nodes_ints.append(groups_names.get_loc(node))
        nodes_strs = nodes
    else:
        nodes_ints = nodes
        nodes_strs = [groups_names[node] for node in nodes]

    adata_X = adata
    if use_raw and adata.raw is not None:
        adata_X = adata.raw

    for ikey, key in enumerate(keys):
        x = []
        for igroup, group in enumerate(nodes_ints):
            idcs = np.arange(adata.n_obs)[
                adata.obs[groups_key].values == nodes_strs[igroup]
            ]
            if len(idcs) == 0:
                raise ValueError(
                    'Did not find data points that match '
                    f'`adata.obs[{groups_key!r}].values == {str(group)!r}`. '
                    f'Check whether `adata.obs[{groups_key!r}]` '
                    'actually contains what you expect.'
                )
            idcs_group = np.argsort(
                adata.obs['dpt_pseudotime'].values[
                    adata.obs[groups_key].values == nodes_strs[igroup]
                ]
            )
            idcs = idcs[idcs_group]
            values = (
                adata.obs[key].values if key in adata.obs_keys() else adata_X[:, key].X
            )[idcs]
            x += (values.A if issparse(values) else values).tolist()
            if ikey == 0:
                groups += [group] * len(idcs)
                x_tick_locs.append(len(x))
                for anno in annotations:
                    series = adata.obs[anno]
                    if is_categorical_dtype(series):
                        series = series.cat.codes
                    anno_dict[anno] += list(series.values[idcs])
        if n_avg > 1:
            x = moving_average(x)
            if ikey == 0:
                for key in annotations:
                    if not isinstance(anno_dict[key][0], str):
                        anno_dict[key] = moving_average(anno_dict[key])                 
        
        
        if normalize_to_zero_one:
            x = np.clip(x, np.quantile(x, 0.00), np.quantile(x, 0.95))
            x -= np.min(x)
            x /= np.max(x)
        X.append(x)
        
        if not as_heatmap:
            ax.plot(x[xlim[0] : xlim[1]], label=key)
        if ikey == 0:
            for igroup, group in enumerate(nodes):
                if len(groups_names) > 0 and group not in groups_names:
                    label = groups_names[group]
                else:
                    label = group
                x_tick_labels.append(label)
                
    X = np.asarray(X).squeeze()
    #import scipy.stats
    #X = scipy.stats.zscore(X)
    # X -= np.min(X, axis=0)
    # X /= np.max(X, axis=0)
    
    if as_heatmap:
        img = ax.imshow(X, aspect='auto', interpolation='nearest', cmap=color_map)
        if show_yticks:
            ax.set_yticks(range(len(X)))
            ax.set_yticklabels(keys, fontsize=ytick_fontsize)
        else:
            ax.set_yticks([])
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.tick_params(axis='both', which='both', length=0)
        ax.grid(False)
        if show_colorbar:
            pl.colorbar(img, ax=ax)
        left_margin = 0.2 if left_margin is None else left_margin
        pl.subplots_adjust(left=left_margin)
    else:
        left_margin = 0.4 if left_margin is None else left_margin
        if len(keys) > 1:
            pl.legend(
                frameon=False,
                loc='center left',
                bbox_to_anchor=(-left_margin, 0.5),
                fontsize=legend_fontsize,
            )
    xlabel = groups_key
    if not as_heatmap:
        ax.set_xlabel(xlabel)
        pl.yticks([])
        if len(keys) == 1:
            pl.ylabel(keys[0] + ' (a.u.)')
    else:
        import matplotlib.colors

        # groups bar
        ax_bounds = ax.get_position().bounds
        groups_axis = pl.axes(
            (
                ax_bounds[0],
                ax_bounds[1] - ax_bounds[3] / len(keys),
                ax_bounds[2],
                ax_bounds[3] / len(keys),
            )
        )
        groups = np.array(groups)[None, :]
        groups_axis.imshow(
            groups,
            aspect='auto',
            interpolation="nearest",
            cmap=matplotlib.colors.ListedColormap(
                # the following line doesn't work because of normalization
                # adata.uns['paga_groups_colors'])
                palette_groups[np.min(groups).astype(int) :],
                N=int(np.max(groups) + 1 - np.min(groups)),
            ),
        )
        if show_yticks:
            groups_axis.set_yticklabels(['', xlabel, ''], fontsize=ytick_fontsize)
        else:
            groups_axis.set_yticks([])
        groups_axis.set_frame_on(False)
        if show_node_names:
            ypos = (groups_axis.get_ylim()[1] + groups_axis.get_ylim()[0]) / 2
            x_tick_locs = _sc_utils.moving_average(x_tick_locs, n=2)
            for ilabel, label in enumerate(x_tick_labels):
                groups_axis.text(
                    x_tick_locs[ilabel],
                    ypos,
                    x_tick_labels[ilabel],
                    fontdict=dict(
                        horizontalalignment='center',
                        verticalalignment='center',
                    ),
                )
        groups_axis.set_xticks([])
        groups_axis.grid(False)
        groups_axis.tick_params(axis='both', which='both', length=0)
        # further annotations
        y_shift = ax_bounds[3] / len(keys)
        for ianno, anno in enumerate(annotations):
            if ianno > 0:
                y_shift = ax_bounds[3] / len(keys) / 2
            anno_axis = pl.axes(
                (
                    ax_bounds[0],
                    ax_bounds[1] - (ianno + 2) * y_shift,
                    ax_bounds[2],
                    y_shift,
                )
            )
            arr = np.array(anno_dict[anno])[None, :]
            if anno not in color_maps_annotations:
                color_map_anno = (
                    'Vega10' if is_categorical_dtype(adata.obs[anno]) else 'Greys'
                )
            else:
                color_map_anno = color_maps_annotations[anno]
            img = anno_axis.imshow(
                arr,
                aspect='auto',
                interpolation='nearest',
                cmap=color_map_anno,
            )
            if show_yticks:
                anno_axis.set_yticklabels(['', anno, ''], fontsize=ytick_fontsize)
                anno_axis.tick_params(axis='both', which='both', length=0)
            else:
                anno_axis.set_yticks([])
            anno_axis.set_frame_on(False)
            anno_axis.set_xticks([])
            anno_axis.grid(False)
    if title is not None:
        ax.set_title(title, fontsize=title_fontsize)
    if show is None and not ax_was_none:
        show = False
    else:
        show = settings.autoshow if show is None else show
    _utils.savefig_or_show('paga_path', show=show, save=save)
    if return_data:
        df = pd.DataFrame(data=X.T, columns=keys)
        df['groups'] = moving_average(groups)  # groups is without moving average, yet
        if 'dpt_pseudotime' in anno_dict:
            df['distance'] = anno_dict['dpt_pseudotime'].T
        return ax, df if ax_was_none and not show else df
    else:
        return ax if ax_was_none and not show else None

