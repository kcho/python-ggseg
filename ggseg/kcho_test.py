import sys
sys.path.append('/Users/kc244/python-ggseg')
import ggseg
import os.path as op
import matplotlib.pyplot as plt
from glob import glob

import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.patches import Patch

from typing import List

def get_cortex_dict():
    ofc = ['parsorbitalis', 'medialorbitofrontal', 'lateralorbitofrontal']
    mpfc = ['caudalanteriorcingulate', 'rostralanteriorcingulate',
            'superiorfrontal']
    lpfc = ['parstriangularis', 'rostralmiddlefrontal', 'frontalpole',
            'parsopercularis']
    smc = ['precentral', 'caudalmiddlefrontal', 'postcentral', 'paracentral']
    pc = ['inferiorparietal', 'supramarginal', 'precuneus',
          'posteriorcingulate', 'isthmuscingulate', 'superiorparietal']
    mtc = ['entorhinal', 'parahippocampal', 'fusiform']
    ltc = ['transversetemporal', 'superiortemporal', 'bankssts',
           'inferiortemporal', 'middletemporal', 'temporalpole']
    occ = ['pericalcarine', 'lingual', 'lateraloccipital', 'cuneus']

    tmp_dict = {'OFC': ofc, 'MPFC': mpfc,
                'LPFC': lpfc, 'SMC': smc,
                'MTC': mtc, 'LTC': ltc,
                'PC': pc, 'OCC': occ}

    return tmp_dict


def mark_cortex(cortex_to_highlight: List[str], data_dk: dict = None,
                val: float = 1):

    
    cortex_dict = get_cortex_dict()
    if data_dk is None:
        data_dk = {}

    for cortex, labels in cortex_dict.items():
        if cortex in cortex_to_highlight:
            for label in labels:
                for side in 'left', 'right':
                    if f'{label}_right' in data_dk.keys():
                        data_dk[f'{label}_{side}'] += val
                    else:
                        data_dk[f'{label}_{side}'] = val
        # else:
            # for label in labels:
                # data_dk[f'{label}_left'] = 1
                # data_dk[f'{label}_right'] = 1

    return data_dk


def mark_cortex_with_f_value(cortex_to_highlight: List[str], data_dk: dict = None,
                val: float = 1):

    
    cortex_dict = get_cortex_dict()
    if data_dk is None:
        data_dk = {}

    for cortex, labels in cortex_dict.items():
        if cortex in cortex_to_highlight:
            for label in labels:
                for side in 'left', 'right':
                    if f'{label}_right' in data_dk.keys():
                        data_dk[f'{label}_{side}'] += val
                    else:
                        data_dk[f'{label}_{side}'] = val

    return data_dk


def render_data(data, wd, cmap, ax, edgecolor, vminmax, hatch=None,
        facealpha=1, alpha=1,
        fill=None, title='', white_bg=False):
    # For every region with a provided value, we draw a patch with the color
    # matching the normalized scale
    cmap, norm = ggseg._get_cmap_(cmap, data.values(), vminmax=vminmax)
    for k, v in data.items():
        fp = op.join(wd, k)
        if op.isfile(fp):
            p = open(fp).read()
            codes, verts = ggseg._svg_parse_(p)
            path = Path(verts, codes)

            c = cmap(norm(v))


            c_new = (c[0], c[1], c[2], facealpha)
            print(c_new)
            if white_bg:
                c_new = 'white'

            if fill is None:
                ax.add_patch(patches.PathPatch(path, facecolor=c_new,
                             edgecolor=edgecolor, lw=1, hatch=hatch, alpha=alpha))
            if fill is False:
                ax.add_patch(patches.PathPatch(path, facecolor=c_new, fill=False,
                             edgecolor=edgecolor, lw=1, hatch=hatch, alpha=alpha))
            else:
                ax.add_patch(patches.PathPatch(path, facecolor=c_new,
                             edgecolor=edgecolor, lw=1, hatch=hatch, alpha=alpha,
                             fill=fill))
        else:
            print('%s not found' % fp)
            pass

    # return Patch(facecolor=cmap(norm(1)), edgecolor=edgecolor, hatch=hatch,
                # label=title, alpha=alpha)
    return Patch(facecolor=c_new, edgecolor=edgecolor, hatch=hatch,
                label=title, alpha=alpha)

def render_data_wo_norm(data, wd, c, ax, edgecolor,
        hatch=None, alpha=1,
        fill=None, title=''):
    # For every region with a provided value, we draw a patch with the color
    # matching the normalized scale
    # cmap, norm = ggseg._get_cmap_(cmap, data.values(), vminmax=vminmax)
    for k, v in data.items():
        fp = op.join(wd, k)
        if op.isfile(fp):
            p = open(fp).read()
            codes, verts = ggseg._svg_parse_(p)
            path = Path(verts, codes)
            # c = cmap(norm(v))
            if fill is None:
                ax.add_patch(patches.PathPatch(path, facecolor=c,
                             edgecolor=edgecolor, lw=1,
                             hatch=hatch,
                             alpha=alpha))
            else:
                ax.add_patch(patches.PathPatch(path, facecolor=c,
                             edgecolor=edgecolor, lw=1, hatch=hatch, alpha=alpha,
                             fill=fill))
        else:
            print('%s not found' % fp)
            pass

    return Patch(facecolor=c, edgecolor=edgecolor, hatch=hatch, label=title)


def plot_dk_ifw_project():
    '''Plot DK for iFW project result visualization'''
    draw_baseline()


def draw_info(conf: dict):
    background = 'w'
    edgecolor = 'k'
    figsize = (15, 15)
    title = ''
    fontsize = 15
    bordercolor = 'w'

    wd = op.join(op.dirname(ggseg.__file__), 'data', 'dk')

    # A figure is created by the joint dimensions of the whole-brain outlines
    whole_reg = ['lateral_left', 'medial_left',
                 'lateral_right', 'medial_right']

    files = [open(op.join(wd, e)).read() for e in whole_reg]

    ax = ggseg._create_figure_(
            files, figsize, background,
            title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    ggseg._render_regions_(files, ax, bordercolor, edgecolor)

    # legend
    legend_elements = []
    # ifw_patch, volume_patch, overlap_patch]
    for modality, modality_dict in conf.items():
        if modality != 'overlap':
            patch = render_data(conf[modality]['data'], wd,
                                conf[modality]['cmap'], ax,
                                conf[modality]['edgecolor'],
                                conf[modality]['vminmax'],
                                hatch=conf[modality]['hatch'],
                                alpha=conf[modality]['alpha'],
                                title=conf[modality]['title'])
            legend_elements.append(patch)

        else:
            if conf['ifw']['hatch'] is None:
                overlap_hatch = conf['volume']['hatch']
            elif conf['volume']['hatch'] is None:
                overlap_hatch = conf['ifw']['hatch']
            elif conf['ifw']['hatch'] is None and \
                    conf['volume']['hatch'] is None:
                overlap_hatch = None
            else:
                overlap_hatch = conf['volume']['hatch'] + \
                                conf['ifw']['hatch']

            patch = render_data(conf['overlap']['data'], wd,
                                conf['overlap']['cmap'], ax,
                                conf['overlap']['edgecolor'],
                                conf['overlap']['vminmax'],
                                hatch=overlap_hatch,
                                alpha=conf['overlap']['alpha'],
                                title=conf['overlap']['title'])
            legend_elements.append(patch)

    # DKT regions with no provided values are rendered in gray
    data_regions = list(conf['ifw']['data'].keys()) + \
                   list(conf['volume']['data'].keys())
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    ggseg._render_regions_(files, ax, 'gray', edgecolor)

    ax.legend(handles=legend_elements,
              loc='center', prop={'size': 18})


def draw_iFW_group_time_interaction():
    volume_diff = ['MPFC', 'LTC', 'SMC']
    ifw_diff = ['MPFC', 'LTC', 'SMC', 'PC', 'OCC']
    overlap = [x for x in volume_diff if x in ifw_diff]

    volume_data = mark_cortex(volume_diff)
    ifw_data = mark_cortex(ifw_diff)
    overlap_data = mark_cortex(overlap)

    conf = {'ifw': {
                'data': ifw_data,
                'cmap': 'Blues',
                'edgecolor': 'k',
                'vminmax': [0, 4],
                'hatch': 'O.',
                'alpha': 0.8,
                'title': 'CHR-P > CHR-NP'},
            'volume': {
                'data': volume_data,
                'cmap': 'Blues',
                'edgecolor': 'k',
                'vminmax': [0, 2],
                'hatch': '/',
                'alpha': 1,
                'title': 'CHR-P > HC'},
            'overlap': {
                'data': overlap_data,
                'cmap': 'Blues',
                'edgecolor': 'k',
                'vminmax': [0, 1],
                'hatch': None,
                'alpha': 0.5,
                'title': 'Overlap between the regions'}}
    draw_info(conf)
    plt.savefig(f'ifw_group_time.png', dpi=300)
    plt.show()


def draw_volume_group_effect():
    volume_diff = ['MTC']
    ifw_diff = ['MTC', 'LTC']
    overlap = [x for x in volume_diff if x in ifw_diff]

    volume_data = mark_cortex(volume_diff)
    ifw_data = mark_cortex(ifw_diff)
    overlap_data = mark_cortex(overlap)

    conf = {'ifw': {
                'data': ifw_data,
                'cmap': 'Reds',
                'edgecolor': 'k',
                'vminmax': [0, 4],
                'hatch': '//',
                'alpha': 0.8,
                'title': 'CHR-P > HC'},
            'volume': {
                'data': volume_data,
                'cmap': 'Reds',
                'edgecolor': 'k',
                'vminmax': [0, 2],
                'hatch': None,
                'alpha': 1,
                'title': 'CHR-NP > HC'},
            'overlap': {
                'data': overlap_data,
                'cmap': 'Reds',
                'edgecolor': 'k',
                'vminmax': [0, 1],
                'hatch': None,
                'alpha': 0.5,
                'title': 'Overlap between the regions'}}
    draw_info(conf)
    plt.savefig(f'volume_group_effect.png', dpi=300)
    plt.show()


def draw_baseline():
    volume_baseline_diff = ['LTC', 'MTC']
    ifw_baseline_diff = ['MPFC', 'LTC', 'PC', 'OCC']
    overlap = [x for x in volume_baseline_diff if x in ifw_baseline_diff]

    volume_baseline_data = mark_cortex(volume_baseline_diff)
    ifw_baseline_data = mark_cortex(ifw_baseline_diff)
    overlap_baseline_data = mark_cortex(overlap)

    # conf = {'ifw': {
                # 'data': ifw_baseline_data,
                # 'cmap': 'Blues',
                # 'edgecolor': 'k',
                # 'vminmax': [0, 1.3],
                # 'hatch': 'O.',
                # 'alpha': 0.8,
                # 'title': 'Significantly increased iFW'},
            # 'volume': {
                # 'data': volume_baseline_data,
                # 'cmap': 'Reds',
                # 'edgecolor': 'k',
                # 'vminmax': [0, 1.2],
                # 'hatch': None,
                # 'alpha': 1,
                # 'title': 'Significantly reduced volume'},
            # 'overlap': {
                # 'data': overlap_baseline_data,
                # 'cmap': 'Reds',
                # 'edgecolor': 'k',
                # 'vminmax': [0, 4],
                # 'hatch': None,
                # 'alpha': 1,
                # 'title': 'Overlap between the regions'}}

    conf = {'ifw': {
                'data': ifw_baseline_data,
                'cmap': 'Blues',
                'edgecolor': 'k',
                'vminmax': [0, 4],
                'hatch': 'O.',
                'alpha': 0.8,
                'title': 'Significantly increased iFW'},
            'volume': {
                'data': volume_baseline_data,
                'cmap': 'Reds',
                'edgecolor': 'k',
                'vminmax': [0, 1.2],
                'hatch': None,
                'alpha': 1,
                'title': 'Significantly reduced volume'},
            'overlap': {
                'data': overlap_baseline_data,
                'cmap': 'Reds',
                'edgecolor': 'k',
                'vminmax': [0, 3],
                'hatch': None,
                'alpha': 0.5,
                'title': 'Overlap between the regions'}}
    draw_info(conf)
    plt.savefig(f'baseline_test.png', dpi=300)
    plt.show()


def draw_followup():
    volume_longitudinal_diff = ['OFC', 'LPFC', 'MPFC', 'LTC',
                                'SMC', 'PC', 'OCC']
    ifw_longitudinal_diff = ['MPFC', 'LTC', 'SMC', 'PC', 'OCC']
    overlap = [x for x in volume_longitudinal_diff 
            if x in ifw_longitudinal_diff]

    volume_longitudinal_data = mark_cortex(volume_longitudinal_diff)
    ifw_longitudinal_data = mark_cortex(ifw_longitudinal_diff)
    overlap_longitudinal_data = mark_cortex(overlap)

    conf = {'ifw': {
                'data': ifw_longitudinal_data,
                'cmap': 'Greens',
                'edgecolor': 'k',
                'vminmax': [0, 1.8],
                'hatch': 'O.',
                'alpha': 0.8,
                'title': 'Significantly increased iFW'},
            'volume': {
                'data': volume_longitudinal_data,
                'cmap': 'Reds',
                'edgecolor': 'k',
                'vminmax': [0, 1.8],
                'hatch': None,
                'alpha': 1,
                'title': 'Significantly reduced volume'},
            'overlap': {
                'data': overlap_longitudinal_data,
                'cmap': 'YlOrRd_r',
                'edgecolor': 'k',
                'vminmax': [0, 1.3],
                'hatch': None,
                'alpha': 0.5,
                'title': 'Overlap between the regions'}}
    # conf = {'ifw': {
                # 'data': ifw_longitudinal_data,
                # 'cmap': 'Blues',
                # 'edgecolor': 'k',
                # 'vminmax': [0, 1.3],
                # 'hatch': 'O.',
                # 'alpha': 0.8,
                # 'title': 'Significantly increased iFW'},
            # 'volume': {
                # 'data': volume_longitudinal_data,
                # 'cmap': 'Reds',
                # 'edgecolor': 'k',
                # 'vminmax': [0, 1.2],
                # 'hatch': None,
                # 'alpha': 1,
                # 'title': 'Significantly reduced volume'},
            # 'overlap': {
                # 'data': overlap_longitudinal_data,
                # 'cmap': 'Reds',
                # 'edgecolor': 'k',
                # 'vminmax': [0, 3],
                # 'hatch': None,
                # 'alpha': 0.5,
                # 'title': 'Overlap between the regions'}}

    draw_info(conf)
    plt.show()
    # plt.savefig(f'followup_test.png', dpi=300)


def test_flow():
    draw_volume_group_effect()
    # draw_iFW_group_time_interaction()
    # draw_baseline()
    # draw_followup()


def test_draw_labels():
    from matplotlib.pyplot import cm
    import numpy as np
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap

    color = iter(cm.rainbow(np.linspace(0, 1, 8)))

    data_dict = {}
    roi_list = ['OFC', 'LPFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC', 'OCC']

    for num, roi in enumerate(roi_list):
        data_dict[roi] = {}
        data_dict[roi]['num'] = num
        data_dict[roi]['color'] = next(color)
        data_dict[roi]['data'] = mark_cortex([roi])


    background = 'w'
    edgecolor = 'k'
    cmap = 'coolwarm'
    vminmax = [0, 1]
    figsize = (15, 15)
    title = 'test'
    fontsize = 15
    bordercolor = 'w'

    
    wd = op.join(op.dirname(ggseg.__file__), 'data', 'dk')

    # A figure is created by the joint dimensions of the whole-brain outlines
    whole_reg = ['lateral_left', 'medial_left',
                 'lateral_right', 'medial_right']

    files = [open(op.join(wd, e)).read() for e in whole_reg]

    # for color_t in ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired','Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn','RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r','gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean','ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']:
    ax = ggseg._create_figure_(
            files, figsize, background,
            title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    # ggseg._render_regions_(files, ax, bordercolor, edgecolor)

    # DKT regions with no provided values are rendered in gray
    data_regions = roi_list
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    # ggseg._render_regions_(files, ax, 'gray', edgecolor)

    volume_patches = []
    for roi, data in data_dict.items():
        # for var, value in data.items():
        # print(data['data'])
        volume_patch = render_data_wo_norm(
                data['data'], wd, data['color'],
                ax, edgecolor, alpha=0.6,
                title=roi, fill=False)
        volume_patches.append(volume_patch)


    # # A colorbar is added
    # _add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, ylabel)
    # cmap, norm = ggseg._get_cmap_('Blues_r', [0, 1], vminmax=[0, 1.1])

    # legend_elements = [volume_patch]
    legend_elements = volume_patches
    ax.legend(
            bbox_to_anchor=(0.47, 0.5),
            handles=legend_elements,
            loc='center',
            prop={'size': 18})

    plt.savefig(f'roi_visualization_legend.png', dpi=300)
    plt.show()


def test_baseline_ofer_grant():
    '''Group effect baseline posthoc'''

    all_rois = ['OFC', 'LPFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC', 'OCC']
    # set data
    volume_baseline = ['LTC']
    ifw_baseline = ['MPFC', 'LTC', 'PC', 'OCC']

    all_selected = volume_baseline + ifw_baseline
    # rest = [x for x in all_rois if x not in all_selected]
    rest = all_rois
    print(rest)
    # sys.exit()

    # overlap = [x for x in volume_baseline_diff if x in ifw_baseline_diff]

    # set hash
    ifw_increase_hash = '+O.'
    ifw_increase_hash = '.O'
    ifw_increase_alpha = 0.2
    ifw_increase_color = 'Blues'
    # ifw_increase_color = 'White'
    volume_reduction_hash = None
    volume_reduction_alpha = 0.7
    volume_reduction_color = 'Reds'

    # set group settings
    chrp_alpha = 1
    chrnp_alpha = 1

    
    volume_baseline_data = mark_cortex(volume_baseline)
    ifw_baseline_data = mark_cortex(ifw_baseline)
    rest_data = mark_cortex(rest)


    # data_dk = mark_cortex(['OFC', 'LPFC', 'MPFC', 'LTC', 'SMC', 'PC', 'OCC'])
    # data_dk_2 = mark_cortex(['MPFC', 'LTC', 'SMC', 'PC', 'OCC'], val=-1)
    background = 'w'
    edgecolor = 'k'
    cmap = 'coolwarm'
    vminmax = [-1, 1]
    figsize = (15, 15)
    title = 'test'
    fontsize = 15
    bordercolor = 'w'

    
    wd = op.join(op.dirname(ggseg.__file__), 'data', 'dk')

    # A figure is created by the joint dimensions of the whole-brain outlines
    whole_reg = ['lateral_left', 'medial_left',
                 'lateral_right', 'medial_right']

    files = [open(op.join(wd, e)).read() for e in whole_reg]

    # for color_t in ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired','Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn','RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r','gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean','ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']:
    ax = ggseg._create_figure_(
            files, figsize, background,
            title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    ggseg._render_regions_(files, ax, bordercolor, edgecolor)

    # ifw_patch = render_data(ifw_baseline_data, wd, 'Blues',
                # ax, edgecolor, [0, 1.3], hatch='O.', alpha=0.8,
                # title='Baseline iFW: (CHR-NP & CHR-P) > HC')
    # volume_patch = render_data(volume_baseline_reduction_CHRNP_HC_data,
            # wd, volume_reduction_color,
                # ax, edgecolor, [0, 1.2],
                # facealpha=volume_reduction_alpha,
                # alpha=chrnp_alpha,
                # hatch=volume_reduction_hash,
                # title='Baseline volume: CHR-NP < HC')


    # rest_patch = render_data(rest_data, wd,
            # 'Blues',
            # ax,
            # edgecolor, [0, 100],
            # facealpha=1,
            # alpha=1)

    # try without filling it
    ifw_patch = render_data(ifw_baseline_data, wd,
            ifw_increase_color,
            ax,
            edgecolor, [0, 100],
            facealpha=1,
            hatch=ifw_increase_hash,
            alpha=1,
            title='Baseline iFW: CHR > HC')

    volume_patch = render_data(volume_baseline_data,
            wd, volume_reduction_color,
                ax, edgecolor, [0, 3],
                alpha=1,
                facealpha=volume_reduction_alpha,
            hatch=ifw_increase_hash,
                title='Baseline volume: CHR < HC')
                # hatch=volume_reduction_hash,


    # try without filling it
    # ifw_patch = render_data(ifw_baseline_reduction_CHRNP_HC_data, wd,
            # ifw_increase_color,
            # ax, edgecolor, [0, 1.3],
                # facealpha=ifw_increase_alpha,
            # hatch=ifw_increase_hash,
            # alpha=chrp_alpha,
            # white_bg=True,
            # title='Baseline iFW: CHR-NP > HC')

    # # render_data(overlap_baseline_data, wd, 'YlGn_r',
                # # ax, edgecolor, [0, 3], hatch='O.', alpha=0.7)

    # color_t = 'Reds'
    # overlap_patch = render_data(overlap_baseline_data, wd, color_t,
                # ax, edgecolor, [0, 3], hatch=none, alpha=0.5,
                # title='overlap between the regions')

    # render_data(overlap_baseline_data, wd, color_t,
                # ax, edgecolor, [0, 3], hatch='O.', alpha=0.5, fill=False)

    # DKT regions with no provided values are rendered in gray
    # data_regions = list(volume_baseline_data.keys()) + \
            # list(ifw_baseline_data.keys())
    data_regions = list(rest_data.keys())
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    # ggseg._render_regions_(files, ax, 'gray', edgecolor)
    ggseg._render_regions_(files, ax, 'white', edgecolor)

    # # A colorbar is added
    # _add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, ylabel)

    cmap, norm = ggseg._get_cmap_('Blues_r', [0, 1], vminmax=[0, 1.1])
    # legend_elements = [ifw_patch, volume_patch, overlap_patch]

    # ax.legend(handles=legend_elements, loc='center', prop={'size': 18})

    # plt.savefig(f'baseline_posthoc_ofer_grant.png', dpi=300)
    plt.show()


def test_baseline_posthoc():
    '''Group effect baseline posthoc'''

    # set data
    volume_baseline_reduction_CHRP_HC = ['LTC', 'MTC']
    volume_baseline_reduction_CHRNP_HC = ['MTC']
    ifw_baseline_reduction_CHRP_HC = ['MPFC', 'LTC', 'PC', 'OCC']
    ifw_baseline_reduction_CHRNP_HC = ['MPFC', 'LTC', 'PC', 'OCC']

    # set hash
    ifw_increase_hash = '+O.'
    ifw_increase_hash = '+'
    ifw_increase_alpha = 0.2
    ifw_increase_color = 'Blues'
    # ifw_increase_color = 'White'
    volume_reduction_hash = None
    volume_reduction_alpha = 0.7
    volume_reduction_color = 'Reds'

    # set group settings
    chrp_alpha = 1
    chrnp_alpha = 1

    
    volume_baseline_reduction_CHRP_HC_data = mark_cortex(volume_baseline_reduction_CHRP_HC)
    volume_baseline_reduction_CHRNP_HC_data = mark_cortex(volume_baseline_reduction_CHRNP_HC)
    ifw_baseline_reduction_CHRP_HC_data = mark_cortex(ifw_baseline_reduction_CHRP_HC)
    ifw_baseline_reduction_CHRNP_HC_data = mark_cortex(ifw_baseline_reduction_CHRNP_HC)


    # data_dk = mark_cortex(['OFC', 'LPFC', 'MPFC', 'LTC', 'SMC', 'PC', 'OCC'])
    # data_dk_2 = mark_cortex(['MPFC', 'LTC', 'SMC', 'PC', 'OCC'], val=-1)
    background = 'w'
    edgecolor = 'k'
    cmap = 'coolwarm'
    vminmax = [-1, 1]
    figsize = (15, 15)
    title = 'test'
    fontsize = 15
    bordercolor = 'w'

    
    wd = op.join(op.dirname(ggseg.__file__), 'data', 'dk')

    # A figure is created by the joint dimensions of the whole-brain outlines
    whole_reg = ['lateral_left', 'medial_left',
                 'lateral_right', 'medial_right']

    files = [open(op.join(wd, e)).read() for e in whole_reg]

    # for color_t in ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired','Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn','RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r','gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean','ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']:
    ax = ggseg._create_figure_(
            files, figsize, background,
            title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    ggseg._render_regions_(files, ax, bordercolor, edgecolor)

    # ifw_patch = render_data(ifw_baseline_data, wd, 'Blues',
                # ax, edgecolor, [0, 1.3], hatch='O.', alpha=0.8,
                # title='Baseline iFW: (CHR-NP & CHR-P) > HC')
    volume_patch = render_data(volume_baseline_reduction_CHRP_HC_data,
            wd, volume_reduction_color,
                ax, edgecolor, [0, 1.2],
                alpha=chrp_alpha,
                facealpha=volume_reduction_alpha,
                hatch=volume_reduction_hash,
                title='Baseline volume: CHR-P < HC')

    volume_patch = render_data(volume_baseline_reduction_CHRNP_HC_data,
            wd, volume_reduction_color,
                ax, edgecolor, [0, 1.2],
                facealpha=volume_reduction_alpha,
                alpha=chrnp_alpha,
                hatch=volume_reduction_hash,
                title='Baseline volume: CHR-NP < HC')


    # try without filling it
    ifw_patch = render_data(ifw_baseline_reduction_CHRP_HC_data, wd,
            ifw_increase_color,
            ax, edgecolor, [0, 1.3],
                facealpha=ifw_increase_alpha,
            hatch=ifw_increase_hash,
            alpha=chrnp_alpha,
            white_bg=True,
            title='Baseline iFW: CHR-P > HC')

    # try without filling it
    ifw_patch = render_data(ifw_baseline_reduction_CHRNP_HC_data, wd,
            ifw_increase_color,
            ax, edgecolor, [0, 1.3],
                facealpha=ifw_increase_alpha,
            hatch=ifw_increase_hash,
            alpha=chrp_alpha,
            white_bg=True,
            title='Baseline iFW: CHR-NP > HC')

    # # render_data(overlap_baseline_data, wd, 'YlGn_r',
                # # ax, edgecolor, [0, 3], hatch='O.', alpha=0.7)

    # color_t = 'Reds'
    # overlap_patch = render_data(overlap_baseline_data, wd, color_t,
                # ax, edgecolor, [0, 3], hatch=None, alpha=0.5,
                # title='Overlap between the regions')

    # render_data(overlap_baseline_data, wd, color_t,
                # ax, edgecolor, [0, 3], hatch='O.', alpha=0.5, fill=False)

    # DKT regions with no provided values are rendered in gray
    data_regions = list(volume_baseline_reduction_CHRP_HC_data.keys()) + \
            list(volume_baseline_reduction_CHRNP_HC_data.keys()) + \
            list(ifw_baseline_reduction_CHRP_HC_data.keys()) + \
            list(ifw_baseline_reduction_CHRNP_HC_data.keys())
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    ggseg._render_regions_(files, ax, 'gray', edgecolor)

    # # A colorbar is added
    # _add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, ylabel)

    cmap, norm = ggseg._get_cmap_('Blues_r', [0, 1], vminmax=[0, 1.1])
    # legend_elements = [ifw_patch, volume_patch, overlap_patch]

    # ax.legend(handles=legend_elements, loc='center', prop={'size': 18})

    plt.savefig(f'baseline_posthoc.png', dpi=300)
    # plt.show()


def test_baseline():
    volume_baseline_diff = ['LTC', 'MTC']
    ifw_baseline_diff = ['MPFC', 'LTC', 'PC', 'OCC']
    overlap = [x for x in volume_baseline_diff if x in ifw_baseline_diff]

    volume_baseline_data = mark_cortex(volume_baseline_diff)
    ifw_baseline_data = mark_cortex(ifw_baseline_diff)
    overlap_baseline_data = mark_cortex(overlap)


    # data_dk = mark_cortex(['OFC', 'LPFC', 'MPFC', 'LTC', 'SMC', 'PC', 'OCC'])
    # data_dk_2 = mark_cortex(['MPFC', 'LTC', 'SMC', 'PC', 'OCC'], val=-1)
    background = 'w'
    edgecolor = 'k'
    cmap = 'coolwarm'
    vminmax = [-1, 1]
    figsize = (15, 15)
    title = 'test'
    fontsize = 15
    bordercolor = 'w'

    
    wd = op.join(op.dirname(ggseg.__file__), 'data', 'dk')

    # A figure is created by the joint dimensions of the whole-brain outlines
    whole_reg = ['lateral_left', 'medial_left',
                 'lateral_right', 'medial_right']

    files = [open(op.join(wd, e)).read() for e in whole_reg]

    # for color_t in ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired','Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn','RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r','gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean','ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']:
    ax = ggseg._create_figure_(
            files, figsize, background,
            title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    ggseg._render_regions_(files, ax, bordercolor, edgecolor)

    ifw_patch = render_data(ifw_baseline_data, wd, 'Blues',
                ax, edgecolor, [0, 1.3], hatch='O.', alpha=0.8,
                title='Significantly increased iFW')

    volume_patch = render_data(volume_baseline_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], alpha=1,
                title='Significantly reduced volume')

    # render_data(overlap_baseline_data, wd, 'YlGn_r',
                # ax, edgecolor, [0, 3], hatch='O.', alpha=0.7)

    color_t = 'Reds'
    overlap_patch = render_data(overlap_baseline_data, wd, color_t,
                ax, edgecolor, [0, 3], hatch=None, alpha=0.5,
                title='Overlap between the regions')

    render_data(overlap_baseline_data, wd, color_t,
                ax, edgecolor, [0, 3], hatch='O.', alpha=0.5, fill=False)

    # DKT regions with no provided values are rendered in gray
    data_regions = list(volume_baseline_data.keys()) + \
            list(ifw_baseline_data.keys())
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    ggseg._render_regions_(files, ax, 'gray', edgecolor)

    # # A colorbar is added
    # _add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, ylabel)

    cmap, norm = ggseg._get_cmap_('Blues_r', [0, 1], vminmax=[0, 1.1])
    legend_elements = [ifw_patch, volume_patch, overlap_patch]

    ax.legend(handles=legend_elements, loc='center', prop={'size': 18})

    plt.savefig(f'baseline.png', dpi=300)
    # plt.show()


def test_longitudinal():
    volume_longitudinal_diff = ['OFC', 'LPFC', 'MPFC', 'LTC',
                                'SMC', 'PC', 'OCC']
    ifw_longitudinal_diff = ['MPFC', 'LTC', 'SMC', 'PC', 'OCC']
    ifw_longitudinal_diff_chrp = ['MPFC', 'LTC', 'SMC']
    overlap = [x for x in volume_longitudinal_diff
            if x in ifw_longitudinal_diff]

    volume_longitudinal_data = mark_cortex(volume_longitudinal_diff)
    ifw_longitudinal_data = mark_cortex(ifw_longitudinal_diff)
    ifw_longitudinal_chrp_data = mark_cortex(ifw_longitudinal_diff_chrp)
    overlap_longitudinal_data = mark_cortex(overlap)


    # data_dk = mark_cortex(['OFC', 'LPFC', 'MPFC', 'LTC', 'SMC', 'PC', 'OCC'])
    # data_dk_2 = mark_cortex(['MPFC', 'LTC', 'SMC', 'PC', 'OCC'], val=-1)
    background = 'w'
    edgecolor = 'k'
    cmap = 'coolwarm'
    vminmax = [-1, 1]
    figsize = (15, 15)
    title = 'test'
    fontsize = 15
    bordercolor = 'w'

    wd = op.join(op.dirname(ggseg.__file__), 'data', 'dk')

    # A figure is created by the joint dimensions of the whole-brain outlines
    whole_reg = ['lateral_left', 'medial_left',
                 'lateral_right', 'medial_right']

    files = [open(op.join(wd, e)).read() for e in whole_reg]
    ax = ggseg._create_figure_(
            files, figsize, background,
            title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    ggseg._render_regions_(files, ax, bordercolor, edgecolor)

    render_data(ifw_longitudinal_data, wd, 'Blues',
                ax, edgecolor, [0, 1.3], hatch='O.', alpha=0.8)

    render_data(volume_longitudinal_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], alpha=1)


    render_data(overlap_longitudinal_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch=None, alpha=0.5)

    render_data(overlap_longitudinal_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch='O.', alpha=0.5, fill=False)

    render_data(ifw_longitudinal_chrp_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch='+O.', alpha=0.5)


    # DKT regions with no provided values are rendered in gray
    data_regions = list(volume_longitudinal_data.keys()) + \
            list(ifw_longitudinal_data.keys())
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    ggseg._render_regions_(files, ax, 'gray', edgecolor)

    # # A colorbar is added
    # _add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, ylabel)

    # plt.savefig('followup_chrp.png', dpi=300)

    plt.show()


def test_longitudinal_new():
    volume_longitudinal_diff = ['OFC', 'LPFC', 'MPFC', 'LTC',
                                'SMC', 'PC', 'OCC']
    ifw_longitudinal_diff = ['MPFC', 'LTC', 'SMC', 'PC', 'OCC']
    ifw_longitudinal_diff_chrp = ['MPFC', 'LTC', 'SMC']
    overlap = [x for x in volume_longitudinal_diff
            if x in ifw_longitudinal_diff]

    volume_longitudinal_data = mark_cortex(volume_longitudinal_diff)
    ifw_longitudinal_data = mark_cortex(ifw_longitudinal_diff)
    ifw_longitudinal_chrp_data = mark_cortex(ifw_longitudinal_diff_chrp)
    overlap_longitudinal_data = mark_cortex(overlap)


    # data_dk = mark_cortex(['OFC', 'LPFC', 'MPFC', 'LTC', 'SMC', 'PC', 'OCC'])
    # data_dk_2 = mark_cortex(['MPFC', 'LTC', 'SMC', 'PC', 'OCC'], val=-1)
    background = 'w'
    edgecolor = 'k'
    cmap = 'coolwarm'
    vminmax = [-1, 1]
    figsize = (15, 15)
    title = 'test'
    fontsize = 15
    bordercolor = 'w'

    wd = op.join(op.dirname(ggseg.__file__), 'data', 'dk')

    # A figure is created by the joint dimensions of the whole-brain outlines
    whole_reg = ['lateral_left', 'medial_left',
                 'lateral_right', 'medial_right']

    files = [open(op.join(wd, e)).read() for e in whole_reg]
    ax = ggseg._create_figure_(
            files, figsize, background,
            title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    ggseg._render_regions_(files, ax, bordercolor, edgecolor)

    # render_data(ifw_longitudinal_data, wd, 'Blues',
                # ax, edgecolor, [0, 1.3], hatch='O.', alpha=0.8)

    render_data(volume_longitudinal_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], alpha=1)


    # render_data(overlap_longitudinal_data, wd, 'Reds',
                # ax, edgecolor, [0, 1.2], hatch=None, alpha=0.5)

    # render_data(overlap_longitudinal_data, wd, 'Reds',
                # ax, edgecolor, [0, 1.2], hatch='O.', alpha=0.5, fill=False)

    # render_data(ifw_longitudinal_chrp_data, wd, 'Reds',
                # ax, edgecolor, [0, 1.2], hatch='+O.', alpha=0.5)


    # edit here
    edgecolor='white'
    # render_data(ifw_longitudinal_data, wd, 'Blues',
                # ax, edgecolor, [0, 1.3], hatch='O.', alpha=0.8)

    render_data(overlap_longitudinal_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch=None, alpha=1, fill=None)

    render_data(overlap_longitudinal_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch='O', alpha=1, fill=None)

    render_data(ifw_longitudinal_chrp_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch='O.', alpha=1, fill=None)


    edgecolor='black'
    # render_data(ifw_longitudinal_data, wd, 'Blues',
                # ax, edgecolor, [0, 1.3], hatch='O.', alpha=0.8)

    render_data(overlap_longitudinal_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch=None, alpha=1, fill=False)

    render_data(overlap_longitudinal_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch=None, alpha=1, fill=False)

    render_data(ifw_longitudinal_chrp_data, wd, 'Reds',
                ax, edgecolor, [0, 1.2], hatch=None, alpha=1, fill=False)


    # for hatch
    # render_data(overlap_longitudinal_data, wd, 'Reds',
                # ax, edgecolor, [0, 1.2], hatch='O', white_bg=True,
                # alpha=1, fill=None)

    # DKT regions with no provided values are rendered in gray
    data_regions = list(volume_longitudinal_data.keys()) + \
            list(ifw_longitudinal_data.keys())
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    ggseg._render_regions_(files, ax, 'gray', edgecolor)

    # # A colorbar is added
    # _add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, ylabel)

    # plt.savefig('followup_chrp.png', dpi=300)

    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 0.5
    plt.savefig('followup_new.png', dpi=300)
    # plt.show()


def tmp():
    data_dk = {'bankssts_left': 2.48,
               'caudalanteriorcingulate_left': 2.228,
               'caudalmiddlefrontal_left': 2.618,
               'cuneus_left': 2.137,
               'entorhinal_left': 3.332,
               'fusiform_left': 2.605,
               'inferiorparietal_left': 2.436,
               'inferiortemporal_left': 2.782,
               'isthmuscingulate_left': 2.138,
               'lateraloccipital_left': 2.258,
               'lateralorbitofrontal_left': 2.582,
               'lingual_left': 2.333,
               'medialorbitofrontal_left': 2.321,
               'middletemporal_left': 2.949,
               'parahippocampal_left': 2.731,
               'paracentral_left': 2.433,
               'parsopercularis_left': 2.561,
               'parsorbitalis_left': 2.715,
               'parstriangularis_left': 2.463,
               'pericalcarine_left': 1.978,
               'postcentral_left': 2.213,
               'posteriorcingulate_left': 2.206,
               'precentral_left': 2.73,
               'precuneus_left': 2.238,
               'rostralanteriorcingulate_left': 2.632,
               'rostralmiddlefrontal_left': 2.406,
               'superiorfrontal_left': 2.638,
               'superiorparietal_left': 2.244,
               'superiortemporal_left': 2.656,
               'supramarginal_left': 2.507,
               'frontalpole_left': 2.579,
               'temporalpole_left': 3.652,
               'transversetemporal_left': 2.567,
               'insula_left': 2.869,
               'bankssts_right': 2.579,
               'caudalanteriorcingulate_right': 2.501,
               'caudalmiddlefrontal_right': 2.649,
               'cuneus_right': 2.265,
               'entorhinal_right': 2.448,
               'fusiform_right': 2.602,
               'inferiorparietal_right': 2.424,
               'inferiortemporal_right': 2.609,
               'isthmuscingulate_right': 2.127,
               'lateraloccipital_right': 2.381,
               'lateralorbitofrontal_right': 2.533,
               'lingual_right': 2.424,
               'medialorbitofrontal_right': 2.266,
               'middletemporal_right': 2.741,
               'parahippocampal_right': 2.741,
               'paracentral_right': 2.45,
               'parsopercularis_right': 2.521,
               'parsorbitalis_right': 2.505,
               'parstriangularis_right': 2.442,
               'pericalcarine_right': 2.02,
               'postcentral_right': 2.171,
               'posteriorcingulate_right': 2.4,
               'precentral_right': 2.654,
               'precuneus_right': 2.348,
               'rostralanteriorcingulate_right': 2.541,
               'rostralmiddlefrontal_right': 2.362,
               'superiorfrontal_right': 2.642,
               'superiorparietal_right': 2.265,
               'superiortemporal_right': 2.587,
               'supramarginal_right': 2.459,
               'frontalpole_right': 2.551,
               'temporalpole_right': 3.486,
               'transversetemporal_right': 2.714,
               'insula_right': 2.994}

    ggseg.plot_dk(data_dk, background='w', edgecolor='k', cmap='hot')
    ggseg.plot_dk(data_dk, background='w', edgecolor='k', cmap='hot',
                  vminmax=[2.5, 3.5])
