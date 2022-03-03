from typing import List


def _get_cortex_dict() -> dict:
    '''Returns dictionary of {eight cortical lobe: labels}'''
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

    tmp_dict = {'OFC': ofc, 'MPFC': mpfc, 'LPFC': lpfc, 'SMC': smc,
                'MTC': mtc, 'LTC': ltc, 'PC': pc, 'OCC': occ}

    return tmp_dict


def mark_cortex(cortex_to_highlight: List[str],
                data_dk: dict = None,
                val: float = 1) -> dict:
    '''Return dict of {label: value} for a list of bilateral cortical lobes

    Key Arguments:
        - cortex_to_highlight: list of cortical lobes to highlight, list of
                               str. Required.
                               ['OFC', 'MPFC', 'LPFC', 'SMC',
                                'PC', 'MTC', 'LTC', 'OCC']
        - data_dk: dictionary of {label: value}, None by default.
        - val: value to highlight the given cortex with, float.

    Returns:
        - dictionary of {label: value}, None by default.
    '''
    cortex_dict = _get_cortex_dict()
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


def plot_lobes(lobes_dict: dict, hatches: list, cmap_name: str):
    '''Plot lobes'''

    data_dict = dict(zip(
        [x for x in lobes_dict.keys()],
        [mark_cortex(y) for y in lobes_dict.values()],
        ))

    background = 'w'
    edgecolor = 'k'
    figsize = (15, 15)
    title = 'test'
    fontsize = 15
    bordercolor = 'w'

    import os.path as op
    import sys
    sys.path.append('/Users/kc244/tmp/python-ggseg')
    import ggseg
    from glob import glob
    import matplotlib.pyplot as plt

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


    cmap, norm = ggseg._get_cmap_(cmap_name, 1, vminmax=[0, 1])
    for num, (key, data) in enumerate(data_dict.items()):
        hatch = hatches[num]
        if num == 0:
            ggseg._render_data_(data, wd,
                                cmap, norm, ax, edgecolor, hatch=hatch)
        else:
            # iFW
            edgecolor='white'
            ggseg._render_data_(data, wd, cmap,
                                norm, ax, edgecolor,
                                hatch=hatch, fill=False, alpha=1)

            # iFW cortex boundary
            edgecolor='black'
            ggseg._render_data_(data, wd, cmap,
                                norm, ax, edgecolor,
                                hatch=None, fill=False, alpha=1)

    # DKT regions with no provided values are rendered in gray
    def flatten(t):
        return [item for sublist in t for item in sublist]
    data_regions = list(
            set(flatten([list(x.keys()) for x in data_dict.values()])))
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    ggseg._render_regions_(files, ax, 'gray', edgecolor)

    # # A colorbar is added
    ggseg._add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, 'ha')

    # plt.savefig('followup_chrp.png', dpi=300)

    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 0.5
    import matplotlib.pyplot as plt
    plt.show()


def test_longitudinal_new():
    volume_longitudinal_diff = ['OFC', 'LPFC', 'MPFC', 'LTC',
                                'SMC', 'PC', 'OCC']
    ifw_longitudinal_diff = ['MPFC', 'LTC', 'SMC', 'PC', 'OCC']
    ifw_longitudinal_diff_chrp = ['MPFC', 'LTC', 'SMC']


    lobes_dict = {'volume': volume_longitudinal_diff,
                  'ifw': ifw_longitudinal_diff,
                  'ifw_chrp': ifw_longitudinal_diff_chrp}

    hatches = [None, 'O', 'O.']
    plot_lobes(lobes_dict, hatches, 'Reds')


