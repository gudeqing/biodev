import os
from functools import partial
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from plotly.offline import plot as plt
import plotly.io as pio
plt = partial(plt, auto_open=False)


def get_color_pool(n):
    # https://plot.ly/ipython-notebooks/color-scales/
    import colorlover
    if n <= 8:
        if n <= 3:
            n = 3
        return colorlover.scales[str(n)]['qual']['Set2']
    if n <= 12:
        return colorlover.scales[str(n)]['qual']['Set3']

    import random
    random.seed(666)

    def get_random_color(pastel_factor=0.5):
        return [(x + pastel_factor) / (1.0 + pastel_factor) for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]]

    def color_distance(c1, c2):
        return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])

    def generate_new_color(existing_colors, pastel_factor=0.5):
        max_distance = None
        best_color = None
        for i in range(0, 100):
            color = get_random_color(pastel_factor=pastel_factor)
            # exclude some colors
            if np.absolute(np.array(color) - np.array([1, 1, 1])).sum() < 0.08:
                continue
            if not existing_colors:
                return color
            best_distance = min([color_distance(color, c) for c in existing_colors])
            if not max_distance or best_distance > max_distance:
                max_distance = best_distance
                best_color = color
        return best_color

    color_pool = []
    for i in range(0, n):
        color_pool.append(generate_new_color(color_pool, pastel_factor=0.9))
    color_pool = [(int(x * 255), int(y * 255), int(z * 255)) for x, y, z in color_pool]
    color_pool = sorted(color_pool, key=lambda x: (x[0], x[1], x[2]))
    return colorlover.to_rgb(color_pool)


def get_marker_pool(n):
    maker_pool = [
        0, 'circle', 100, 'circle-open', 200, 'circle-dot', 300,
        'circle-open-dot', 1, 'square', 101, 'square-open', 201,
        'square-dot', 301, 'square-open-dot', 2, 'diamond', 102,
        'diamond-open', 202, 'diamond-dot', 302,
        'diamond-open-dot', 3, 'cross', 103, 'cross-open', 203,
        'cross-dot', 303, 'cross-open-dot', 4, 'x', 104, 'x-open',
        204, 'x-dot', 304, 'x-open-dot', 5, 'triangle-up', 105,
        'triangle-up-open', 205, 'triangle-up-dot', 305,
        'triangle-up-open-dot', 6, 'triangle-down', 106,
        'triangle-down-open', 206, 'triangle-down-dot', 306,
        'triangle-down-open-dot', 7, 'triangle-left', 107,
        'triangle-left-open', 207, 'triangle-left-dot', 307,
        'triangle-left-open-dot', 8, 'triangle-right', 108,
        'triangle-right-open', 208, 'triangle-right-dot', 308,
        'triangle-right-open-dot', 9, 'triangle-ne', 109,
        'triangle-ne-open', 209, 'triangle-ne-dot', 309,
        'triangle-ne-open-dot', 10, 'triangle-se', 110,
        'triangle-se-open', 210, 'triangle-se-dot', 310,
        'triangle-se-open-dot', 11, 'triangle-sw', 111,
        'triangle-sw-open', 211, 'triangle-sw-dot', 311,
        'triangle-sw-open-dot', 12, 'triangle-nw', 112,
        'triangle-nw-open', 212, 'triangle-nw-dot', 312,
        'triangle-nw-open-dot', 13, 'pentagon', 113,
        'pentagon-open', 213, 'pentagon-dot', 313,
        'pentagon-open-dot', 14, 'hexagon', 114, 'hexagon-open',
        214, 'hexagon-dot', 314, 'hexagon-open-dot', 15,
        'hexagon2', 115, 'hexagon2-open', 215, 'hexagon2-dot',
        315, 'hexagon2-open-dot', 16, 'octagon', 116,
        'octagon-open', 216, 'octagon-dot', 316,
        'octagon-open-dot', 17, 'star', 117, 'star-open', 217,
        'star-dot', 317, 'star-open-dot', 18, 'hexagram', 118,
        'hexagram-open', 218, 'hexagram-dot', 318,
        'hexagram-open-dot', 19, 'star-triangle-up', 119,
        'star-triangle-up-open', 219, 'star-triangle-up-dot', 319,
        'star-triangle-up-open-dot', 20, 'star-triangle-down',
        120, 'star-triangle-down-open', 220,
        'star-triangle-down-dot', 320,
        'star-triangle-down-open-dot', 21, 'star-square', 121,
        'star-square-open', 221, 'star-square-dot', 321,
        'star-square-open-dot', 22, 'star-diamond', 122,
        'star-diamond-open', 222, 'star-diamond-dot', 322,
        'star-diamond-open-dot', 23, 'diamond-tall', 123,
        'diamond-tall-open', 223, 'diamond-tall-dot', 323,
        'diamond-tall-open-dot', 24, 'diamond-wide', 124,
        'diamond-wide-open', 224, 'diamond-wide-dot', 324,
        'diamond-wide-open-dot', 25, 'hourglass', 125,
        'hourglass-open', 26, 'bowtie', 126, 'bowtie-open', 27,
        'circle-cross', 127, 'circle-cross-open', 28, 'circle-x',
        128, 'circle-x-open', 29, 'square-cross', 129,
        'square-cross-open', 30, 'square-x', 130, 'square-x-open',
        31, 'diamond-cross', 131, 'diamond-cross-open', 32,
        'diamond-x', 132, 'diamond-x-open', 33, 'cross-thin', 133,
        'cross-thin-open', 34, 'x-thin', 134, 'x-thin-open', 35,
        'asterisk', 135, 'asterisk-open', 36, 'hash', 136,
        'hash-open', 236, 'hash-dot', 336, 'hash-open-dot', 37,
        'y-up', 137, 'y-up-open', 38, 'y-down', 138,
        'y-down-open', 39, 'y-left', 139, 'y-left-open', 40,
        'y-right', 140, 'y-right-open', 41, 'line-ew', 141,
        'line-ew-open', 42, 'line-ns', 142, 'line-ns-open', 43,
        'line-ne', 143, 'line-ne-open', 44, 'line-nw', 144,
        'line-nw-open'
    ]
    return sorted([x for x in maker_pool if type(x) == int])[: n]


def draw(fig: go.Figure, prefix='', outdir=os.getcwd(), formats=('html',),
         height:int=None, width:int=None, scale=3, desc=None):
    for format in formats:
        out_name = os.path.join(outdir, '{}.{}'.format(prefix, format))
        if format == 'html':
            plt(fig, filename=out_name)
        else:
            if format in ['svg', 'pdf']:
                scale = 1
            pio.write_image(fig, format=format, file=out_name, height=height, width=width, scale=scale)
        if desc:
            with open(out_name+'.describe.txt', 'w') as fw:
                fw.write(desc+'\n')


def redefine_gene_category(gene_type_pair_df):
    category = {
        'IG': ['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_V_gene'],
        'TCR': ['TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_D_gene'],
        'pseudogene': ['.*pseudogene'],
        'miRNA': ['miRNA'],
        'ncRNA': ['Mt_rRNA', 'Mt_tRNA', 'misc_RNA', 'rRNA', 'scRNA', 'snRNA',
                  'snoRNA', 'ribozyme', 'sRNA', 'scaRNA', 'vaultRNA'],
        'lncRNA': ['3prime_overlapping_ncRNA', 'antisense', 'antisense_RNA', 'bidirectional_promoter_lncRNA',
                   'lincRNA', 'macro_lncRNA', ' non_coding', 'processed_transcript',
                   'sense_intronic', 'sense_overlapping'],
        # To be Experimentally Confirmed.
        'TEC': ['TEC'],
        'protein_coding': ['protein_coding']
    }

    type_dict = dict()
    for k, v in category.items():
        for each in v:
            type_dict[each] = k

    gene_type_pair_df['GeneType'] = 'unknown'
    for gene in gene_type_pair_df.index:
        detail_type = gene_type_pair_df.loc[gene, 'gene_type']
        if detail_type.endswith('pseudogene'):
            gene_type_pair_df.loc[gene, 'GeneType'] = 'pseudogene'
        elif detail_type in type_dict:
            gene_type_pair_df.loc[gene, 'GeneType'] = type_dict[detail_type]
    result = gene_type_pair_df.loc[:, ['GeneType']]
    result.columns = ['gene_type']
    return result


def detected_gene_type_stat(expr, gene_type=None, lower=2, marker_genes=None, gene2symbol=None):
    exp_table = pd.read_csv(expr, header=0, index_col=0, sep=None, engine='python')
    if gene_type is None:
        gene_type = "/nfs2/database/gencode_v29/gene_type.txt"
    if marker_genes:
        marker_genes = set(x.strip() for x in open(marker_genes))
    if gene2symbol:
        gene2symbol = dict(x.strip().split()[:2] for x in open(gene2symbol))
    gene_type = pd.read_csv(gene_type, header=0, index_col=0, sep=None, engine='python')
    gene_type.columns = ['gene_type']
    gene_type = redefine_gene_category(gene_type)
    exp_table = exp_table.join(gene_type)
    exp_table.to_csv(expr+'.AddType.csv')
    stat_data = list()
    for sample in exp_table.columns[:-1]:
        expressed_genes = exp_table.loc[exp_table[sample] >= lower]
        if marker_genes:
            detected_markers = set(expressed_genes.index) & marker_genes
            if gene2symbol:
                detected_markers = [gene2symbol[x] if x in gene2symbol else x for x in detected_markers]
            # print(f'Detected {len(detected_markers)} in {sample}')
            with open('detected.marker.genes.xls', 'a') as f:
                f.write(f'{sample}\t{len(detected_markers)}\t{detected_markers}\n')
        stat_dict = expressed_genes.groupby('gene_type').size().to_dict()
        stat_dict['sample'] = sample
        stat_data.append(stat_dict)
    stat_data = pd.DataFrame(stat_data).set_index('sample').fillna(0)
    stat_data = stat_data.transpose()
    stat_data = stat_data.sort_values(list(stat_data.columns), ascending=False)
    stat_data = stat_data.transpose()
    stat_data.to_csv(expr+'.TypeStat.csv')
    # plot
    df = stat_data
    order = sorted(df.index)
    df = df.loc[order]
    # print(df.head())
    colors = get_color_pool(df.shape[1])
    color_dict = dict(zip(df.columns, colors))
    data = [go.Bar(x=df.index, y=df[x] / df.sum(axis=1), name=x, marker=dict(color=color_dict[x])) for x in df.columns]
    layout = go.Layout(
        title="Gene type distribution",
        # xaxis=dict(title='Sample'),
        barmode='stack',
    )
    outdir = os.path.dirname(expr)
    fig = go.Figure(data=data, layout=layout)
    prefix = "GeneTypeDistribution"
    draw(fig, prefix=prefix, outdir=outdir)
    return stat_data


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['detected_gene_type_stat'])



