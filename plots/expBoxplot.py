import os
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np


def expr_box_plot(expr_matrix, sample_group=None, x_col='gene', xlabel=None, prefix='',
                  log_transform=False, log_additive=1):
    data = pd.read_csv(expr_matrix, header=0, index_col=0, sep=None, engine='python')
    if log_transform:
        data = np.log2(data + log_additive)
    data.index.name = 'Gene'
    samples = data.columns
    data = data.reset_index('Gene')
    data = data.melt(id_vars=['Gene'], var_name='Sample', value_name='Expression')
    if sample_group:
        group = pd.read_csv(sample_group, header=0, index_col=0, sep=None, engine='python')
    else:
        group = pd.DataFrame({'All': {k:k for k in samples}})
    group.index.name = 'Sample'
    group_names = group.columns
    group.reset_index('Sample', inplace=True)
    data = data.merge(group, on='Sample')
    for name in group_names:
        hue = None if sample_group is None else name
        if x_col.lower() == 'gene':
            ax = sns.boxplot(x='Gene', y='Expression', hue=hue, data=data)
        else:
            ax = sns.boxplot(x='Sample', y='Expression', hue=hue, data=data)
        ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=6, rotation=90)
        if xlabel:
            ax.set(xlabel=xlabel)
        plt.savefig(f'{prefix}{name}.boxplot.png', dpi=300, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['expr_box_plot'])
