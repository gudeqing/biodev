import os
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns


def expr_box_plot(expr_matrix, sample_group, x_col='gene', xlabel=None, prefix=''):
    data = pd.read_csv(expr_matrix, header=0, index_col=0, sep=None, engine='python')
    data.index.name = 'Gene'
    data = data.reset_index('Gene')
    data = data.melt(id_vars=['Gene'], var_name='Sample', value_name='Expression')
    group = pd.read_csv(sample_group, header=0, index_col=0, sep=None, engine='python')
    group.index.name = 'Sample'
    group_names = group.columns
    group.reset_index('Sample', inplace=True)
    data = data.merge(group, on='Sample')
    print(data.head())
    for name in group_names:
        if x_col.lower() == 'gene':
            ax = sns.boxplot(x='Gene', y='Expression', hue=name, data=data)
        else:
            ax = sns.boxplot(x='Sample', y='Expression', hue=name, data=data)
        ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=6, rotation=90)
        if xlabel:
            ax.set(xlabel=xlabel)
        plt.savefig(f'{prefix}{name}.boxplot.png', dpi=300, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['expr_box_plot'])
