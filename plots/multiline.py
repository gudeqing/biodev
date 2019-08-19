import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")
import pandas as pd
import numpy as np
import itertools
from scipy.stats import zscore
from scipy.stats import gmean
# plt.style.use('ggplot')


def plot(data, out='multiline.png'):
    colors = ['b', 'r', 'g', 'c', 'm', 'y', 'grey', 'pink', 'darkgrey']
    line_style = ['-', '-.', ':']
    cmbs = list(itertools.product(colors, line_style))
    # max_data = data.values.max()
    # min_data = data.values.min()
    # print(max_data)
    # step = (max_data-min_data)/(data.shape[0]+2)
    # print(step)
    plt.figure(figsize=(7, 8))
    for ind, index in enumerate(data.index):
        plt.plot([-1, 1], data.loc[index], label=index, color=cmbs[ind][0], linestyle=cmbs[ind][1],
                 marker='o', markeredgecolor=cmbs[ind][0], markerfacecolor='w',
                 markersize=5)
        # plt.annotate(index, xy=(x_coor[-1], data.loc[index][-1]), xytext=(x_coor[-1]+0.1, max_data-step*ind),
        #              fontsize=5, arrowprops=dict(arrowstyle='->', lw=0.3))
    plt.xticks([-1, 1], data.columns)
    plt.xlim(-1.1, 1+1.5)
    plt.legend(loc=1, fontsize='xx-small')
    plt.savefig(out, dpi=300, bbox_inches='tight')


def multiline_geometric_mean(data, sample_group, out='multiline.png', group_order=None):
    group_df = pd.read_csv(sample_group, index_col=0, header=0, sep=None, engine='python')
    group_names = set(group_df.iloc[:, 0])
    group_names = sorted(list(group_names))
    data = pd.read_csv(data, index_col=0, header=0, sep=None, engine='python')
    geometric = pd.DataFrame()
    for name in group_names:
        target_samples = list(group_df.loc[group_df.iloc[:, 0]==name].index)
        tmp = data[target_samples]
        tmp = tmp.transpose().apply(gmean).transpose()
        geometric[name] = tmp
    # print(geometric.mean(axis=1))
    centered = geometric.sub(geometric.mean(axis=1), axis=0)
    # print(centered)
    plot(centered, out=out)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['multiline_geometric_mean'])
