import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")
import pandas as pd
import itertools
from scipy.stats import gmean
# plt.style.use('ggplot')


def _plot_multiline(data, out='multiline.png', annotate_at_end=False):
    colors = ['b', 'r', 'g', 'c', 'm', 'y', 'grey', 'pink', 'darkgrey']
    line_style = ['-', '-.', ':']
    cmbs = list(itertools.product(colors, line_style))
    plt.figure(figsize=(7, 8))
    for ind, index in enumerate(data.index):
        plt.plot(
            range(data.shape[1]),
            data.loc[index],
            label=index,
            color=cmbs[ind][0],
            linestyle=cmbs[ind][1],
            marker='o',
            markeredgecolor=cmbs[ind][0],
            markerfacecolor='w',
            markersize=5
        )
        if annotate_at_end:
            plt.annotate(
                index,
                xy=(data.shape[1]-1, data.loc[index][-1]),
                fontsize=5
            )
    plt.xticks(range(data.shape[1]), data.columns)
    plt.xlim(-0.2, data.shape[1])
    if not annotate_at_end:
        plt.legend(loc=1, fontsize='xx-small')
    plt.savefig(out, dpi=300, bbox_inches='tight')


def plot_exp_lines(data, index_col:list=0, sample_group=None, annotate_at_end=False,
                   group_order:list=None, out='multiline.png'):
    if sample_group is not None:
        group_df = pd.read_csv(sample_group, index_col=0, header=0, sep=None, engine='python')
        group_names = set(group_df.iloc[:, 0])
        if group_order is None:
            group_names = sorted(list(group_names))
        else:
            group_names = group_order
        if type(data) == str:
            data = pd.read_csv(data, index_col=index_col, header=0, sep=None, engine='python')
        centered = data.sub(data.mean(axis=1), axis=0)
        mean_centered = pd.DataFrame()
        for name in group_names:
            target_samples = list(group_df.loc[group_df.iloc[:, 0]==name].index)
            mean_centered[name] = centered[target_samples].mean(axis=1)
        # print(geometric.mean(axis=1))
        plot_data = mean_centered
    else:
        if type(data) == str:
            plot_data = pd.read_csv(data, index_col=index_col, header=0, sep=None, engine='python')
        else:
            plot_data = data

    _plot_multiline(plot_data, out=out, annotate_at_end=annotate_at_end)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['plot_exp_lines'])

