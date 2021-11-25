# %%
import os
import matplotlib
# matplotlib.use('agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
# import statsmodels.api as sm
matplotlib.rcParams['hatch.linewidth'] = 0.1

"""
基于matplotlib实现分组boxplot，并和seaborn.boxplot对比 
"""


def pair_boxplot(df: pd.DataFrame, x_col, y_col, hue, ax=plt.subplot()):
    def get_limit(series, k=3):
        """get cutoff for clip of outliers"""
        q1 = series.quantile(q=0.25)
        q3 = series.quantile(q=0.75)
        iqr = q3 - q1
        upper = q3 + 3 * iqr
        lower = q1 - 3 * iqr
        return upper, lower

    def get_ticks(major_tick_num, group_num):
        major_tick = range(1, major_tick_num+1)
        minor_ticks_list = []
        interval_width = 0.1
        box_width = (1-interval_width*2)/group_num
        for i in range(group_num):
            minor_ticks_list += [[(m-0.5+interval_width) +
                                  box_width*i for m in major_tick]]
        return major_tick, minor_ticks_list, box_width

    def get_data(var_lst, group_lst):
        # group by x_col
        box_data_dict = dict()
        d1 = df.groupby(x_col)
        for v in var_lst:
            d2 = d1.get_group(v).groupby(hue)
            for g in group_lst:
                if not g in d2.groups:
                    box_data = []
                else:
                    box_data = d2.get_group(g)[y_col].values
                box_data_dict.setdefault(g, []).append(box_data)
        return box_data_dict

    def get_colors(n):
        return plt.get_cmap('tab10').colors[:n]

    lengend_patches = []
    var_num = len(set(df[x_col]))
    group_num = len(set(df[hue]))
    var_lst = sorted(set(df[x_col]))
    group_lst = sorted(set(df[hue]))
    major_ticks, minor_ticks_lst, box_width = get_ticks(var_num, group_num)
    box_data_dict = get_data(var_lst, group_lst)
    colors = get_colors(group_num)
    for pos, group, color in zip(minor_ticks_lst, box_data_dict, colors):
        res = ax.boxplot(
            box_data_dict[group],
            positions=pos,
            widths=box_width,
            patch_artist=True
        )
        lengend_patches.append(res['boxes'][0])
        for box in res['boxes']:
            box.set_fc(color)
        for median_line in res['medians']:
            median_line.set_color('k')
    ax.legend(lengend_patches, group_lst, loc='best')
    ax.set_xticks(major_ticks)
    ax.set_xlim(minor_ticks_lst[0][0]-box_width - 0.005, 
                minor_ticks_lst[-1][-1]+box_width + 0.005)
    ax.set_xticklabels(var_lst, rotation=0)


if __name__ == '__main__':
    fig, axes = plt.subplots(2, 1)
    df = pd.DataFrame(np.random.rand(10, 5), columns=['a', 'b', 'c', 'd', 'e'])
    df2 = df.melt()
    df2['g'] = list('AABCABCAAC'*5)
    df2['k'] = list('XXXYYYZZZX'*5)
    pair_boxplot(df2, x_col='g', y_col='value', hue='k', ax=axes[0])
    sns.boxplot(x='g', y='value', hue='k', data=df2, ax=axes[1])


# %%
