# %%
import os
import matplotlib
# matplotlib.use('agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame
import scipy.stats as stats
import seaborn as sns
# import statsmodels.api as sm
# from statsmodels.stats.diagnostic import lillifors
matplotlib.rcParams['hatch.linewidth'] = 0.1
color_pool = plt.get_cmap('tab10').colors
# %matplotlib inline
# %config InlineBackend.figure_format = 'svg'


def bar(df, x_col='Sample', y_col='count', group_col='stage', group_col2='tumor_type', 
        ax=None, title=None, out=None, labelsize=6, labelrotation=90,
        ordered_group=None, ordered_group2=None, 
        hatch_lst=(None, '***', '///', '+++',  '...', 'xxx', 'ooo')):
    """
    :param df:
    :param x_col: column name for data source of bar label
    :param y_col: column name for data source of bar height
    :param group_col: column name for data source of bar color group. 柱子颜色分组信息
    :param group_col2: column name for data source of bar hatch(decorative pattern) group. 柱子里的花纹分组信息，可以赋值None
    :param ax: 
    :param title:
    :param out:
    :param labelsize: xtick label size
    :param labelrotation: xtick label rotation angle
    :param ordered_group: ordered group color, 决定柱子的顺序
    :param ordered_group2: ordered group hatch，与hatch_lst一一对应，决定哪些分组对应什么纹理
    :param hatch_lst: 柱子的花纹选项
    :return:
    """
    ax = plt.subplots()[1] if ax is None else ax
    # color_pool = plt.get_cmap('Set2').colors
    ordered_group = sorted(set(df[group_col])) if ordered_group is None else ordered_group
    color_dict = dict(zip(ordered_group, color_pool))
    colors = [color_dict[x] if x in color_dict else "gray" for x in df[group_col]]
    
    # 画柱子
    ax.grid(axis='y', linestyle='--', linewidth=0.5, alpha=0.5)
    bars = ax.bar(df[x_col], df[y_col], color=colors, align='center', width=0.7)
    
    # 设置柱子的面子颜色和纹理
    hatch_dict = dict()
    if group_col2 is not None:
        ordered_group2 = sorted(set(df[group_col2])) if ordered_group2 is None else ordered_group2
        hatch_dict = dict(zip(ordered_group2, hatch_lst))
        for ind, patch in enumerate(bars.patches):
            patch.set_hatch(hatch_dict[df[group_col2][ind]])
            patch.set_facecolor(color_dict[df[group_col][ind]])

    # customise legend
    group_count = df[group_col].value_counts()
    legend_patches = [mpatches.Patch(color=color_dict[x], label=x+f'({group_count[x]})') for x in ordered_group]
    if group_col2 is not None:
        group2_count = df[group_col2].value_counts()
        legend_patches += [mpatches.Patch(hatch=hatch_dict[x], label=x+f'({group2_count[x]})', facecolor='w', edgecolor='gray') for x in ordered_group2]
    ncol = 1 if len(legend_patches) <= 4 else 2
    ax.legend(handles=legend_patches, loc='upper right', ncol=ncol)

    # fine tune plot
    ax.xaxis.set_tick_params(labelsize=labelsize, labelrotation=labelrotation)
    if title:
        ax.set_title(title)
    ax.set_xlim(-1, df.shape[0])
    if out:
        plt.tight_layout()
        plt.savefig(f'{out}', dpi=300)
    return ax


def test_normality(data):
    # 正太分布检验
    if len(data) < 50:
        p_value = stats.normaltest(data)[1]
        return 'normalTest', p_value
    if len(data) < 300:
        p_value = stats.shapiro(data)[1]
        return "shapiro", p_value
    if len(data) >= 300:
        p_value = stats.kstest(data, 'norm')[1]
        return "kstest", p_value


def displot(d1, d2, label1, label2, ax=None, kind='kde'):
    s1 = pd.Series(d1)
    s2 = pd.Series(d2)
    # color_pool = plt.get_cmap('Set2').colors
    ax = plt.subplots()[1] if not ax else ax
    if kind == 'kde':
        ax = s1.plot.kde(label=label1, ax=ax, alpha=0.7, color=color_pool[0])
        ax = s2.plot.kde(label=label2, ax=ax, alpha=0.7, color=color_pool[1])
        ax.tick_params(labelsize='small')
        ax.legend(fontsize='small', loc='best')
        ax.set_title("Kernel density estimation")
    elif kind == 'hist':
        ax = s1.plot.hist(label=label1, ax=ax, alpha=0.7, color=color_pool[0])
        ax = s2.plot.hist(label=label2, ax=ax, alpha=0.6, color=color_pool[1])
        ax.yaxis.tick_right()
        ax.tick_params(labelsize='small')
        ax.legend(fontsize='small', loc='best')
        ax.set_title("Histogram plot")
    elif kind == 'box':
        def get_limit(series, k=3):
            q1 = series.quantile(q=0.25)
            q3 = series.quantile(q=0.75)
            iqr = q3 - q1
            upper = q3 + 3 * iqr
            lower = q1 - 3 * iqr
            return upper, lower
        # d1 = s1.clip(*get_limit(s1))
        # d2 = s2.clip(*get_limit(s1))
        d1, d2 = s1, s2
        df = pd.DataFrame({'value': list(d1)+list(d2), 'group': [label1]*len(d1) + [label2]*len(d2)})
        sns.boxplot(data=df, x='group', y='value', ax=ax, palette=color_pool)
        sns.swarmplot(data=df, x='group', y='value', ax=ax, color='0.25')
    else:
        raise Exception(f'unsupported kind {kind}')
    return ax


def pvalue_plot(d1, d2, label1, label2, axes=None):
    if axes is None:
        fig, axes = plt.subplots(2, 2, tight_layout=True)
        axes = [x for y in axes for x in y]

    def format_pvalue(x):
        return f'{x:.3e}' if x < 0.0001 else round(x, 4)

    # color_pool = plt.get_cmap('Set2').colors
    # 正态性检验： get pvalue
    method1, pvalue1 = test_normality(d1)
    method2, pvalue2 = test_normality(d2)

    # 正态性检验：Generates a probability plot
    for data, label, method, pvalue, ax, color in zip(
        [d1, d2], 
        [label1, label2],
        [method1, method2],
        [pvalue1, pvalue2],
        axes[2:], color_pool
        ):
        (osm, osr), (slope, intercept, r) = stats.probplot(data, dist="norm", plot=None)
        ax.plot(osm, osr, 'o', osm, slope*osm + intercept, 'r--', markerfacecolor=color)
        ax.annotate(
            f'{method}: pvalue={format_pvalue(pvalue)} and ' + "$R^2=%1.4f$" % r,
            xy=(0.05, 0.9), xycoords='axes fraction', fontsize=7
        )
        ax.set_title(f'Probability Plot of group {label}')

    # 方差齐性检验
    s, pvalue = stats.levene(d1, d2, center='median')

    # 两组均值差异：t检验
    t, t_test_pvalue = stats.ttest_ind(d1, d2, equal_var=(pvalue <= 0.05 or False))

    # 两组非参数检验
    s, rank_test_pvalue = stats.mannwhitneyu(d1, d2)

    # distribution plot：
    displot(d1, d2, label1, label2, ax=axes[0], kind='box')

    # plot table, 不使用rowLabels是为了避免空间占用，它占用的是axis外部的空间
    data = [
        ['T-Test', format_pvalue(t_test_pvalue)],
        ['Mann Whitney U Test', format_pvalue(rank_test_pvalue)],
        [f'"{label1}" Normality Test', format_pvalue(pvalue1)],
        [f'"{label2}" Normality Test', format_pvalue(pvalue2)],
    ]
    axes[1].set_title('P-value Table')
    table = axes[1].table(
        cellText=data, 
        colLabels=['Method', 'P-value'],
        cellLoc='center', loc='center', fontsize=8,
        colColours=['lightblue', 'lightblue'],
        colWidths=[0.62, 0.38]
    )
    axes[1].set_axis_off()
    # 调整table row height
    for i in [1, 2, 3, 4]:
        for j in [0, 1]:
            table[i, j].set_height(.2)
    return axes


def count_distribution(df, prefix='MHC-I_neoantigen_count', y_col='count', log2=True,
                       group_col='tumor_type', label1='adeno', label2='squamous'):
    # data process
    if log2:
        df[y_col] = np.log2(df[y_col])
    d1 = df[df[group_col] == label1][y_col]
    d2 = df[df[group_col] == label2][y_col]
    # layout setting
    fig = plt.figure(constrained_layout=True, figsize=(8, 9))
    gs = GridSpec(1+2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])
    stat_axes = []
    for i in [1, 2]:
        for j in [0, 1]:
            stat_axes.append(fig.add_subplot(gs[i, j]))

    # drawing
    bar(df, ax=ax1, x_col='Sample', y_col='count', group_col='stage', group_col2='tumor_type')
    pvalue_plot(d1, d2, label1=label1, label2=label2, axes=stat_axes)
    # save
    plt.savefig(f'{prefix}.stats.pdf', dpi=300)
    # plt.show()
    return fig


def cell_fraction_boxplot(file, hue='tumor_type'):
    # stats
    df = pd.read_csv(file)
    df['stage_group'] = ['I-II' if x in ['I', 'II'] else 'III-IV' for x in df['stage']]
    
    def get_data(x_col, hue, y_col):
        # group by x_col
        var_lst = sorted(set(df[x_col]))
        group_lst = sorted(set(df[hue]))
        box_data_dict = dict()
        d1 = df.groupby(x_col)
        for v in var_lst:
            d2 = d1.get_group(v).groupby(hue)
            for g in group_lst:
                if not g in d2.groups:
                    box_data = []
                else:
                    box_data = d2.get_group(g)[y_col].values
                box_data_dict.setdefault(v, dict())[g] = box_data
        return box_data_dict
    
    gd = get_data(x_col='cell_type', hue=hue, y_col='cell_fraction')
    pvalues = dict()
    man_pvalues = dict()
    for cell, value_dict in gd.items():
        print(value_dict.values())
        pvalues[cell] = stats.ttest_ind(*value_dict.values())
        man_pvalues[cell] = stats.mannwhitneyu(*value_dict.values())
    print(pvalues)
        
    # boxplot with table
    fig, ax = plt.subplots()
    sns.boxplot(data=df, x='cell_type', y='cell_fraction', hue_order=sorted(set(df[hue])),
                hue=hue, ax=ax, linewidth=1.5)
    ax.xaxis.set_tick_params(rotation=90)
    pvalues = [pvalues[x.get_text()][1] for x in ax.xaxis.get_majorticklabels()]
    man_pvalues = [man_pvalues[x.get_text()][1] for x in ax.xaxis.get_majorticklabels()]
    pvalue_lst = [f'{x:.3f}' for x in pvalues]
    man_pvalue_lst = [f'{x:.3f}' for x in man_pvalues]
    cell_text = pd.DataFrame({'ttest-pvalue': pvalue_lst, 'ranksum-pvalue': man_pvalue_lst}).T.values
    # print(cell_text)
    cmap = plt.get_cmap('RdYlBu')
    # 截取色谱中间的颜色
    mid_cmap = cmap(np.linspace(0.2, 0.8, 100))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('mycmap', mid_cmap)
    # cmap = plt.get_cmap('bwr')
    cell_colors = [[cmap(x) for x in pvalues], [cmap(x) for x in man_pvalues]]
    cell_colors = [[cmap(x) for x in pvalues], cmap(man_pvalues)]
    table = ax.table(
        cellText=cell_text,
        rowLabels=['ttest-pvalue', 'ranksum-pvalue'], 
        # colLabels=ax.get_xticklabels(), 
        cellColours=cell_colors,
        cellLoc='center', loc='top',
        # fontsize=8,
        alpha=0.5
    )
    # ax.set_xticklabels([])
    plt.subplots_adjust(left=0.2, top=0.8)
    
    # swarmplot
    fig, ax = plt.subplots(1, 1)
    sns.swarmplot(data=df, x='cell_type', y='cell_fraction', hue_order=sorted(set(df[hue])),
                  hue=hue, dodge=True, ax=ax, size=3)
    ax.xaxis.set_tick_params(rotation=90)
    
 
def cell_fraction_stackbar(file, group_field='tumor_type'):
    df = pd.read_csv(file)
    df['stage_group'] = ['I-II' if x in ['I', 'II'] else 'III-IV' for x in df['stage']]
    # stackbar
    fig, ax = plt.subplots(1, 1)
    df2 = df.pivot(columns='cell_type', index='Sample', values='cell_fraction')
    # get order of stack and get order of sample
    col_order = df2.iloc[0, :].sort_values(ascending=False).index
    df2['cellfraction_sum'] = df2.sum(axis=1)
    tumor_type_dict = dict(zip(df['Sample'], df[group_field]))
    stage_dict = dict(zip(df['Sample'], df['stage_group']))
    df2[group_field] = [tumor_type_dict[x] for x in df2.index]
    row_order = df2.sort_values(by=[group_field, 'cellfraction_sum'], ascending=False).index
    df2 = df2.loc[row_order, col_order]
    ax = df2.plot.bar(stacked=True, rot=90, ax=ax)
    first_legend = ax.legend(bbox_to_anchor=(1.02, 0.65), title='cell_type')
    ax.xaxis.set_tick_params(labelsize=5)

    # add tumor type group bar
    ymin, ymax = ax.get_ylim()
    group_bar_height = ymax*0.05
    labels = [tumor_type_dict[x] for x in row_order]
    color_dict = dict(zip(set(labels), plt.get_cmap('Accent').colors))
    colors = [color_dict[x] for x in labels]
    new_ymax = ymax+group_bar_height
    ax.set_ylim(ymin, new_ymax)
    bars = ax.bar(range(len(labels)), bottom=ymax, height=new_ymax, color=colors)
    assigned_labels = set()
    assigned_patches = []
    for label, patch in zip(labels, bars.patches):
        if label not in assigned_labels:
            patch.set_label(label)
            assigned_labels.add(label)
            assigned_patches.append(patch)
    ax.axhline(ymax, color='gray')
    second_legend = ax.legend(handles=assigned_patches, bbox_to_anchor=(1.02, 1), title=group_field)
    # 查看源码似乎不可能获得legend的宽度

    # stage info
    ymin, ymax = ax.get_ylim()
    labels = [stage_dict[x] for x in row_order]
    color_dict = dict(zip(set(labels), plt.get_cmap('Paired').colors))
    colors = [color_dict[x] for x in labels]
    new_ymax = ymax+group_bar_height
    ax.set_ylim(ymin, new_ymax)
    bars = ax.bar(range(len(labels)), bottom=ymax, height=new_ymax, color=colors)
    assigned_labels = set()
    assigned_patches = []
    for label, patch in zip(labels, bars.patches):
        if label not in assigned_labels:
            patch.set_label(label)
            assigned_labels.add(label)
            assigned_patches.append(patch)     
    third_legend = ax.legend(handles=assigned_patches, bbox_to_anchor=(1.02+0.5, 1), title='clinic_stage')
    # as first legend will be replaced by second one, we need to add it back
    ax.add_artist(first_legend)
    ax.add_artist(second_legend)
    # add grid
    ax.xaxis.grid(color='gray', linestyle='--', linewidth=0.2)
    return df2


def neocount_vs_cellfraction(count_file, cellfraction_file, hue=None, out='neocount_vs_cellfraction.pdf'):
    neo_count = pd.read_csv(count_file, index_col=0)
    cf = pd.read_csv(cellfraction_file)
    cf = cf.pivot(columns='cell_type', index='Sample', values='cell_fraction')
    df = neo_count.join(cf)
    df['log2Count'] = np.log2(df['count'])
    fig, axes = plt.subplots(cf.shape[1])
    for ind, cell_type in enumerate(cf.columns):
        if not hue:
            j = sns.jointplot(data=df[['log2Count', cell_type, 'tumor_type']], 
                              ax=axes[ind], x='log2Count', y=cell_type, kind='reg')
            r, p = stats.spearmanr(df['count'], df[cell_type])
            j.ax_joint.annotate(
                'r={:f}\nspearman_pval={:f}'.format(r,p),
                xy=(0.05, 0.9),
                xycoords='axes fraction'
            )
        else:
            j = sns.jointplot(data=df[['log2Count', cell_type, 'tumor_type']], ax=axes[ind],
                              x='log2Count', y=cell_type, hue=hue)
    fig.savefig(out, dpi=300, bbox_inches='tight')

    
if __name__ == '__main__':
    # for mhc_type, file in zip(
    #     ['MHC-I', 'MHC-II', 'MHC-both'],
    #     ['MHC_I.count.csv', 'MHC_II.count.csv', 'MHC_I_II.count.csv']
    #     ):
    #     df = pd.read_csv(file)
    #     df['stage_group'] = ['I-II' if x in ['I', 'II'] else 'III-IV' for x in df['stage']]
    #     count_distribution(df, prefix=f'{mhc_type}.adeno_vs_squamous', y_col='count', 
    #                        group_col='tumor_type', label1='adeno', label2='squamous')
    #     count_distribution(df, prefix=f'{mhc_type}.I-II_vs_III-IV', y_col='count', 
    #                        group_col='stage_group', label1='I-II', label2='III-IV')
    # cell_fraction_boxplot('./immuneCellFraction.csv')
    # cell_fraction_boxplot('./immuneCellFraction.csv', hue='stage_group')
    # cell_fraction_stackbar('./immuneCellFraction.csv')
    neocount_vs_cellfraction('./MHC_I.count.csv', './immuneCellFraction.csv')
    neocount_vs_cellfraction('./MHC_II.count.csv', './immuneCellFraction.csv', hue='tumor_type')
    

# %%
