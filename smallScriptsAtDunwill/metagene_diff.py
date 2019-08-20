import os
import pandas as pd
import scipy.stats as stats
from functools import partial
import statistics
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")
import itertools
from scipy.stats import gmean


def read_metagene_group(group_txt, format='gene2group'):
    if format != 'gene2group':
        result_dict = dict()
        with open(group_txt) as f:
            for line in f:
                metagene, genes = line.strip().split('\t', 1)
                genes = [y.strip() for x in genes.split() for y in x.strip().split(',')]
                if metagene in result_dict:
                    raise Exception(f"{metagene} duplicated, please check it!")
                else:
                    result_dict[metagene] = genes
    else:
        result_dict = read_sample_group(group_txt)
    return result_dict


def read_sample_group(group_txt):
    group_df = pd.read_csv(group_txt, sep=None, engine='python', header=0, index_col=0)
    group_dict = dict()
    for scheme in group_df:
        tmp_dict = dict(list(group_df.loc[:, [scheme]].groupby(scheme)))
        for group, df_val in tmp_dict.items():
            if df_val.shape[0] == group_df.shape[0]:
                raise Exception('In column of {}, only one group was found!'.format(scheme))
            group_dict[group] = sorted(df_val.index)
    return group_dict


def read_compare_info(cmp_info, group_dict):
    with open(cmp_info) as f:
        cmp_list = list()
        error_names = list()
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            tmp_ctrl, tmp_test = line.strip().split()
            if tmp_ctrl not in group_dict:
                error_names.append(tmp_ctrl)
            if tmp_test not in group_dict:
                error_names.append(tmp_test)
            cmp_list.append((tmp_ctrl, tmp_test))
        if error_names:
            print(f'Please be aware that Each group name of {error_names} is not in group dict!')
    cmp_list = sorted(list(set(cmp_list)))
    return cmp_list


def score_metagene(exp_matrix, metagene_dict:dict, score='geometric_mean'):
    exp_df = pd.read_csv(exp_matrix, header=0, index_col=0, sep=None, engine='python')
    score_df = pd.DataFrame(columns=exp_df.columns)
    valid_member_list = list()
    for metagene, genes in metagene_dict.items():
        intersect = set(genes) & set(exp_df.index)
        if intersect:
            valid_members = ','.join(intersect)
            valid_member_list.append(valid_members)
            meta_exp = exp_df.loc[intersect]
            if score == 'geometric_mean':
                meta_score = meta_exp.apply(stats.gmean, axis=0)
            else:
                meta_score = meta_exp.mean(axis=0)
            score_df.loc[metagene] = meta_score
        else:
            print(f'None member of metagene "{metagene}" is in expression matrix')
    score_df['members'] = valid_member_list
    score_df.set_index('members', append=True, inplace=True)
    return score_df


def diff_test(score_df, group_dict, cmp_list, method='mannwhitneyu', equal_var=True, prefix=''):
    for ctrl, test in cmp_list:
        if not(ctrl in group_dict and test in group_dict):
            print(f'skip {ctrl} and {test} for both or one of them are/is not in group info dict!')
            continue
        ctrl_samples = group_dict[ctrl]
        test_samples = group_dict[test]
        print(ctrl_samples)
        print(test_samples)
        ctrl_num = len(ctrl_samples)
        target_data = score_df[ctrl_samples+test_samples]
        centered = target_data.sub(target_data.mean(axis=1), axis=0)
        mean_centered = pd.DataFrame()
        mean_centered[test] = centered[test_samples].mean(axis=1)
        mean_centered[ctrl] = centered[ctrl_samples].mean(axis=1)

        if method == 'ranksums':
            test_func = stats.ranksums
        elif method == 'mannwhitneyu':
            test_func = stats.mannwhitneyu
        elif method == "wilcoxon":
            test_func = stats.wilcoxon
        elif method == 'ttest_ind':
            if equal_var:
                test_func = stats.ttest_ind
            else:
                test_func = partial(stats.ttest_ind, equal_var=False)
        else:
            raise Exception(f'{method} is not supported! '
                            f'Choose one of [ranksums, mannwhitneyu, wilcoxon, ttest_ind')
        test_df = pd.DataFrame()
        test_df['pvalue'] = target_data.apply(lambda x:test_func(x[:ctrl_num], x[ctrl_num:])[1], axis=1)
        ctrl_median_exp = target_data[ctrl_samples].apply(statistics.median, axis=1)
        test_median_exp = target_data[test_samples].apply(statistics.median, axis=1)
        test_df[ctrl+'_median'] = ctrl_median_exp
        test_df[test+'_median'] = test_median_exp
        test_df['median_log2FC'] = test_median_exp - ctrl_median_exp
        test_df['regulation'] = 'down'
        test_df.loc[test_df['median_log2FC']>0, 'regulation'] = 'up'
        test_df = test_df.join(target_data)
        test_df.sort_values(by='pvalue', inplace=True)
        test_df.to_csv(f'{prefix}{ctrl}_vs_{test}.metagene.{method}.xls', sep='\t')
        # plot
        mean_centered = mean_centered.loc[test_df.index]
        mean_centered.index = [x[0] for x in mean_centered.index]
        plot_exp_lines(mean_centered, out=f'{prefix}{ctrl}_vs_{test}.metagene.mean_centered.png')


def metagene_diff(exp_matrix, sample_group, metagene_group, compare, prefix='',
                  metagene_group_format='gene2group', score='geometric_mean',
                  method='mannwhitneyu', equal_var=True, box_x='gene',
                  box_xlabel='MetaGene'):
    """
    给定表达矩阵, 样本分组信息和基因分组信息(metagene) 以及比较信息, 对meta基因做差异检验
    :param exp_matrix: 表达矩阵文件路径
    :param sample_group: 样本分组信息路径, 第一行为header
    :param metagene_group: 基因分组信息,也即metagene信息文件路径. 文件无需header, 内容有两种格式可选:
       (1) 第一列为metagene名,第二列为metagene的members，用逗号隔开 (2)第一列为基因名，第二列为metagene名
    :param compare: 分组比较信息，第一列为对照组，第二列为实验组/测试组
    :param prefix: 输出文件的前缀
    :param metagene_group_format: gene2group表示第二种格式, group2gene表示第三种格式
    :param score: 计算metagene score的方式, 默认geometric_mean, 即计算members的几何平均值，否则为算术平均值
    :param method: 检验方法, 默认'mannwhitneyu', 还有[ranksums, wilcoxon(配对比较时可选用), ttest_ind]可选
    :param equal_var: method为ttest_ind时可以指定的参数,默认做假设等方差检验
    :param box_x: 画box图的参数 横轴以基因为单位(默认), 另可指定为sample, 此时横轴以sample为单位
    :param box_xlabel: 画box图的参数 横轴的label, 默认为MetaGene
    :return:
    """
    metagene_group_dict = read_metagene_group(metagene_group, format=metagene_group_format)
    metagene_score_df = score_metagene(exp_matrix, metagene_group_dict, score=score)
    box_data = metagene_score_df.copy()
    box_data.index = [x for x, _ in box_data.index]
    expr_box_plot(box_data, sample_group, x_col=box_x, xlabel=box_xlabel, prefix=prefix)
    sample_group_dict = read_sample_group(sample_group)
    compare_list = read_compare_info(compare, sample_group_dict)
    print(compare_list)
    diff_test(metagene_score_df, sample_group_dict, compare_list,
              method=method, equal_var=equal_var, prefix=prefix)


def _plot_multiline(data, out='multiline.png', annotate_at_end=False):
    colors = ['b', 'r', 'g',  'y', 'c', 'pink', 'm', 'grey', 'darkgrey']
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


def expr_box_plot(expr_matrix, sample_group, x_col='gene', xlabel=None, prefix=''):
    if type(expr_matrix) == str:
        data = pd.read_csv(expr_matrix, header=0, index_col=0, sep=None, engine='python')
    else:
        data = expr_matrix
    data.index.name = 'Gene'
    data = data.reset_index('Gene')
    data = data.melt(id_vars=['Gene'], var_name='Sample', value_name='Expression')
    group = pd.read_csv(sample_group, header=0, index_col=0, sep=None, engine='python')
    group.index.name = 'Sample'
    group_names = group.columns
    group.reset_index('Sample', inplace=True)
    data = data.merge(group, on='Sample')
    # print(data.head())
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
    xcmds.xcmds(locals(), include=['metagene_diff'])

