import os
from bokeh.models import ColumnDataSource, ColorBar, NumeralTickFormatter
from bokeh.models import Range1d
from bokeh.plotting import figure, save, output_file
from bokeh.layouts import gridplot
from bokeh.io import export_png
import bokeh.palettes as bp
from bokeh.transform import linear_cmap
import pandas as pd
import numpy as np
from colorpool import colorpool
from scipy import stats


def refExpScatter(exp_matrix, ref_samples:list, ref_mode='mean', ref_name='ref', id2symbol=None, corr_method='pearson',
                  log_base=2, log_additive=1, group_gene_files:list=None, scatter_ncols=2, mad_ncols=4,
                  prefix='refExpScatter', mad_ymax:float=None):
    """
    输入表达矩阵，指定参考样本，计算其他样本与其的相关性，并以scatter和MAD图的形式呈现
    :param exp_matrix: 表达矩阵
    :param ref_samples: 参考样列表
    :param ref_mode: 如mean表示计算参考样的平均值作为参考表达量
    :param ref_name: 参考样的名称
    :param id2symbol: 第一列为geneid而第二列为symbol的文件
    :param corr_method: 相关系数计算方法
    :param log_base: log base default 2
    :param log_additive: pseudo count added before log
    :param group_gene_files: 每个文件记录一组基因id, 方便对基因进行分组并标注不同的颜色
    :param scatter_ncols: 散点图的列数
    :param mad_ncols: MAD图的列数
    :param prefix: 输出文件prefix
    :param mad_ymax: MAD图的y坐标最大值
    :return:
    """
    exp_df = pd.read_csv(exp_matrix, index_col=0, header=0, sep=None, engine='python')
    if log_base == 2:
        exp_df = np.log2(exp_df + log_additive)
    elif log_base == 10:
        exp_df = np.log10(exp_df + log_additive)
    elif log_base == 1:
        pass
    else:
        raise Exception(f'{log_base} not supported! use 1(mean no log), 2, 10')
    id2symbol = dict(x.strip().split()[:2] for x in open(id2symbol)) if id2symbol else dict()
    if ref_mode == 'mean':
        ref_exp_list = exp_df.loc[:, ref_samples].mean(axis=1)
    else:
        ref_exp_list = exp_df.loc[:, ref_samples].median(axis=1)
    grouped_gene_list = list()
    names = ['gene']
    group_gene_files = [] if group_gene_files is None else group_gene_files
    for group_file in group_gene_files:
        names.append(os.path.basename(group_file))
        grouped_gene_list.append([x.strip().split()[0] for x in open(group_file)])
    if grouped_gene_list:
        not_grouped_genes = set(exp_df.index) - set(x for y in grouped_gene_list for x in y)
        grouped_gene_list.insert(0, not_grouped_genes)
    else:
        grouped_gene_list = [list(exp_df.index)]
    group_and_color = list(zip(grouped_gene_list, colorpool.get_color_pool(len(grouped_gene_list)), names))
    # plot scatter
    test_samples = [x for x in exp_df.columns if x not in ref_samples]
    exp_df['symbols'] = [id2symbol[x] if x in id2symbol else x for x in exp_df.index]
    plots = list()
    for each in test_samples:
        plot_data = exp_df.loc[:, [each]]
        plot_data[ref_name] = ref_exp_list
        plot_data = plot_data[plot_data[ref_name]*plot_data[each]>0]
        plot_data = plot_data.round(4)
        plot_data['symbols'] = exp_df['symbols']
        corr, pval = correlation_function(corr_method)(plot_data[each], plot_data[ref_name])
        pval = format(pval, '.2e') if pval < 0.001 else round(pval, 4)
        corr = round(corr, 3)
        p = figure(
            title="{} vs {}({}) {}_Corr={} Pval={}".format(each, ref_mode, ref_name,corr_method, corr, pval),
            # tools="wheel_zoom,reset,hover",
            tooltips=[
                ('x', '@{}'.format(ref_name)),
                ('y', '@{}'.format(each)),
                ('gene', '@symbols' if id2symbol else '@index'),
            ]
        )
        for group, color, name in group_and_color:
            target_index = list(set(group) & set(plot_data.index))
            source_data = plot_data.loc[target_index, :]
            if source_data.shape[0] == 0:
                continue
            source = ColumnDataSource(source_data)
            p.scatter(
                x=ref_name,
                y=each,
                # line_color=mapper,
                color=color,
                fill_alpha=0.2,
                size=5,
                legend=name,
                source=source
            )
        p.xaxis.axis_label = 'log{}(expr+{}) of {}'.format(log_base, log_additive, ref_name)
        p.yaxis.axis_label = 'log{}(expr+{}) of {}'.format(log_base, log_additive, each)
        plots.append(p)

    fig = gridplot(plots, ncols=scatter_ncols)
    output_file(prefix+'.scatter.html')
    save(fig)
    # export_png(fig, prefix+'.scatter.png')

    # plot Variation in
    # gene expression as a function of gene expression level across sample replicates
    upper_list = []
    plots = list()
    for each in test_samples:
        plot_data = exp_df.loc[:, [each]]
        plot_data[ref_name] = ref_exp_list
        plot_data = plot_data[plot_data[ref_name]*plot_data[each]>0]
        plot_data[each] = (plot_data[each] - plot_data[ref_name]).abs()/plot_data[ref_name]
        plot_data = plot_data.round(4)
        if mad_ymax is None:
            describe = plot_data[each].describe()
            upper = describe['75%'] + 2 * (describe['75%'] - describe['25%'])
        else:
            upper = mad_ymax
        upper_list.append(upper)
        plot_data[each][plot_data[each] > upper] = upper
        plot_data['symbols'] = exp_df['symbols']
        p = figure(
            title="{} vs {}({})".format(each, ref_mode, ref_name),
            # tools="wheel_zoom,reset,hover",
            tooltips=[
                ('x', '@{}'.format(ref_name)),
                ('y', '@{}'.format(each)),
                ('gene', '@symbols' if id2symbol else '@index'),
            ]
        )
        for group, color, name in group_and_color:
            source_data = plot_data.loc[set(group) & set(plot_data.index)]
            if source_data.shape[0] == 0:
                continue
            source = ColumnDataSource(source_data)
            p.scatter(
                x=ref_name,
                y=each,
                # line_color=mapper,
                color=color,
                fill_alpha=0.2,
                size=5,
                legend=name,
                source=source
            )
        p.xaxis.axis_label = 'log{}(expr+{}) of {}'.format(log_base, log_additive, ref_name)
        p.yaxis.axis_label = '|expr - {}_expr| / {}_expr'.format(ref_mode, ref_mode)
        plots.append(p)
        if mad_ymax is not None:
            p.y_range = Range1d(0, mad_ymax)
    else:
        fig = gridplot(plots, ncols=mad_ncols, sizing_mode='stretch_width')
        output_file(prefix+'.MDA.html')
        save(fig)


def correlation_function(func_name):
    func_name = func_name.lower()
    if func_name == "pearson":
        return stats.pearsonr
    elif func_name == "spearman":
        return stats.spearmanr
    elif func_name == "biserial":
        return stats.pointbiserialr
    elif func_name == 'kendall':
        return stats.kendalltau
    else:
        raise Exception('{} is not in [pearson, spearman, biserial, kendall]')

if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['refExpScatter'])
