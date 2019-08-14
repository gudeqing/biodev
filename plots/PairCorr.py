import os
from scipy import stats
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")
import pandas as pd
import numpy as np


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


def corr_annotate(x, y, method='pearson', font_size=9, x_pos=0.1, y_pos=0.9, **kwargs):
    func = correlation_function(method)
    r, pvalue = func(x, y)
    pvalue = format(pvalue, '.2e') if pvalue < 0.001 else round(pvalue, 4)
    ax = plt.gca()
    ax.annotate("{}.r={:.2f}\npvalue={}".format(method, r, pvalue), size=font_size,
                xy=(x_pos, y_pos), xycoords=ax.transAxes)


def pair_corr(data, hue=None, hue_order=None, palette=None, vars=None, top:int=None,
              x_vars=None, y_vars=None, diag_kind='hist',
              height=2.5, aspect=1, dropna=True, size=None,
              group=None, corr_method='pearson',
              log_base=2, log_additive=1,
              prefix='PairCorr', fig_format='png'):
    if type(data) == str:
        data = pd.read_csv(data, header=0, index_col=0, sep=None, engine='python')

    group_dict = dict()
    if type(group) == str:
        tmp_dict = dict(x.strip().split()[:2] for x in open(group))
        for k, v in tmp_dict.items():
            group_dict.setdefault(v, list())
            group_dict[v].append(k)
    elif type(group) == dict:
        group_dict = group
    elif group is None:
        pass
    else:
        raise Exception('value of group should be a group file or a python dict!')

    if not group_dict:
        group_dict = {'all': data.columns}

    for group_name, samples in group_dict.items():
        tdata = data[samples]
        tdata = tdata[tdata.apply(lambda x: sum(x > 0), axis=1) >= len(samples)]
        tdata.to_csv(f'{group_name}.comm.data.xls', sep='\t')
        if top is not None:
            mean_expr = tdata.mean(axis=1).sort_values(ascending=False)
            tdata = tdata.loc[mean_expr.index[:top]]
        if log_base == 2:
            tdata = np.log2(tdata+log_additive)
        elif log_base == 10:
            tdata = np.log10(tdata+log_additive)
        elif log_base in [0, 1]:
            print('no log transformation')
        else:
            raise Exception(f'log base {log_base} is not supported!')
        gridplot = sns.PairGrid(tdata, hue=hue, hue_order=hue_order, palette=palette,
                                vars=vars, x_vars=x_vars, y_vars=y_vars, height=height,
                                aspect=aspect, dropna=dropna, size=size)
        if diag_kind == 'hist':
            gridplot.map_diag(plt.hist)
        elif diag_kind == 'scatter':
            gridplot.map_diag(sns.scatterplot)
        elif diag_kind == 'kde':
            gridplot.map(sns.kdeplot)
        else:
            gridplot.map_diag(sns.scatterplot)
        for ax, col in zip(np.diag(gridplot.axes), tdata.columns):
            ax.set_xlabel(col)
        # lower
        gridplot.map_lower(sns.scatterplot)
        gridplot.map_lower(sns.regplot, scatter=False)
        gridplot.map_lower(corr_annotate, method=corr_method)

        # upper
        gridplot.map_upper(corr_annotate, method=corr_method, x_pos=0.15, y_pos=0.6, font_size=14)

        plt.savefig(prefix + f'.{group_name}.common{tdata.shape[0]}genes.{fig_format}', dpi=300)
        plt.close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pair_corr'])
