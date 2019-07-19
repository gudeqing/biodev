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


def corr_annotate(x, y, method='pearson', **kwargs):
    func = correlation_function(method)
    r, pvalue = func(x, y)
    pvalue = format(pvalue, '.2e') if pvalue < 0.001 else round(pvalue, 4)
    ax = plt.gca()
    ax.annotate("{}.r={:.2f}\npvalue={}".format(method, r, pvalue), size=8,
                xy=(.1, .9), xycoords=ax.transAxes)


def pair_corr(data, hue=None, hue_order=None, palette=None, vars=None,
              x_vars=None, y_vars=None, kind='scatter', diag_kind='auto',
              markers=None, height=2.5, aspect=1, dropna=True, plot_kws=None,
              diag_kws=None, grid_kws=None, size=None,
              group=None, corr_method='pearson', log_base=2, log_additive=1,
              prefix='PairCorr', fig_format='pdf'):
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
    else:
        raise Exception('value of group should be a group file or a python dict!')

    if not group_dict:
        group_dict = {'all': data.columns}

    for group_name, samples in group_dict.items():
        tdata = data[samples]
        tdata = tdata[tdata.apply(lambda x: sum(x > 0), axis=1) >= len(samples)]
        if log_base == 2:
            tdata = np.log2(tdata+log_additive)
        elif log_base == 10:
            tdata = np.log10(tdata+log_additive)
        elif log_base in [0, 1]:
            print('no log transformation')
        else:
            raise Exception(f'log base {log_base} is not supported!')
        gridplot = sns.pairplot(tdata, hue=hue, hue_order=hue_order, palette=palette,
                                vars=vars, x_vars=x_vars, y_vars=y_vars, kind=kind,
                                diag_kind=diag_kind, markers=markers, height=height,
                                aspect=aspect, dropna=dropna, plot_kws=plot_kws,
                                diag_kws=diag_kws, grid_kws=grid_kws, size=size)
        gridplot.map_upper(sns.regplot, scatter=False)
        gridplot.map_upper(corr_annotate, method=corr_method)
        plt.savefig(prefix + f'.{group_name}.common{tdata.shape[0]}genes.{fig_format}')
        plt.close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pair_corr'])
