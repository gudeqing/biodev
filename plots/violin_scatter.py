import os
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns


def violin(df, data_col, group_cols:list, hue_cols:list=None, index_col=None, split=False, orient=None,
                exchange_xy=False, out=None, scale='width', style='darkgrid', target_index=None, inner=None):
    sns.set(style=style)
    sns.set(font_scale=0.5)
    init_inner = inner
    # output name
    if out is None:
        if hue_cols is not None:
            out = '{data_col}.{group_cols}.{hue_cols}.violin.pdf'.format(
                data_col=data_col, group_cols='_'.join(group_cols), hue_cols='_'.join(hue_cols)
            )
        else:
            out = '{data_col}.{group_cols}.violin.pdf'.format(
                data_col=data_col, group_cols='_'.join(group_cols)
            )
    # prepare data
    if type(df) is str:
        data = pd.read_csv(df, index_col=0, header=0, sep=None, engine='python')
    else:
        data = df
    if index_col is not None:
        data.set_index(index_col, inplace=True)
    if target_index is not None:
        targets = [x.strip().split()[0] for x in open(target_index)]
        data = data.loc[targets]
    data.columns = [x.strip() for x in data.columns]
    # plot
    if hue_cols is None:
        hue_cols = [None] * len(group_cols)
    if len(group_cols) > 1:
        fig, axes = plt.subplots(len(group_cols), 1)
        for ind, group_col, hue_col in zip(range(len(group_cols)), group_cols, hue_cols):
            print(data.groupby(group_col).size())
            print(data.groupby(group_col).size().mean())
            if init_inner is None:
                if data.groupby(group_col).size().mean() > 50:
                    inner = 'quartile'
                else:
                    inner = 'stick'
            if hue_col is not None:
                if hue_col.lower() == 'none':
                    hue_col = None
            if not exchange_xy:
                ax = sns.violinplot(x=group_col, y=data_col, data=data, hue=hue_col, scale=scale, orient=orient,
                                    width=0.8, linewidth=0.5, inner=inner, split=split, ax=axes[ind])
            else:
                ax = sns.violinplot(x=data_col, y=group_col, data=data, hue=hue_col, scale=scale, orient=orient,
                                    width=0.8, linewidth=0.5, inner=inner, split=split, ax=axes[ind])
            plt.setp(ax.collections, linewidth=0.3)
    else:
        for ind, group_col, hue_col in zip(range(len(group_cols)), group_cols, hue_cols):
            if inner is None:
                if data.groupby(group_col).size().mean() > 50:
                    inner = 'quartile'
                else:
                    inner = 'stick'
            if hue_col is not None:
                if hue_col.lower() == 'none':
                    hue_col = None
            if not exchange_xy:
                ax = sns.violinplot(x=group_col, y=data_col, data=data, hue=hue_col, scale=scale, orient=orient,
                                    width=0.8, linewidth=0.5, inner=inner, split=split)
            else:
                ax = sns.violinplot(x=data_col, y=group_col, data=data, hue=hue_col, scale=scale, orient=orient,
                                    width=0.8, linewidth=0.5, inner=inner, split=split)
            plt.setp(ax.collections, linewidth=0.3)

    plt.savefig(out, dpi=300)


def scatter(df, data_col, group_cols:list, hue_cols:list=None, style_cols:list=None, index_col=None, exchange_xy=False,
            plot_style='whitegrid', target_index=None, out=None, legend='brief', alpha:float=None):
    sns.set_style(style=plot_style)
    sns.set(font_scale=0.5)
    # output name
    if out is None:
        out_name = 'scatter.' + data_col
        out_name += '.' + '-'.join(group_cols)
        if hue_cols is not None:
            out_name += '.' + '-'.join(hue_cols)
        if style_cols is not None:
            out_name += '.' + '-'.join(style_cols)
        out_name += '.png'
        print('Output:', out_name)
    else:
        out_name = out
    # prepare data
    if type(df) is str:
        data = pd.read_csv(df, index_col=0, header=0, sep=None, engine='python')
    else:
        data = df
    if index_col is not None:
        data.set_index(index_col, inplace=True)
    if target_index is not None:
        targets = [x.strip().split()[0] for x in open(target_index)]
        data = data.loc[targets]
    data.columns = [x.strip() for x in data.columns]
    # plot
    if hue_cols is None:
        hue_cols = [None] * len(group_cols)
    if style_cols is None:
        style_cols = [None] * len(group_cols)

    axes = None
    if len(group_cols) > 1:
        fig, axes = plt.subplots(len(group_cols), 1)
    if alpha is None:
        alpha='auto'

    for ind, group_col, hue_col, style_col in zip(range(len(group_cols)), group_cols, hue_cols, style_cols):
        print(data.groupby(group_col).size())
        print(data.groupby(group_col).size().mean())
        if hue_col is not None:
            if hue_col.lower() == 'none':
                hue_col = None
        if style_col is not None:
            if style_col.lower() == 'none':
                style_col = None
        if axes is not None:
            tmp_ax = axes[ind]
        else:
            tmp_ax = None
        if not exchange_xy:
            x_data = group_col
            y_data = data_col
        else:
            y_data = group_col
            x_data = data_col
        ax = sns.scatterplot(x=x_data, y=y_data, data=data, hue=hue_col, style=style_col, ax=tmp_ax,
                             legend=legend, alpha=alpha)
        plt.setp(ax.collections, linewidth=0.3)

    plt.savefig(out_name, dpi=300)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['violin', 'scatter'])
