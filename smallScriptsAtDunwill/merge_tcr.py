import pandas as pd
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns


def merge_vj_matrix(file_list:list, column_sep='D', out='merged.counts.csv', group_info:str=None):
    table = pd.read_csv(file_list[0], index_col=[0, 1, 2, 3], header=0, sep=None, engine='python')
    for each in file_list[1:]:
        each_table = pd.read_csv(each, index_col=[0, 1, 2, 3], header=0)
        table = table.join(each_table, how='outer')
    table.columns = [x.strip().split(column_sep, 1)[0] for x in table.columns]
    table = table.fillna(0)
    table.index.name = 'SampleID'
    group_df = pd.read_csv(group_info, index_col=0, header=0, sep=None, engine='python')
    table = group_df.transpose().append(table, sort=False)
    # write out result
    table.to_csv(out)
    return table


def merge_metric_matrix(file_list:list, column_sep='D', out='merged.metric.csv', group_info:str=None):
    table = pd.read_csv(file_list[0], index_col=None, header=0, sep=None, engine='python')
    for each in file_list[1:]:
        each_table = pd.read_csv(each, index_col=0, header=0)
        table = table.append(each_table, sort=False)
    samples = [x.strip().split(column_sep, 1)[0] for x in table['Sample_Name']]
    table = table.iloc[:, 4:]
    table.index = samples
    table.index.name = 'SampleID'
    group_df = pd.read_csv(group_info, index_col=0, header=0, sep=None, engine='python')
    table = group_df.join(table, how='right', sort=False)
    table.columns = [x.strip() for x in table.columns]
    table.to_csv(out)
    return table


def violin_plot(df, data_col, group_cols:list, hue_cols:list=None, index_col=None,
                out=None, scale='width', style='darkgrid', target_index=None, inner=None):
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
            ax = sns.violinplot(x=group_col, y=data_col, data=data, hue=hue_col, scale=scale,
                                width=0.8, linewidth=0.5, inner=inner, split=False, ax=axes[ind])
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
            ax = sns.violinplot(x=group_col, y=data_col, data=data, hue=hue_col, scale=scale,
                                width=0.8, linewidth=0.5, inner=inner, split=False)
            plt.setp(ax.collections, linewidth=0.3)

    plt.savefig(out, dpi=300)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), exclude=['pd'])

