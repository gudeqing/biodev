import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns; sns.set()
import pandas as pd
sns.set(font_scale=0.5)


def heatmap(data, vmin:float=None, vmax:float=None, cmap:str='RdYlBu_r', center:float=None, robust=False, label_size=8,
            pvalue:str=None, annot=True, fmt='', linewidths:float=0, linecolor='white', cbar=True,
            annot_fontsize=6, square=False, xticklabels='auto', yticklabels='auto', out='heatmap.png', dpi=300):
    data = pd.read_csv(data, header=0, index_col=0, sep=None, engine='python')
    if annot:
        annot = data.round(3).applymap(str)
        if pvalue is not None:
            pvalues = pd.read_csv(pvalue, header=0, index_col=0, sep=None, engine='python')
            annot += '(p=' + pvalues.round(3).applymap(str) + ')'
    ax = sns.heatmap(data, vmin=vmin, vmax=vmax, cmap=cmap, center=center, robust=robust, annot_kws={"size": annot_fontsize},
                    annot=annot, fmt=fmt, linewidths=linewidths, linecolor=linecolor,
                    cbar=cbar, square=square, xticklabels=xticklabels, yticklabels=yticklabels)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=label_size)
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()


# heatmap.__doc__ = sns.heatmap.__doc__


if __name__ == '__main__':
    from xcmds.xcmds import xcmds
    xcmds(locals(), include=['heatmap'])

