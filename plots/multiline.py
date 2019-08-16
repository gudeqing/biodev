import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")
import pandas as pd
import numpy as np
import itertools
from scipy.stats import zscore
# plt.style.use('ggplot')


def plot(data, out='multiline.png'):
    colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
    line_style = ['-', '--', '-.', ':']
    cmbs = list(itertools.product(colors, line_style))
    max_data = data.values.max()
    min_data = data.values.min()
    print(max_data)
    step = (max_data-min_data)/(data.shape[0]+2)
    print(step)
    x_coor = [x + 0.1 for x in range(data.shape[1])]
    for ind, index in enumerate(data.index):
        plt.plot(x_coor, data.loc[index], label=index, color=cmbs[ind][0], linestyle=cmbs[ind][1])
        # plt.annotate(index, xy=(x_coor[-1], data.loc[index][-1]), xytext=(x_coor[-1]+0.1, max_data-step*ind),
        #              fontsize=5, arrowprops=dict(arrowstyle='->', lw=0.3))
    plt.xticks(x_coor, data.columns)
    plt.xlim(0, x_coor[-1]+0.5)
    plt.legend(loc=1, fontsize='xx-small')
    plt.savefig(out, dpi=300)


def read_data(data, target_cols:list, out='multiline.png'):
    data = pd.read_csv(data, index_col=0, header=0, sep=None, engine='python')
    data = data[target_cols]
    print(data.head())
    plot(data, out=out)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['read_data'])
