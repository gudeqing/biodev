import os
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gaussian_kde

d2 = pd.DataFrame()

fig, ax = plt.subplots()
fig.set_size_inches(w=3, h=12)
ax.grid(color='gray', linestyle='--', linewidth=0.5)
ax.set_xlim(0.4, 1.02)
ax.yaxis.tick_right()
ax.xaxis.tick_top()
outlier_style = dict(markerfacecolor='blue', markeredgecolor='blue', marker='D', markersize=1)
ax.boxplot(d2, vert=False, flierprops=outlier_style)
# ax.set_xticklabels(ax.get_xmajorticklabels(), fontdict={'fontsize':5})
ax.set_yticklabels(d2.index, fontdict={'fontsize':5})
ax.tick_params('x', labelsize=7)
plt.savefig(f'boxplot.pdf', bbox_inches='tight')
plt.close()

# ax = axes[i].twinx()
# density = gaussian_kde(s1)
# xs = range(s1.min()-2, s1.max()+2)
# ax.plot(xs, density(xs), label='Normal', alpha=0.7, color='darkgreen')
# density = gaussian_kde(s2)
# xs = range(s2.min() - 2, s2.max() + 2, )
# ax.plot(xs, density(xs), label='Tumor', alpha=0.7, color='darkorange')
fig, ax = plt.subplots()
data = np.random.normal(size=100)
ax.hist(data)
from scipy.stats import gaussian_kde
s1 = pd.Series(data)
density = gaussian_kde(s1)
xs = np.linspace(s1.min()-2, s1.max()+2, 100)
ax2 = ax.twinx()
ax2.plot(xs, density(xs), label='Normal', alpha=0.7, color='darkgreen')
plt.savefig(f'hist_density.pdf', bbox_inches='tight')


def hist_density(data, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    ax.hist(data)
    s1 = pd.Series(data)
    density = gaussian_kde(s1)
    xs = np.linspace(s1.min() - 2, s1.max() + 2, 100)
    ax2 = ax.twinx()
    ax2.plot(xs, density(xs), label='Normal', alpha=0.7, color='darkgreen')
