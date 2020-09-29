import os
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

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


