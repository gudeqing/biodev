import numpy as np
from matplotlib import pyplot as plt


def volcano_plot(data, pvalue_cutoff=0.05, fc_cutoff=2.0, xlimit=(-3, 3), ylimit=(-0.5, 5), fig_name='volcano plot',
                 fig_size=(7, 8), fig_title='', dpi=300, delimiter='\t', markers=('v', '^', '.'), mark_points=None,
                 xlabel='log2(Fold Change)', ylabel='-log10(pvalue)'):
    """
    :param data: array or text file without column header, second column contains p-values, and the
    other column contains log2 of fold change,negative value denotes down-regulation.
    :param pvalue_cutoff: the result of ttest, default value is 0.05.
    :param fc_cutoff: fold change cutoff, default value is 0.05
    :return: a figure, plotting the negative log of the p-value on the y-axis (usually base 10);
    and the x-axis is the log of the fold change between the two conditions. negative value represents
    down-regulation, and positive value represents up-regulation.
    """
    if data not in locals().keys():
        print('>>> I will get data from:', data)
        data = np.genfromtxt(data, skip_header=1, comments=None, dtype=float, delimiter=delimiter)
    else:
        data = np.array(data, dtype=float)
    # distribute data
    fc_cutoff = np.log2(fc_cutoff)
    p_values = -np.log10(data[:, 1])
    p_values[p_values > ylimit[1]] = ylimit[1]*0.99  # to make points visible in fig
    fc_values = data[:, 0]
    fc_values[fc_values > xlimit[1]] = xlimit[1]*0.99
    fc_values[fc_values < xlimit[0]] = xlimit[0]*0.99
    point_number = len(p_values)
    # scatter plot
    plt.figure(fig_name, fig_size)
    plt.subplot(111, axisbg='w')
    colors = np.tile('#CDC5BF', point_number)
    colors[np.logical_and(p_values >= -np.log10(pvalue_cutoff), fc_values >= fc_cutoff)] = 'r'
    colors[np.logical_and(p_values >= -np.log10(pvalue_cutoff), fc_values <= -fc_cutoff)] = 'g'
    # plot
    plt.scatter(fc_values[colors == 'r'], p_values[colors == 'r'], marker=markers[1], c=colors[colors == 'r'], alpha=0.6, edgecolors='r', label='up')
    plt.scatter(fc_values[colors == 'g'], p_values[colors == 'g'], marker=markers[0], c=colors[colors == 'g'], alpha=0.7, edgecolors='g', label='down')
    plt.scatter(fc_values[colors=='#CDC5BF'], p_values[colors=='#CDC5BF'], marker=markers[2], c=colors[colors=='#CDC5BF'], alpha=0.5, edgecolors=None, label='--')
    # mark points
    if mark_points:
        for point in mark_points:
            plt.annotate(point[0], xy=(point[1],point[2]))
    # set limit axis of x and y
    axes = plt.gca()
    axes.set_xlim(xlimit)
    axes.set_ylim(ylimit)
    # p-value cutoff line
    plt.axhline(- np.log10(pvalue_cutoff), c='b', lw=1, ls='--')
    # right fold change cutoff line
    plt.axvline(fc_cutoff, c='b', lw=1, ls='--')
    # left fold change cutoff line
    plt.axvline(-fc_cutoff, c='b', lw=1, ls='--')
    # label, title...
    plt.xlabel("$"+xlabel+"$")
    plt.ylabel("$"+ylabel+"$")
    plt.legend()
    plt.title(fig_title, {'fontsize':16})
    plt.grid()
    plt.tight_layout()
    plt.savefig(fig_name, dpi=dpi)
    plt.close('all')

if __name__ == "__main__":
    #volcano_plot('log2fc_pvalue.txt', pvalue_cutoff=0.05, fc_cutoff=1.5, fig_size=(7, 8), fig_name='volcano_plot',
    #             fig_title='S vs B', dpi=300, delimiter='\t', mark_points=[('test', 2, 3),])

    volcano_plot('log2fc_fdr.txt', pvalue_cutoff=0.05, fc_cutoff=1.5, fig_size=(7, 8), fig_name='volcano_plot',
            fig_title='S vs B', dpi=300, delimiter='\t', ylabel='log10(fdr)')

