import os
from bokeh.models import ColumnDataSource, ColorBar, NumeralTickFormatter
from bokeh.plotting import figure, save, output_file
from bokeh.layouts import gridplot
import bokeh.palettes as bp
from bokeh.transform import linear_cmap
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import venn
from upsetplot import from_contents
from upsetplot import plot



def plotTopExpGenes(exp_matrix, id2symbol=None, top=50, controls=None, ncols=2,
                    control_name='MT', out_name="TopExpGenes.html",
                    venn_list:list=None, venn_names:list=None):
    controls = [x.strip().split()[0] for x in open(controls)] if controls else []
    id2symbol = dict(x.strip().split()[:2] for x in open(id2symbol)) if id2symbol else dict()
    df = pd.read_csv(exp_matrix, sep='\t', header=0, index_col=0)
    df.index.name = 'id'
    # plot for each sample
    plots = list()
    top_dict = dict()
    detected_gene_dict = dict()
    for sample in df.columns:
        data = df.loc[:, [sample]].sort_values(by=sample, ascending=False)
        detected_gene_dict[sample] = list(data[data[sample]>0].index)
        data['percent'] = data[sample]/data[sample].sum()
        plot_data = data.iloc[:top].copy()
        top_dict[sample] = set(plot_data.index)
        total_percent = plot_data['percent'].sum()
        plot_data['order'] = range(plot_data.shape[0])
        plot_data['symbols'] = [id2symbol[x] if x in id2symbol else x for x in plot_data.index]
        plot_data['marker'] = ['*' if x in controls else "circle" for x in plot_data.index]
        p = figure(
            title="Top {} account for {:.2%} of total in {}".format(top, total_percent, sample),
            tools="wheel_zoom,reset,hover",
            tooltips=[
                ('x', '@percent{0.00%}'),
                ('y', '@symbols'),
            ]
        )
        # min_val = plot_data['percent'].min()
        # max_val = plot_data['percent'].max()
        # mapper = linear_cmap(field_name='percent', palette=bp.Oranges[8], low=min_val, high=max_val)
        p.xaxis.axis_label = '% of Total(~{})'.format(int(data[sample].sum()))
        for marker in ['*', 'circle']:
            source_data = plot_data[plot_data['marker']==marker]
            if source_data.shape[0] == 0:
                continue
            source = ColumnDataSource(source_data)
            p.scatter(
                x='percent',
                y='order',
                # line_color=mapper,
                # color=mapper,
                fill_alpha=0.2,
                size=10,
                marker=marker,
                legend='non-{}'.format(control_name) if marker=="circle" else control_name,
                source=source
            )
        p.yaxis.ticker = list(range(plot_data.shape[0]))
        p.yaxis.major_label_overrides = dict(
            zip(range(plot_data.shape[0]),
                (id2symbol[x] if x in id2symbol else x for x in plot_data.index)
                )
        )
        p.xaxis[0].formatter = NumeralTickFormatter(format="0.00%")
        plots.append(p)
    else:
        fig = gridplot(plots, ncols=ncols)
        output_file(out_name)
        save(fig)
    # plot venn
    if venn_list is None:
        if len(top_dict) <= 6:
            venn.venn(top_dict, cmap="tab10")
            plt.savefig('venn.pdf')
    else:
        if len(venn_list) == 1 and ',' not in venn_list[0]:
            with open(venn_list[0]) as f:
                group_dict = dict(x.strip().split()[:2] for x in f)
                tmp_dict = dict()
                for k, v in group_dict.items():
                    tmp_dict.setdefault(v, set())
                    tmp_dict[v].add(k)
            venn_list = []
            venn_names = []
            for k, v in tmp_dict.items():
                venn_list.append(','.join(v))
                venn_names.append(k)
        if venn_names is None:
            venn_names = []
            for group in venn_list:
                venn_names.append(group.replace(',', '-'))
        for group, name in zip(venn_list, venn_names):
            groups = group.split(',')
            tmp_dict = {x: y for x, y in top_dict.items() if x in groups}
            if 2<= len(tmp_dict) <= 6:
                venn.venn(tmp_dict, cmap="tab10", fmt="{size}\n{percentage:.2f}%", fontsize=9)
                plt.savefig('top{}.{}.venn.pdf'.format(top, name))
            else:
                print('venn for {}?'.format(groups))
                print('venn only support 2-6 sets')
    plt.close()
    # intersection plot
    if venn_list is None:
        plot(from_contents(detected_gene_dict), sum_over=False, sort_categories_by=None, show_counts=True)
        plt.savefig('all.cmbVenn.pdf')
        plt.close()
    else:
        for group, name in zip(venn_list, venn_names):
            groups = group.split(',')
            tmp_dict = {x: y for x,y in detected_gene_dict.items() if x in groups}
            plot(from_contents(tmp_dict), sum_over=False, sort_categories_by=None, show_counts=True)
            plt.savefig('all.{}.cmbVenn.pdf'.format(name))
            plt.close()



if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['plotTopExpGenes'])
