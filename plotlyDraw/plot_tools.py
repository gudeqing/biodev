import os
from functools import partial
import colorlover
from plotly import tools
import plotly.graph_objs as go
from plotly.offline import plot as plt
plt = partial(plt, auto_open=False)
import pandas as pd


def gene_body_coverage(files):
    layout = go.Layout(title="geneBodyCoverage")
    all_fig = go.Figure(layout=layout)
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=None, index_col=0)
        fig = go.Figure(layout=go.Layout(title='{}'.format(sample)))
        normalize_y = data.iloc[1, :]/data.iloc[1, :].max()
        fig.add_scatter(x=data.iloc[0, :], y=normalize_y, name=sample)
        all_fig.add_scatter(x=data.iloc[0, :], y=normalize_y, name=sample)
        plt(fig, filename='{}.geneBodyCoverage.html'.format(sample))
    plt(all_fig, filename='samples.geneBodyCoverage.html')


def fragment_length(files, outdir, min_len=50, max_len=600):
    layout = go.Layout(
        title="Fragment length distribution",
        xaxis=dict(title='Fragment length'),
        yaxis=dict(title='Probability'),
    )
    all_fig = go.Figure(layout=layout)
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=0, index_col=0)
        data = data.loc[(data['frag_median'] >= min_len) & (data['frag_median'] <= max_len)]
        fig = go.Figure(layout=go.Layout(
            title='{}'.format(sample),
            xaxis=dict(title='Fragment length'),
            yaxis=dict(title='probability'),
        ))
        fig.add_histogram(x=data["frag_median"], histnorm='probability', name=sample)
        all_fig.add_histogram(x=data["frag_median"], histnorm='probability', name=sample)
        plt(fig, filename="{}.fragmentLengthDistribution.html".format(sample))
    plt(all_fig, filename="samples.fragmentLengthDistribution.html")


def inner_distance(files, outdir, min_dist=-250, max_dist=250):
    """
    抽样1000000得到的统计结果，第三列的值可能是：
    readPairOverlap
    sameTranscript=No,dist=genomic
    sameTranscript=Yes,nonExonic=Yes,dist=genomic
    sameTranscript=Yes,sameExon=No,dist=mRNA
    sameTranscript=Yes,sameExon=Yes,dist=mRNA
    """
    groups = [
        'sameTranscript=No,dist=genomic',
        'sameTranscript=Yes,nonExonic=Yes,dist=genomic',
        'readPairOverlap',
        'sameTranscript=Yes,sameExon=No,dist=mRNA',
        'sameTranscript=Yes,sameExon=Yes,dist=mRNA',
    ]
    layout = go.Layout(
        title="Inner distance distribution",
        xaxis=dict(title='Inner Distance'),
        yaxis=dict(title='Probability'),
    )
    all_fig = go.Figure(layout=layout)
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=None, index_col=0)
        values = [data[2][data[2]==x].shape[0] for x in groups]
        norm_values = [x/sum(values) for x in values]
        percents = dict(zip(groups, norm_values))
        data = data[(data[1] >= min_dist) & (data[1] <= max_dist)]
        fig = go.Figure(layout=go.Layout(
                title='{}'.format(sample),
                xaxis=dict(title='Inner Distance'),
                yaxis=dict(title='Frequency'),
        ))
        for each in groups:
            target = data[1][data[2]==each]
            name = "{}({:.2%})".format(each, percents[each])
            fig.add_histogram(x=target, name=name)
        all_fig.add_histogram(x=data[1][data[2]!="sameTranscript=No,dist=genomic"], histnorm='probability', name=sample)
        plt(fig, filename="{}.InnerDistanceDistribution.html".format(sample))
    html_file = os.path.join(outdir, 'samples.InnerDistanceDistribution.html')
    plt(all_fig, filename=html_file)


def read_distribution(files, outdir):
    all_data = list()
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, index_col=0, skiprows=4, skipfooter=1, sep='\s+')
        data = data.loc[:, 'Tag_count']
        all_data.append(data)
        fig = go.Figure(layout=go.Layout(
            title='{}'.format(sample),
        ))
        fig.add_pie(labels=data.index, values=data)
        plt(fig, filename="{}.ReadDistribution.html".format(sample))
    df = pd.concat(all_data, axis=1).T
    data = [go.Bar(x=df.index, y=df[x], name=x) for x in df.columns]
    layout = go.Layout(
        title="Read distribution",
        xaxis=dict(title='Sample'),
        barmode='stack'
    )
    fig = go.Figure(data=data, layout=layout)
    html_file = os.path.join(outdir, 'samples.ReadDistribution.html')
    plt(fig, filename=html_file)


def read_duplication(pos_dup_files, outdir=os.getcwd(), max_dup=500):
    traces = list()
    for each in pos_dup_files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=0, index_col=None)
        data = data[data.iloc[:, 0] <= max_dup]
        trace = go.Scatter(x=data.iloc[:, 0], y=data.iloc[:, 1], name=sample, mode='markers')
        traces.append(trace)
        layout = go.Layout(
            title=sample,
            xaxis=dict(title='Occurrence'),
            yaxis=dict(title='UniqReadNumber', type='log'),
        )
        fig = go.Figure([trace], layout=layout)
        out_name = os.path.join(outdir, "{}.ReadDuplication.html".format(sample))
        plt(fig, filename=out_name)

    layout = go.Layout(
        title="Read duplication",
        xaxis=dict(title='Occurrence'),
        yaxis=dict(title='UniqReadNumber', type='log'),
    )
    fig = go.Figure(traces, layout=layout)
    out_name = os.path.join(outdir, "samples.ReadDuplication.html")
    plt(fig, filename=out_name)


def exp_saturation(exp_files, outdir=os.getcwd(), outlier_limit=5):
    all_fig = tools.make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Q1',
            'Q2',
            'Q3',
            'Q4',
        ),
        shared_xaxes=False
    )
    color_pool = colorlover.scales['{}'.format(len(exp_files)+1)]['div']['RdYlBu']

    for exp_file, sample_color in zip(exp_files, color_pool):
        sample = os.path.basename(exp_file).split('.', 1)[0]
        data = pd.read_table(exp_file, header=0, index_col=0)
        # plot deviation
        describe = data['100'].describe()
        regions = [
            (describe['min'], describe['25%']),
            (describe['25%'], describe['50%']),
            (describe['50%'], describe['75%']),
            (describe['75%'], describe['max']),
        ]
        errors = data.apply(lambda column: column / data.loc[:, '100'], axis=0)
        errors = ((errors - 1)*100).abs()
        plot_data = list()
        for lower, upper in regions:
            tmp = errors[(data['100'] >= lower) & (data['100'] <= upper)]
            plot_data.append(tmp)

        fig = tools.make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Q1: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[0]),
                'Q2: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[1]),
                'Q3: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[2]),
                'Q4: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[3]),
            ),
            shared_xaxes=False
        )
        for ind, each in enumerate(plot_data):
            x = 1 if (ind+1) / 2 <= 1 else 2
            y = 1 if (ind+1) % 2 == 1 else 2
            for col in each.columns:
                upper_limit = each[col].describe()['50%']*outlier_limit
                y_data = each[col][each[col] <= upper_limit]
                # box = go.Box(y=y_data, name=col, showlegend=False, boxpoints=False)
                box = go.Box(y=y_data, name=col, showlegend=False)
                fig.append_trace(box, x, y)
            line = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines', name='median')
            fig.append_trace(line, x, y)
            if ind != 0:
                line2 = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines',
                                   legendgroup=sample, name=sample, line=dict(color=sample_color), showlegend=False)
            else:
                line2 = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines',
                                   legendgroup=sample, name=sample, line=dict(color=sample_color))
            all_fig.append_trace(line2, x, y)
        out_name = os.path.join(outdir, "{}.TPM.Saturation.html".format(sample))
        plt(fig, filename=out_name)

    # plot all sample together
    out_name = os.path.join(outdir, "samples.TPM.Saturation.html")
    plt(all_fig, filename=out_name)


