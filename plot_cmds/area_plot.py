import os
from plotly import tools
import plotly.graph_objs as go
from plotly.offline import plot as plt
from functools import partial
plt = partial(plt, auto_open=False)
import pandas as pd


def area_plot(files:list, prefix='area'):
    traces = []
    labels = []
    x_list = [0]
    tickvals = []
    for ind, each_file in enumerate(files):
        df = pd.read_csv(each_file, index_col=0, sep='\t')
        x_list = [x + max(x_list) for x in range(df.shape[0])]
        labels += list(df.index)
        tickvals += x_list
        y_list = df.iloc[:, 0]
        trace = go.Scatter(
            x=x_list,
            y=y_list,
            fill='tozeroy',
            fillcolor='lightgrey',
            line=dict(
                color='lightgrey'
            ),
            text="lightgrey",
            hoverinfo='text'
        )
        traces.append(trace)
    print(labels)
    layout = go.Layout(
        xaxis=dict(
            tickvals=tickvals,
            ticktext=labels,
            showticklabels=True,
            dtick=1,
        )
    )

    fig = go.Figure(data=traces, layout=layout)
    plt(fig, filename=prefix+'.html')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['area_plot'])
