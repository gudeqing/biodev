import numpy as np
import pandas as pd
from bokeh.plotting import figure, save, output_file
from bokeh.models import Range1d
from bokeh.models import ColumnDataSource, CDSView, GroupFilter, \
    HoverTool, LabelSet, Legend, CustomJS, TapTool, OpenURL,Div
from bokeh.models.annotations import Title
from bokeh.events import ButtonClick
from bokeh.layouts import gridplot
from bokeh.palettes import Spectral6
from bokeh.transform import factor_cmap


def bar(df, x, y, legend=None, color_col=None, output='bar.html', title='bar plot', palette=None):
    # fruits = ['Apples', 'Pears', 'Nectarines', 'Plums', 'Grapes', 'Strawberries']
    # counts = [5, 3, 4, 2, 4, 6]
    # source = ColumnDataSource(data=dict(fruits=fruits, counts=counts))
    source = ColumnDataSource(data=df)
    fig = figure(x_range=x, toolbar_location=None, title=title)
    fig.vbar(
        x=x,
        top=y,
        width=0.7,
        source=source,
        legend_field=x if legend is None else legend,
        line_color='white',
        fill_color=factor_cmap(color_col or x, palette=palette, factors=sorted(set(df[color_col or x]))) if palette else None
    )
    fig.xgrid.grid_line_color = None
    fig.y_range.start = 0
    # fig.y_range.end = None
    fig.legend.orientation = "vertical"
    fig.legend.location = "top_right"
    output_file(output, title=title)
    save(fig)
