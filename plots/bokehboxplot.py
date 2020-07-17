import numpy as np
import pandas as pd
from bokeh.plotting import figure, save, output_file
from bokeh.models import Range1d
from bokeh.models import ColumnDataSource, CDSView, GroupFilter, \
    HoverTool, LabelSet, Legend, CustomJS, TapTool, OpenURL,Div
from bokeh.models.annotations import Title
from bokeh.events import ButtonClick
from bokeh.layouts import gridplot


def boxplot(df:pd.DataFrame, x, y, title=None, y_range=None, dot=False):
    df = df[[y, x]]
    df.columns = ['score', 'group']
    # find the quartiles and IQR for each category
    groups = df.groupby('group', sort=False)
    cats = list(groups.groups.keys())
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    # find the outliers for each category
    def outliers(group):
        cat = group.name
        return group[(group.score > upper.loc[cat]['score']) | (group.score < lower.loc[cat]['score'])]['score']

    out = groups.apply(outliers).dropna()
    # prepare outlier data for plotting, we need coordinates for every outlier.
    outx = []
    outy = []
    if not out.empty:
        for keys in out.index:
            outx.append(keys[0])
            outy.append(out.loc[keys[0]].loc[keys[1]])

    p = figure(
        tools="",
        width=250,
        plot_height=250,
        background_fill_color="#efefef",
        x_range=cats,
        y_range=y_range,
        toolbar_location='left',
        title=title
    )
    # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    upper.score = [min([x,y]) for (x,y) in zip(list(qmax.loc[:,'score']),upper.score)]
    lower.score = [max([x,y]) for (x,y) in zip(list(qmin.loc[:,'score']),lower.score)]

    # stems
    p.segment(cats, upper.score, cats, q3.score, line_color="black")
    p.segment(cats, lower.score, cats, q1.score, line_color="black")

    # boxes
    box_width = 0.7
    p.vbar(cats, box_width, q2.score, q3.score, fill_color="#E08E79", line_color="black")
    p.vbar(cats, box_width, q1.score, q2.score, fill_color="#3B8686", line_color="black")

    # whiskers (almost-0 height rects simpler than segments)
    p.rect(cats, lower.score, 0.2, 0.01, line_color="black")
    p.rect(cats, upper.score, 0.2, 0.01, line_color="black")

    # plot dot
    if dot:
        ymin, ymax = df['score'].min(), df['score'].max()
        # print(ymin, ymax)
        step = (ymax - ymin)/50
        y = ymin
        points = []
        cats_x = {v:k+0.5 for k, v in enumerate(cats, start=0)}
        while y <= ymax:
            in_rg = df[(df['score'] < y+step) & (df['score'] >= y)]
            y += step
            if in_rg.shape[0] == 0:
                continue
            for cat in cats:
                d = in_rg[in_rg['group']==cat]
                num = d.shape[0]
                if num >= 1:
                    center = cats_x[cat]
                    scores = list(d['score'])
                    upper_ = upper.loc[cat]['score']
                    lower_ = lower.loc[cat]['score']
                    colors = ['#F38630' if (x > upper_ or x < lower_) else "gray" for x in scores]
                    x_step = box_width/(num+1)
                    odd = num%2 != 0
                    for i in range(num):
                        if i == 0 and odd:
                            points.append([center, scores[0], colors[i]])
                        else:
                            fold = -1*i//2 if i%2 == 0 else i//2 + 1
                            x_coord = center + x_step*fold
                            points.append([x_coord, scores[i], colors[i]])
        # end of while
        # print(points)
        x_coords = [x[0] for x in points]
        y_coords = [x[1] for x in points]
        colors = [x[2] for x in points]
        # p.circle('group', 'score', size=5, color='gray', fill_alpha=0.6, source=df)
        p.circle(x_coords, y_coords, size=3, color=colors, fill_alpha=0.7)

    # outliers
    if not out.empty and not dot:
        p.circle(outx, outy, size=6, color="#F38630", fill_alpha=0.6)


    # p.circle(outx, outy, size=6, color="#F38630", fill_alpha=0.6)

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = "white"
    p.grid.grid_line_width = 2
    p.xaxis.major_label_text_font_size="16px"
    return p


def boxplots(df, x, y, hue=None, out='boxplots.html', ncols=3, y_range=(0, 0), dot=False):
    if type(df) == str:
        df = pd.read_csv(df, sep=None, engine='python').fillna(0)
    y_range = Range1d(*y_range) if all(y_range) else None
    plots = []
    if hue:
        for each in df[hue].unique():
            tmp = df.loc[df[hue]==each]
            p = boxplot(tmp, x, y, title=each, y_range=y_range, dot=dot)
            plots.append(p)
    else:
        plots = [boxplot(df, x, y, dot=dot)]
    fig = gridplot(
        plots,
        toolbar_location='left',
        sizing_mode='stretch_{}'.format('width'),
        ncols=ncols
    )
    output_file(out, title="boxplot")
    save(fig)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['boxplots'])
