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


def plotTopExpGenes(exp_matrix, id2symbol=None, top=50, controls=None, ncols=2, control_name='MT'):
    controls = [x.strip().split()[0] for x in open(controls)] if controls else []
    id2symbol = dict(x.strip().split()[:2] for x in open(id2symbol)) if id2symbol else dict()
    df = pd.read_csv(exp_matrix, sep='\t', header=0, index_col=0)
    df.index.name = 'id'
    # plot for each sample
    plots = list()
    top_dict = dict()
    for sample in df.columns:
        data = df.loc[:, [sample]].sort_values(by=sample, ascending=False)
        data['percent'] = data[sample]/data[sample].sum()
        plot_data = data.iloc[:top].copy()
        top_dict[sample] = set(plot_data.index)
        total_percent = plot_data['percent'].sum()
        plot_data['order'] = range(plot_data.shape[0])
        plot_data['marker'] = ['*' if x in controls else "circle" for x in plot_data.index]
        p = figure(title="Top {} account for {:.2%} of total in {}".format(top, total_percent, sample))
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
        output_file("TopExpGenes.html")
        save(fig)
    # plot venn
    venn.venn(top_dict, cmap="tab10")
    plt.savefig('venn.pdf')


if __name__ == '__main__':
    class Func2Command(object):
        def __init__(self, callable_dict):
            self.call(callable_dict)

        @staticmethod
        def introduce_command(func):
            import argparse
            import inspect
            import json
            import time
            import sys
            if isinstance(func, type):
                description = func.__init__.__doc__
            else:
                description = func.__doc__
            if '-h' not in sys.argv or '--help' in sys.argv or '-help' in sys.argv:
                description = None
            if description:
                _ = [print(x.strip()) for x in description.split('\n') if x.strip()]
                parser = argparse.ArgumentParser(add_help=False)
            else:
                parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=description)
            func_args = inspect.getfullargspec(func)
            arg_names = func_args.args
            if not arg_names:
                func()
                return
            arg_defaults = func_args.defaults
            if not arg_defaults:
                arg_defaults = []
            arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
            sig = inspect.signature(func)
            for arg, value in zip(arg_names, arg_defaults):
                if arg == 'self':
                    continue
                arg_type = sig.parameters[arg].annotation
                if value == 'None':
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=True, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=True, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=True, metavar=arg)
                elif type(value) == bool:
                    if value:
                        parser.add_argument('--'+arg, action="store_false", help='default: True')
                    else:
                        parser.add_argument('--'+arg, action="store_true", help='default: False')
                elif value is None:
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=False, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=False, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=False, metavar='Default:' + str(value))
                else:
                    if arg_type in [list, tuple, set] or (type(value) in [list, tuple, set]):
                        default_value = ' '.join(str(x) for x in value)
                        if type(value) in [list, tuple]:
                            one_value = value[0]
                        else:
                            one_value = value.pop()
                        parser.add_argument('-' + arg, default=value, nargs='+', type=type(one_value),
                                            metavar='Default:'+default_value, )
                    else:
                        parser.add_argument('-' + arg, default=value, type=type(value),
                                            metavar='Default:' + str(value), )
            if func_args.varargs is not None:
                print("warning: *varargs is not supported, and will be neglected! ")
            if func_args.varkw is not None:
                print("warning: **keywords args is not supported, and will be neglected! ")
            args = parser.parse_args().__dict__
            try:
                with open("Argument_detail.json", 'w') as f:
                    json.dump(args, f, indent=2, sort_keys=True)
            except IOError:
                print('Current Directory is not writable, thus argument log is not written !')
            start = time.time()
            func(**args)
            print("total time: {}s".format(time.time() - start))

        def call(self, callable_dict):
            import sys
            excludes = ['introduce_command', 'Func2Command']
            _ = [callable_dict.pop(x) for x in excludes if x in callable_dict]
            if len(callable_dict) >= 2:
                if len(sys.argv) <= 1:
                    print("The tool has the following sub-commands: ")
                    _ = [print(x) for x in callable_dict]
                    return
                sub_cmd = sys.argv[1]
                sys.argv.remove(sub_cmd)

                if sub_cmd in callable_dict:
                    self.introduce_command(callable_dict[sub_cmd])
                else:
                    print('sub-command: {} is not defined'.format(sub_cmd))
            elif len(callable_dict) == 1:
                self.introduce_command(callable_dict.pop(list(callable_dict.keys())[0]))

    callable_dict = {x: y for x, y in locals().items() if callable(y)}
    _ = [callable_dict.pop(x) for x in {'Func2Command'} if x in callable_dict]
    Func2Command(callable_dict)
