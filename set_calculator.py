def run(files:list, exp=None, out=None, has_header=False,
        intersect_only=True, intersect_x_outof_y=1, union_only=False):
    set_number = len(files)
    if len(files) >= 2:
        for ind, each in enumerate(files, start=1):
            exec('s{}=set(open("{}").readlines())'.format(ind, each))
    else:
        import pandas as pd
        table = pd.read_table(files[0], header=0 if has_header else None)
        set_number = table.shape[1]
        for i in range(table.shape[1]):
            exec('s{}=set(table.iloc[:, {}])'.format(i+1, i))

    result = list()
    if exp:
        print("do as you say in exp")
        result = eval(exp)
    elif intersect_x_outof_y > 1:
        print('do intersect_x_outof_y')
        union = eval('|'.join(['s'+str(x) for x in range(1, set_number+1)]))
        result = set()
        for each in union:
            varspace = dict(locals())
            in_times = sum(eval("each in s{}".format(x), varspace) for x in range(1, set_number+1))
            if in_times >= intersect_x_outof_y:
                result.add(each)
    elif union_only:
        print('do union only')
        result = eval('|'.join(['s'+str(x) for x in range(1, set_number+1)]))
    elif intersect_only:
        print('do intersect only')
        result = eval('&'.join(['s'+str(x) for x in range(1, set_number+1)]))
    if not result:
        print('result is empty!')
    else:
        print('result size: {}'.format(len(result)))
    with open(out or 'result.list', 'w') as f:
        _ = [f.write(x) for x in result]


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
                arg_type = sig.parameters[arg].annotation
                if arg == 'self':
                    continue
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
            # print("total time: {}s".format(time.time() - start))

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
