# from fuzzywuzzy import process
from pprint import pprint


def parse_fastq_info(fastq_info_file) -> dict:
    """
    解析fastq输入信息
    :param fastq_info_file:
    format example, at least two column, the third column is optional
    '''
    sample_name | Read1_1_abspath;Read1_2_abspath | Read2_1_abspath;Read2_2_abspath
    sample_name | Read1_abspath | Read2_abspath
    '''
    :return: dict
    """
    fastq_info = dict()
    with open(fastq_info_file) as f:
        for line in f:
            if line.startswith('#') or (not line.strip()):
                pass
            tmp_list = line.strip().split()
            sample, fqs = tmp_list[0], tmp_list[1:]
            fastq_info.setdefault(sample, list())
            read1_list = [x.strip() for x in fqs[0].split(';')]
            fastq_info[sample].append(read1_list)
            if len(fqs) >= 2:
                read2_list = [x.strip() for x in fqs[1].split(';')]
                fastq_info[sample].append(read2_list)
    return fastq_info


def merge_fastq_info(ori: str, again: str, out='new.fastq.info'):
    info = parse_fastq_info(ori)
    info2 = parse_fastq_info(again)
    sample_match_dict = dict()
    for sample in info2.keys():
        match = ''
        for each in info.keys():
            if each[:-4].startswith(sample[:-4]):
                match = each
                break
        if match:
            sample_match_dict[match] = sample
        else:
            raise Exception('{} not found in {}'.format(sample, ori))
    pprint(sample_match_dict)
    final_info = dict()
    for sample in info:
        fastq = info[sample]
        if sample in sample_match_dict:
            sample = sample_match_dict[sample]
            new_fastq = info2[sample]
            for ind in range(len(fastq)):
                fastq[ind] = fastq[ind] + new_fastq[ind]
        final_info[sample] = fastq

    with open(out, 'w') as f:
        for sample, fastq in final_info.items():
            fastq_str = []
            for ind in range(len(fastq)):
                fastq_str.append(';'.join(fastq[ind]))
            f.write('\t'.join([sample] + fastq_str) + '\n')


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
            if isinstance(func, type):
                description = func.__init__.__doc__
            else:
                description = func.__doc__
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
                    parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
                else:
                    if sig.parameters[arg].annotation in [list, tuple]:
                        parser.add_argument('-' + arg, default=value, nargs='+', type=type(value), metavar='Default:' + str(value), )
                    else:
                        parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
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
