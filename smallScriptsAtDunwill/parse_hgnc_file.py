import pandas as pd


class ParseHGNC(object):
    def __init__(self, data):
        """
        :param data: https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit
        """
        target_cols = [
            'HGNC ID',
            'Approved symbol',
            'Approved name',
            'Previous symbols',
            'Synonyms',
            'Ensembl gene ID'
        ]
        self.data = pd.read_csv(data, index_col=0, header=0, sep='\t', usecols=target_cols)

    def withdraw_dict(self):
        result = dict()
        for row in self.data.itertuples():
            if '~withdrawn' in row._1 and 'symbol withdrawn' in row._2:
                key = row._1.split('~withdrawn')[0]
                result.setdefault(key, set())
                result[key].add(row._2.split('see')[1].strip())
        return result

    def ensembl2symbols(self):
        """ensembl_id: [approved_name, approved_symbol, alias, alias2, alias3 ...]"""
        result = dict()
        for row in self.data.itertuples():
            if '~withdrawn' in row._1:
                continue
            if pd.isnull(row._5):
                continue
            result.setdefault(row._5, list())
            if not pd.isnull(row._2):
                result[row._5].append(row._2)
            else:
                result[row._5].append('not_found')
            result[row._5].append(row._1)
            if not pd.isnull(row._3):
                result[row._5] += [x.strip() for x in row._3.split(',')]
            if not pd.isnull(row[4]):
                result[row._5] += [x.strip() for x in row[4].split(',')]

        return result

    def symbol2ensembl(self):
        result = dict()
        for k, v in self.ensembl2symbols().items():
            for sym in v[1:]:
                result.setdefault(sym, set())
                result[sym].add(k)
        return result

    def converting(self, query: list or str, out='symbol2ensembl.txt', known_pair=None, symbol2id=True):
        queries = query if type(query) == list else [x.strip().split()[0] for x in open(query)]
        result = list()
        if symbol2id:
            known_dict = dict(x.strip().split()[:2][::-1] for x in open(known_pair) if x.strip()) if known_pair else dict()
        else:
            known_dict = dict(x.strip().split()[:2] for x in open(known_pair) if x.strip()) if known_pair else dict()
        s2e = self.symbol2ensembl() if symbol2id else self.ensembl2symbols()
        with open(out, 'w') as f:
            not_found = []
            for each in queries:
                if each in known_dict:
                    f.write('{}\t{}\n'.format(known_dict[each], each))
                    result.append(known_dict[each])
                else:
                    not_found.append(each)
            if known_pair:
                print('Success to convert {} genes by querying prior known pair'.format(len(result)))

            # find the remained ones
            withdraw_dict = self.withdraw_dict() if symbol2id else dict()
            failed_ones = []
            for each in not_found:
                if each in s2e:
                    result.append(s2e[each])
                    if len(s2e[each]) > 1:
                        print("{} was found associated with {} genes".format(each, len(s2e[each])))
                    for g in s2e[each]:
                        f.write('{}\t{}\n'.format(g, each))
                elif each in withdraw_dict:
                    print('{} was found in withdraw'.format(each))
                    for new_sym in withdraw_dict[each]:
                        if new_sym in s2e:
                            result.append(s2e[new_sym])
                            for g in s2e[new_sym]:
                                f.write('{}\t{}\n'.format(g, each))
                else:
                    failed_ones.append(each)
                    # print('{} is not found'.format(each))
                    # f.write('{}\t{}\n'.format(each, 'not_found'))
            if not_found:
                print('Success to convert {} genes by query hgnc_custom'.format(len(not_found)-len(failed_ones)))
            if failed_ones:
                print("Failed to query: ")
                print(failed_ones)
        return result


if __name__ == '__main__':
    def converting(query, hgnc_custom="/mnt/resources/HGNC/custom.txt", out='query_result.txt',
                   prior_known_pair=None, symbol2id=False):
        """
        converting ensembl id to symbol or reverse
        :param hgnc_custom: https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit
        :param sym: 待查询的列表文件
        :param out: 输出文件名
        :param prior_known_pair: 已经有的ensembl id 和 symbol对应文件, 包含两列; 如提供，则将优先使用该文件做转换
        :param symbol2id: bool, 如果想把symbol转换为id, 则请用此参数
        """
        object = ParseHGNC(hgnc_custom)
        return object.converting(query=query, symbol2id=symbol2id, out=out, known_pair=prior_known_pair)

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










