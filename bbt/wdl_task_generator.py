import sys
import configparser
# import pandas as pd


def get_task_meta(data, section='tool'):
    if section in data.sections():
        return data[section]
    else:
        return {'author': '?'}


def get_runtime(data, section='runtime'):
    if section not in data.sections():
        return {'docker': 'your_image_id'}
    else:
        return data[section]


def get_parameter_meta(data):
    arg_meta = dict()
    # keys = ['type', 'desc', 'level', 'default', 'range']
    for section in data.sections():
        detail = data[section]
        if {'name', 'type', 'desc'} - set(detail.keys()):
            # 判断当前section不是参数描述
            continue
        meta = arg_meta.setdefault(detail['name'], dict())
        meta['desc'] = detail['desc']

        # get level
        meta['level'] = 'required' if detail['require'] == 'yes' else 'optional'

        # get type
        if detail['is_infile'] == 'yes':
            meta['type'] = 'infile'
        elif detail['input_dir'] == 'yes':
            meta['type'] = 'indir'
        else:
            meta['type'] = detail['type']

        # get ranges
        if 'value_candidates' in detail and detail['value_candidates'] != "none":
            # for old style
            meta['range'] = detail['value_candidates']
        else:
            meta['range'] = ''

        # get default
        if detail['default'] != "none":
            meta['default'] = detail['default']
        else:
            meta['default'] = ''
    return arg_meta


def get_inputs_and_cmd(data):
    inputs = []
    cmd = []
    for section in data.sections():
        arg_info = ''
        detail = data[section]
        if {'name', 'type', 'desc'} - set(detail.keys()):
            # 判断当前section不是参数描述
            continue
        detail.setdefault('multi_values', 'no')
        detail.setdefault('multi_times', 'no')
        detail.setdefault('input_dir', 'no')
        detail.setdefault('out_dir', 'no')
        detail.setdefault('value_candidates', 'none')

        # define type of arg
        if detail['type'] == 'bool':
            arg_info = 'Boolean'
        elif detail['type'] == 'int':
            arg_info = 'Int'
        elif detail['type'] == 'float':
            arg_info = 'Float'
        elif detail['type'] == 'str' and detail['is_infile'] == 'no':
            arg_info = 'String'
        elif detail['type'] == 'str' and (detail['is_infile'] == 'yes' or detail['is_outfile'] == 'yes'):
            arg_info = 'File'

        if detail['multi_values'] == 'yes':
            arg_info = f'Array[{arg_info}]'

        if detail['require'] == 'no' and detail['type'] != 'bool':
            arg_info += '?'

        # add arg name
        arg_info += ' ' + detail['name']

        # get default value
        if detail['type'] == "bool":
            if detail['default'] == 'yes':
                arg_info += ' = true'
            else:
                arg_info += ' = false'
        elif detail['type'] == "str":
            if detail['default'] != 'none':
                arg_info += ' = "{}"'.format(detail['default'])
        else:
            if detail['default'] != 'none':
                arg_info += ' = {}'.format(detail['default'])

        inputs.append(arg_info)

        # format cmd
        if detail['multi_values'] == 'no':
            if detail['prefix'] == 'none':
                cmd += ['~{' + detail['name'] + '}']
            else:
                if detail['type'] != 'bool':
                    cmd += ['~{' + f'"{detail["prefix"]} "' + ' + ' + detail['name'] + '}']
                else:
                    cmd += ['~{' + f'if {detail["name"]} then "{detail["prefix"]} " else ""' + '}']
        else:
            if 'multi_value_sep' not in detail or (not detail['multi_value_sep']):
                detail['multi_value_sep'] = ' '
            delimiter = detail['multi_value_sep']

            if detail['prefix'] == 'none':
                cmd += ['~{sep=' + f'"{delimiter}" ' + detail['name'] + '}']
            else:
                if detail['multi_times'] == 'yes':
                    cmd += ['~{' + 'prefix(' + f'"{detail["prefix"]} ", ' + detail['name'] + ')}']
                else:
                    if detail['type'] != 'bool':
                        prefix = '~{' + f'if defined({detail["name"]}) then "{detail["prefix"]} " else ""' + '}'
                        cmd += [prefix + '~{sep=' + f'"{delimiter}" ' + detail['name'] + '}']
                    else:
                        raise Exception('不支持Bool Array！')
    return inputs, cmd


def get_outputs(data):
    outputs = []
    if 'outputs' not in data.sections():
        for section in data.sections():
            detail = data[section]
            if {'name', 'type', 'desc'} - set(detail.keys()):
                # 判断当前section不是参数描述
                continue
            detail.setdefault('input_dir', 'no')
            detail.setdefault('out_dir', 'no')
            if detail['is_outfile'] == 'yes':
                # 避免输出名和输入名重复，需要改名
                outputs += [f'File o_{detail["name"]} = ' + '~{' + detail["name"] + '}']
            elif detail['out_dir'] == 'yes':
                outputs += [f'Array[File] outputs = glob(".*")']
        if not outputs:
            outputs += [f'Array[File] outputs = glob(".*")']
    else:
        for k, v in data['outputs'].items():
            name = k
            typ, value = v.split(' ', 1)
            outputs += [typ + ' ' + name + ' = ' + value]
            # outputs += [f'{k} = {v}']

    return outputs


def format_wdl_task(data):
    """
    version = 1.0

    task cmd_name {
        input {
            xxx
        }
        command <<<
            xxx
        >>>
        output {
            xxx
        }
        runtime {
            xxx
        }
        meta {
            xxx
        }
        parameter_meta {
            xxx
        }
    }
    :return: *.wdl file
    """
    task_name = data['tool']['name']
    lines = f'task {task_name}' + '{\n'
    # input
    lines += ' '*4 + 'input {\n'
    inputs, cmd = get_inputs_and_cmd(data)
    for line in inputs:
        lines += ' '*8 + line + '\n'
    lines += ' '*8 + '# for runtime\n'
    for k, v in get_runtime(data).items():
        if k in ['cpu', 'time_minutes']:
            typ = 'Int'
        else:
            typ = 'String'
            v = f'"{v}"'
        lines += ' ' * 8 + f'{typ} {k} = {v}' + '\n'
    lines += ' '*4 + '}\n\n'

    # command
    lines += ' ' * 4 + 'command <<<\n'
    lines += ' '*8 + 'set -e \n'
    base_cmd = data['tool'].setdefault('baseCmd', 'xxx')
    lines += ' '*8 + base_cmd + ' \\\n'
    for line in cmd:
        lines += ' '*8 + line + ' \\\n'
    lines = lines[:-2] + '\n'
    lines += ' ' * 4 + '>>>\n\n'

    # output
    lines += ' ' * 4 + 'output {\n'
    for line in get_outputs(data):
        lines += ' ' * 8 + line + '\n'
    lines += ' ' * 4 + '}\n\n'

    # runtime
    lines += ' ' * 4 + 'runtime {\n'
    for k, v in get_runtime(data).items():
        # lines += ' ' * 8 + f'{k}: "{v}"' + '\n'
        lines += ' ' * 8 + f'{k}: {k}' + '\n'
    lines += ' ' * 4 + '}\n\n'

    # task meta
    lines += ' ' * 4 + 'meta {\n'
    for k, v in get_task_meta(data, section='tool').items():
        lines += ' ' * 8 + f'{k}: "{v}"' + '\n'
    lines += ' ' * 4 + '}\n\n'

    # parameter_meta
    lines += ' ' * 4 + 'parameter_meta {\n'
    import json
    for k, v in get_parameter_meta(data).items():
        v_info = '{' + ', '.join(f'{k}: "{v}"' for k, v in v.items()) + '}'
        lines += ' ' * 8 + f'{k}: {v_info}' + '\n'
    lines += ' ' * 4 + '}\n\n'

    # write wdl file
    lines += '}\n'
    with open(f'{task_name}.wdl', 'w', encoding='utf-8') as f:
        f.write('version 1.0\n\n')
        f.write(lines)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python wdl_task_generator.py <cmd.ini>')
    in_file = sys.argv[1]
    # parser = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    parser = configparser.ConfigParser()
    parser.read(in_file, encoding='utf-8')
    if '参数属性说明' in parser.sections():
        parser.remove_section('参数属性说明')
    format_wdl_task(parser)

