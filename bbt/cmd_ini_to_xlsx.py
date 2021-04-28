import sys
import configparser
import pandas as pd

in_file = sys.argv[1]

parser = configparser.ConfigParser()
parser.read(in_file, encoding='utf-8')

target_sections = [x for x in parser.sections() if x not in ['tool', 'cmd', '参数属性说明', 'other_args', ]]
result = dict()
for arg_name in target_sections:
    info = dict()
    detail = parser[arg_name]
    detail.setdefault('multi_values', 'no')
    detail.setdefault('multi_times', 'no')
    detail.setdefault('input_dir', 'no')
    detail.setdefault('out_dir', 'no')
    detail.setdefault('value_candidates', 'none')
    # info['name'] = info['name']
    info['prefix'] = detail['prefix'] if detail['prefix'] != 'none' else ''
    if detail['format'] == 'none':
        info['format'] = ''
    else:
        info['format'] = ','.join(x.strip()[1:].upper() for x in detail['format'].split(','))
    info['gen_rule'] = 'returnByType'
    if detail['prefix'] == 'none' and detail['type'] == 'bool':
        info['gen_rule'] = 'returnBooleanShowFlag'
    elif detail['prefix'] == 'none' and detail['type'] == 'str':
        info['gen_rule'] = 'returnValue'

    if detail['default'] == 'none':
        info['default'] = ''
    else:
        info['default'] = detail['default']
    if detail['require'] == 'yes':
        info['required'] = 1
    else:
        info['required'] = 0
    info['enum_value'] = ''
    info['category'] = 'PARAM'
    info['type'] = 'STRING'
    if detail['type'] == 'int' or detail['type'] == 'float':
        info['type'] = 'NUMBER'
    elif detail['type'] == 'bool':
        info['type'] = 'BOOLEAN'
        info['default'] = 'TRUE' if detail['default'] == 'yes' else 'FALSE'
    elif detail['input_dir'] == 'yes':
        info['type'] = 'DIR'
        info['category'] = 'INPUTS'
    elif detail['value_candidates'] != 'none':
        info['type'] = 'ENUM'
        info['enum_value'] = ','.join(x.strip() for x in detail['value_candidates'].split(','))
    elif detail['type'] == 'str':
        if detail['is_infile'] == 'yes':
            info['category'] = 'INPUTS'
            if detail['multi_values'] == 'yes':
                info['type'] = 'FILE_ARRAY'
            else:
                info['type'] = 'FILE'
        else:
            if detail['multi_values'] == 'yes':
                info['type'] = 'STRING_ARRAY'
    else:
        pass
    info['description'] = detail['desc']
    result[detail['name']] = info


pd.DataFrame(result).T.to_excel(f'{parser["tool"]["name"]}.arg_detail.xlsx')
