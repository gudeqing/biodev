import re
import sys
import json
import pandas as pd

# infile是womtool生成的wdl输入参数json文件，本脚本将其转化成excel说明。
infile = sys.argv[1]
with open(infile) as f:
    arg_dict = json.load(f)
    result = dict()
    for key, value in arg_dict.items():
        if key.endswith(('.cpu', '.memory', '.disks', '.other_parameters')):
            continue
        tmp = result.setdefault(key, dict())
        tmp['name'] = key.split('.', 1)[1]
        if 'File' in value:
            tmp['type'] = 'infile'
        elif 'Directory' in value:
            tmp['type'] = 'indir'
        elif 'String' in value:
            tmp['type'] = 'str'
        elif 'Int' in value:
            tmp['type'] = 'int'
        elif 'Float' in value:
            tmp['type'] = 'float'
        elif 'Boolean' in value:
            tmp['type'] = 'bool'
        else:
            raise Exception('unknown type for {}'.format(key))

        if 'optional' in value:
            tmp['level'] = 'optional'
        else:
            tmp['level'] = 'required'

        if 'Array' in value:
            tmp['array'] = "yes"
        else:
            tmp['array'] = "no"

        tmp['range'] = ''
        if tmp['type'] == 'bool':
            tmp['range'] = ["true", "false"]

        if 'default' in value:
            default = re.search(r'default =(.*)\)', value).groups()[0].replace('"', '').strip()
            tmp['default'] = default
        else:
            tmp['default'] = ''
        if '-(' in tmp['default']:
            tmp['default'] = tmp['default'].replace('(', '').replace(')', '')

        tmp['format'] = ''

        tmp['order'] = 1
        tmp['desc'] = key

result = {x: result[x] for x in sorted(result.keys())}
data = pd.DataFrame(result).T
data.index.name = 'key'
data.to_excel('args.detail.xlsx')

