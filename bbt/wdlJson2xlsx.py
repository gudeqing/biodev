import re
import sys
import json
import pandas as pd

# infile是womtool生成的wdl输入参数json文件，本脚本将其转化成excel说明。
print('usage: python wdlJson2xlsx.py <input.json> [wf.wdl] ')

infile = sys.argv[1]

with open(infile) as f:
    arg_dict = json.load(f)
    result = dict()
    for key, value in arg_dict.items():
        if key.endswith(('.cpu', '.memory', '.disks', '.other_parameters', '.time_minutes', '.docker')):
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

#  获取参数注释
if len(sys.argv) >= 3:
    wdl_file = sys.argv[2]
    with open(wdl_file) as f:
        task_names = re.findall(r'task\s+(.*)\{', f.read())

    with open(wdl_file) as f:
        match = re.findall(r'call\s+(.*)\s+as\s+(.*)\s+\{', f.read())
        print(match)
        if match:
            task_mapping = dict(match)
        else:
            task_mapping = dict()
    task_names = [task_mapping[x] if x in task_mapping else x for x in task_names]

    with open(wdl_file) as f:
        # 最后根据是"\s}\n"定位parameter的结束
        arg_desc = re.findall(r'parameter_meta\s+\{(.*?)\s+\}\n', f.read(), flags=re.S)
        arg_desc = [x.strip() for x in arg_desc]

    desc_dict = dict(zip(task_names, arg_desc))

    for k, v in result.items():
        name = v['name']
        if len(name.split('.', 1)) == 1:
            continue
        task = name.split('.', 1)[0]
        para = name.split('.', 1)[1]
        task_para_desc = [x.strip() for x in desc_dict[task].split('\n')]
        target_desc_detail = [x for x in task_para_desc if x.startswith(para+':')]
        if target_desc_detail:
            desc = re.search(r'desc:\s+"(.*?)"', target_desc_detail[0]).groups()[0]
            v['desc'] = desc
        else:
            print(f'no desc found for {name} !')

data = pd.DataFrame(result).T
data.index.name = 'key'
data.to_excel('args.detail.xlsx')

