import re
import sys
import json
import pandas as pd

# infile是womtool生成的wdl输入参数json文件，本脚本将其转化成excel说明。
if len(sys.argv) < 3:
    print('usage: python wdlJson2xlsx.py <input.json> [wf.wdl] ')
    #  获取参数注释, 要求wdl文件中每一个task定义时都定义parameter_meta，而且parameter_meta中参数的描述中必需包含desc字段，如下：
    #  parameter_meta {
    #         sample: {desc: "sample name", type='xx', other='xxx'}
    #         genome: {desc: "genome fasta file"}

infile = sys.argv[1]
attr_annot = {
    "name": "参数名,将在页面上显示，如star.index_dir；不填则默认通过分割key获得",
    "type": "表示参数类型，支持填写以下类型：str,int, float, bool, infile(表示为输入文件）, indir（表示为输入目录）；不填写则默认为str类型",
    "level": "表示参数是否必需，用‘optional’ 或‘required’表示；不填则默认为optional",
    "array": "表示参数是否接收多个值，用‘yes' 或 'no'表示；不填则默认为no",
    "range": "表示参数的可选内容，如['A','B','C']；不填则表示参数无可选内容",
    "default": "表示参数的默认值, 不填则表示该参数无默认值",
    "format": "表示对输入文件的格式要求说明，例如'vcf,vcf.gz'；不填则默认对格式无要求",
    "order": "表示参数的先后顺序, 数字越小越靠前；不填则默认为1",
    "desc": "表示参数说明，不填则默认显示key填写的内容",
}

key_annot = "表示参数ID,必填项，需跟流程中定义的保持一致,如pipeline.star.index_dir"
new_arg_dict = dict()
with open(infile) as f:
    arg_dict = json.load(f)
    result = dict()
    for key, value in arg_dict.items():
        if key.endswith(('.cpu', '.memory', '.disks', '.other_parameters', '.time_minutes', '.docker')):
            continue
        tmp = result.setdefault(key, dict())
        if len(key.split('.')) > 2:
            tmp['name'] = key.split('.', 1)[1]
        else:
            tmp['name'] = key
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
            # default = re.search(r'default =(.*)\)', value).groups()[0].strip()
            tmp['default'] = default
        else:
            if tmp['array'] == 'yes' and tmp['level'] == 'required':
                tmp['default'] = '[]'
            else:
                tmp['default'] = ''
        if '-(' in tmp['default']:
            tmp['default'] = tmp['default'].replace('(', '').replace(')', '')
        new_arg_dict[key] = tmp['default']
        tmp['format'] = ''
        tmp['order'] = 1
        tmp['desc'] = key
    pipeline_name = key.split('.')[0]

# output new json
with open(infile[:-4] + 'detailed.json', 'w') as f:
    json.dump(new_arg_dict, f, indent=2, sort_keys=True)

result[key_annot] = attr_annot
result = {x: result[x] for x in [key_annot] + sorted(set(result.keys()) - {key_annot})}

#  获取参数注释, 要求wdl文件中每一个task定义时都定义parameter_meta，而且parameter_meta中参数的描述中必需包含desc字段，如下：
#  parameter_meta {
#         sample: {desc: "sample name", type='xx', other='xxx'}
#         genome: {desc: "genome fasta file"}

if len(sys.argv) >= 3:
    wdl_file = sys.argv[2]
    # 获取task 列表
    with open(wdl_file) as f:
        task_names = re.findall(r'task\s+(.*)\{', f.read())

    # 获取task别名和task本名的映射信息
    with open(wdl_file) as f:
        match = re.findall(r'call\s+(.*)\s+as\s+(.*)\s+\{', f.read())
        if match:
            task_mapping = dict(match)
            task_mapping = {v: k for k, v in task_mapping.items()}
        else:
            task_mapping = dict()

    # 提取每个task对应的parameter_meta
    with open(wdl_file) as f:
        # 最后根据是"\s}\n"定位parameter的结束
        arg_desc = re.findall(r'parameter_meta\s+\{(.*?)\s+\}\n', f.read(), flags=re.S)
        arg_desc = [x.strip() for x in arg_desc]

    # 将task和parameter_meta按照顺序对应起来, 因为是按序对应，因此要求每个task都有parameter_meta部分，否则容易错配
    desc_dict = dict(zip(task_names, arg_desc))

    # 获取desc信息
    for k, v in result.items():
        if k == key_annot:
            continue
        name = v['name']
        if name.split('.')[0] == pipeline_name:
            continue
        if len(name.split('.', 1)) == 1:
            continue
        task = name.split('.', 1)[0]
        if task in task_mapping:
            task = task_mapping[task]
        para = name.split('.', 1)[1]
        task_para_desc = [x.strip() for x in desc_dict[task].split('\n')]
        target_desc_detail = [x for x in task_para_desc if x.startswith(para+':')]
        if target_desc_detail:
            desc = re.search(r'desc:\s+"(.*?)"', target_desc_detail[0]).groups()[0]
            v['desc'] = desc
        else:
            print(f'no desc found for {name} !')

# 生成最后的xlsx格式的参数说明文件，用于XDP
data = pd.DataFrame(result).T
data['order'] = [data['order'][0]] + list(range(data.shape[0]-1))
data.index.name = 'key'
data.to_excel('args.detail.xlsx')

