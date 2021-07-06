import os
import json
from uuid import uuid4
from dataclasses import dataclass, field
from typing import Any, List, Dict

"""
设计思路
1. 定义argument 类, argument是tool的某个参数
2. 定义command：[docker] + [tool_dir] + <tool> + <arguments> + [outputs]
# 因为不同参数的选择可能导致不同的outputs，所以定义task时可以重新定义outputs
3. 定义task: command_with_arg_value_defined + <outputs> + <depend>

"""


@dataclass()
class Argument:
    name: str = 'name'
    value: Any = None
    # prefix 可以是如 ’-i '或 'i=', 对于前者, 空格一定要带上
    prefix: str = ''
    # type is one of ['str', 'int', 'float', 'bool', 'infile', 'indir']
    type: str = 'str'
    level: str = 'required'
    # for bool type, default is one of ['false', true']
    default: Any = None
    range: Any = None
    array: bool = False
    delimiter: str = ' '
    # 指示一个参数是否可以多次使用
    multi_times: bool = False
    format: str = None
    order: int = 0
    desc: str = 'This is description of the argument.'

    def __post_init__(self):
        if self.type == 'bool':
            # 对于布尔参数，其一定为必要参数类型,可选范围为True或False
            self.level = 'required'
            self.range = {True, False}
            if not self.default:
                # 对于没有默认值的bool参数，强行赋值为false，该参数默认不参与命令行的形成
                self.default = False
        # 对于有集合候选的参数，先对默认值检查一番
        if type(self.range) == set:
            if self.default not in self.range:
                raise Exception(f'default value is not in {self.range}')


@dataclass()
class RunTime:
    image: str = None
    memory: int = None
    cpu: int = 1
    tool_dir: str = ''
    tool: str = ''


@dataclass()
class Output:
    path: str
    out_id: str = uuid4()
    # type should be one of ['File', 'Directory']
    type: str = 'File'
    # 设计locate 参数用于整理结果目录
    locate: str = "report"

    def __post_init__(self):
        if self.type not in ['File', 'Directory']:
            raise Exception("output type should be one of ['File', 'Directory']")


@dataclass()
class Meta:
    name: str = None
    desc: str = 'This is description of the tool/workflow.'
    author: str = 'unknown'
    source: str = 'source URL for the tool'


@dataclass()
class Command:
    meta: Meta = field(default_factory=Meta)
    runtime: RunTime = field(default_factory=RunTime)
    args: Dict[str, Argument] = field(default_factory=dict)
    # 下面支持用wdl的定义’~‘符号, 当前脚本要求所有命令行的输出皆为文件
    outputs: Dict[str, Output] = field(default_factory=dict)

    def format_cmd(self):
        # 该函数应该在参数被具体赋值后才能调用
        if not self.args:
            raise Exception(f'Command {self.name} has no args !')
        cmd = ''
        if self.runtime.tool_dir and not self.runtime.tool_dir.endswith('/'):
            self.runtime.tool_dir += '/'
        cmd += self.runtime.tool_dir + self.runtime.tool

        for arg_name, arg in self.args.items():
            # 对于极端情况:如果不小心定义了"有默认值的非必须参数"，仅当明确赋值时才参与命令行的形成
            if arg.level == 'required':
                arg_value = arg.value or arg.default
            else:
                arg_value = arg.value

            if not arg_value:
                if arg.level == 'required':
                    raise Exception(f'No value found for {arg_name}!')
                else:
                    # 对于非必须参数，且没有赋值的参数，直接跳过，不参与命令行的形成
                    continue
            # 对于可以接收多个值的参数
            if arg.array:
                arg_value = arg.delimiter.join([str(x) for x in arg_value])
            # 处理bool值参数
            if arg.type == "bool":
                if arg_value:
                    cmd += ' ' + arg.prefix
                else:
                    # 如果说bool类型且value为false，则该参数不参与命令行形成
                    continue
            else:
                cmd += ' ' + arg.prefix + str(arg_value)
        return cmd

    def format_wdl_task(self, outfile=None, wdl_version='development'):
        if not outfile:
            outfile = f'{self.name}.wdl'
        ToWdlTask(self, outfile, wdl_version)


@dataclass()
class Task:
    cmd: Command
    task_id: str = uuid4()
    depend: List[str] = field(default_factory=list)
    outputs: dict = field(default_factory=dict)


@dataclass()
class Workflow:
    meta: Meta = field(default_factory=Meta)
    tasks: Dict[str, Task] = field(default_factory=dict)
    outputs: Dict[str, Output] = field(default_factory=dict)


class ToWdlTask(object):
    def __init__(self, command: Command, outfile, wdl_version='development'):
        self.cmd = command
        self.wdl_version = wdl_version
        self.format_wdl_task(outfile, wdl_version)

    def get_task_meta(self):
        return self.cmd.meta.__dict__.copy()

    def get_runtime(self):
        runtime = self.cmd.runtime.__dict__.copy()
        runtime.pop('tool')
        runtime.pop('tool_dir')
        if self.wdl_version in ["development", "1.0"]:
            runtime['docker'] = runtime.pop('image')
        return runtime

    def get_parameter_meta(self):
        arg_meta = dict()
        for arg_name, detail in self.cmd.args.items():
            detail = detail.__dict__.copy()
            detail.pop('value')
            detail.pop('order')
            detail.pop('format')
            detail.pop('multi_times')
            arg_meta[arg_name] = detail
        return arg_meta

    def get_inputs_and_cmd(self):
        inputs = []
        cmd = []
        for arg_name, detail in self.cmd.args.items():
            arg_info = ''
            detail = detail.__dict__
            # define type of arg
            if detail['type'] == 'bool':
                arg_info = 'Boolean'
            elif detail['type'] == 'int':
                arg_info = 'Int'
            elif detail['type'] == 'float':
                arg_info = 'Float'
            elif detail['type'] == 'str':
                arg_info = 'String'
            elif detail['type'] == 'infile':
                arg_info = 'File'
            elif detail['type'] == 'indir':
                arg_info = 'Directory'

            if detail['array']:
                arg_info = f'Array[{arg_info}]'

            if detail['level'] == 'optional' and detail['type'] != 'bool':
                arg_info += '?'

            # add arg name
            arg_info += ' ' + arg_name

            # get default value
            if detail['type'] == "bool":
                if detail['default']:
                    arg_info += ' = true'
                else:
                    arg_info += ' = false'
            elif detail['type'] == "str":
                if detail['default'] is not None:
                    arg_info += ' = "{}"'.format(detail['default'])
            else:
                if detail['default'] is not None:
                    arg_info += ' = {}'.format(detail['default'])

            inputs.append(arg_info)

            # format cmd
            if not detail['array']:
                if detail['prefix'] == '':
                    cmd += ['~{' + arg_name + '}']
                else:
                    if detail['type'] != 'bool':
                        cmd += ['~{' + f'"{detail["prefix"]}"' + ' + ' + arg_name + '}']
                    else:
                        cmd += ['~{' + f'if {arg_name} then "{detail["prefix"]} " else ""' + '}']
            else:
                delimiter = detail['delimiter']
                if detail['prefix'] == '':
                    cmd += ['~{sep=' + f'"{delimiter}" ' + arg_name + '}']
                else:
                    if detail['multi_times']:
                        cmd += ['~{' + 'prefix(' + f'"{detail["prefix"]} ", ' + arg_name + ')}']
                    else:
                        if detail['type'] != 'bool':
                            prefix = '~{' + f'if defined({arg_name}) then "{detail["prefix"]} " else ""' + '}'
                            cmd += [prefix + '~{sep=' + f'"{delimiter}" ' + arg_name + '}']
                        else:
                            raise Exception('不支持Array[bool] !')
        return inputs, cmd

    def get_outputs(self):
        outputs = []
        for name, v in self.cmd.outputs.items():
            outputs += [v.type + ' ' + name + ' = ' + f'"{v.path}"']
        return outputs

    def format_wdl_task(self, outfile, wdl_version='development'):
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
        data = self.cmd
        lines = f'task {data.meta.name}' + '{\n'
        # input
        lines += ' ' * 4 + 'input {\n'
        inputs, cmd = self.get_inputs_and_cmd()
        for line in inputs:
            lines += ' ' * 8 + line + '\n'
        lines += ' ' * 8 + '# for runtime\n'
        for k, v in self.get_runtime().items():
            if k in ['cpu', 'time_minutes']:
                typ = 'Int'
            else:
                typ = 'String'
                v = f'"{v}"'
            lines += ' ' * 8 + f'{typ} {k} = {v}' + '\n'
        lines += ' ' * 4 + '}\n\n'

        # command
        lines += ' ' * 4 + 'command <<<\n'
        lines += ' ' * 8 + 'set -e \n'
        base_cmd = ''
        tool_dir = data.runtime.tool_dir
        if tool_dir and not tool_dir.endswith('/'):
            tool_dir += '/'
        base_cmd += tool_dir + data.runtime.tool
        lines += ' ' * 8 + base_cmd + ' \\\n'
        for line in cmd:
            lines += ' ' * 8 + line + ' \\\n'
        lines = lines[:-2] + '\n'
        lines += ' ' * 4 + '>>>\n\n'

        # output
        lines += ' ' * 4 + 'output {\n'
        for line in self.get_outputs():
            lines += ' ' * 8 + line + '\n'
        lines += ' ' * 4 + '}\n\n'

        # runtime
        lines += ' ' * 4 + 'runtime {\n'
        for k, v in self.get_runtime().items():
            # lines += ' ' * 8 + f'{k}: "{v}"' + '\n'
            lines += ' ' * 8 + f'{k}: {k}' + '\n'
        lines += ' ' * 4 + '}\n\n'

        # task meta
        lines += ' ' * 4 + 'meta {\n'
        for k, v in self.get_task_meta().items():
            lines += ' ' * 8 + f'{k}: "{v}"' + '\n'
        lines += ' ' * 4 + '}\n\n'

        # parameter_meta
        lines += ' ' * 4 + 'parameter_meta {\n'
        for k, v in self.get_parameter_meta().items():
            v_info = '{' + ', '.join(f'{k}: "{v}"' for k, v in v.items()) + '}'
            lines += ' ' * 8 + f'{k}: {v_info}' + '\n'
        lines += ' ' * 4 + '}\n\n'

        # write wdl file
        lines += '}\n'
        with open(outfile, 'w', encoding='utf-8') as f:
            f.write(f'version {wdl_version}\n\n')
            f.write(lines)


# ---example---
def salmon():
    # 定义一个command
    cmd = Command()
    cmd.meta.name = 'salmon'
    cmd.meta.desc = 'transcript expression quantification'
    cmd.runtime.image = None
    cmd.runtime.memory = 1024
    cmd.runtime.cpu = 2
    cmd.runtime.tool_dir = '/usr/bin/'
    cmd.runtime.tool = 'salmon quant'
    cmd.args['indexDir'] = Argument(prefix='-i ', type='indir', level='required')
    cmd.args['read1'] = Argument(prefix='-1 ', type='infile', level='required')
    cmd.args['read2'] = Argument(prefix='-2 ', type='infile', level='required')
    cmd.args['outDir'] = Argument(prefix='-o ', type='str', level='optional')
    cmd.args['gcBias'] = Argument(prefix='--gcBias', type='bool', default=True)
    cmd.outputs['trans_expr'] = Output(path="~{outdir}" + "/quant.sf", locate='quant')
    return cmd


if __name__ == '__main__':
    wf = Workflow()
    wf.meta.name = 'pipeline'
    wf.meta.desc = 'rna-seq pipeline'
    depend = ['ddd']
    for sample in ['s1', 's2']:
        task_name = 'quant_' + sample
        task = Task(cmd=salmon(), depend=depend)
        task.cmd.format_wdl_task(outfile='salmon.wdl')
        task.cmd.args['read1'].value = sample + '.r1.fq'
        task.cmd.args['read2'].value = sample + '.r2.fq'
        task.cmd.args['outDir'].value = sample
        task.cmd.args['indexDir'].value = 'index/'
        print(task.cmd.format_cmd())
        print('task_id:', task.task_id)
        wf.tasks[task_name] = task
    print(wf)
