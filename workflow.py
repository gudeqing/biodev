# coding=utf-8
import os
import sys
from pyflow import WorkflowRunner
import subprocess
import logging
import time


class WorkFlow(WorkflowRunner):
    def __init__(self, cmd_dict, dependency_dict, cmd_mem_dict=None, cmd_cpu_dict=None):
        """
        :param cmd_dict: A dict such as {"task_name2": cmd2, "task_nam3": cmd3， "task_name4": cmd4}
        :param dependency_dict: A dict likes
            {"task_name2": "task_name3, task_name4", task_name3:"", "task_name4": ""},
            in which key is a task, while value of key is the task(s) that key depends.
            For initial task, its dependency is an empty string.
        :param cmd_mem_dict: A dict which describes each task's memory requirement.
            The memory unit is 'Mb', such as {'task_name': 2048}. Default: 2048*3
        :param cmd_cpu_dict: A dict which describes each task's cpu number needed.
            Key is task_name, value is cpu number. Default: 1
        """
        if not cmd_dict:
            raise Exception("cmd_dict is empty!")
        if not dependency_dict:
            raise Exception("cmd_dict is empty!")
        for task, depend_tasks in dependency_dict.items():
            if depend_tasks:
                all_tasks = depend_tasks.strip().split(',') + [task]
            else:
                all_tasks = [task]
            for each in all_tasks:
                if each not in cmd_dict:
                    raise Exception("{} not in cmd_dict".format(each))
        self.cmd_dict = cmd_dict
        self.dependency_dict = dependency_dict
        if cmd_mem_dict is None:
            cmd_mem_dict = {x: 2048*0.5 for x in self.cmd_dict}
        else:
            for each in self.cmd_dict:
                if each not in cmd_mem_dict:
                    cmd_mem_dict[each] = 2048
        if cmd_cpu_dict is None:
            cmd_cpu_dict = {x : 1 for x in self.cmd_dict}
        else:
            for each in self.cmd_dict:
                if each not in cmd_cpu_dict:
                    cmd_cpu_dict[each] = 1
        self.cmd_mem_dict = cmd_mem_dict
        self.cmd_cpu_dict = cmd_cpu_dict
        self.task_add_order = list()

    def workflow(self):
        added_tasks = set()
        for task, depend_tasks in self.dependency_dict.items():
            if not depend_tasks:
                self.addTask(task, self.cmd_dict[task], memMb=self.cmd_mem_dict[task], nCores=self.cmd_cpu_dict[task])
                self.dependency_dict.pop(task)
                added_tasks.add(task)
                self.task_add_order.append(task)
        while self.dependency_dict:
            the_circle_did_add_task = False
            for task, depend_tasks in self.dependency_dict.items():
                depend_task_list = depend_tasks.strip().split(',')
                if all(x in added_tasks for x in depend_task_list):
                    self.addTask(task, self.cmd_dict[task], memMb=self.cmd_mem_dict[task],
                                 nCores=self.cmd_cpu_dict[task], dependencies=depend_task_list)
                    self.dependency_dict.pop(task)
                    added_tasks.add(task)
                    self.task_add_order.append(task)
                    the_circle_did_add_task = True
            if not the_circle_did_add_task:
                not_added_tasks = self.dependency_dict.keys()
                error = ";".join(not_added_tasks), "--> At least one of the task is not connected with others!"
                raise Exception(error)


def test_workflow(cmd_dict=None, dependency_dict=None):
    if cmd_dict is None and dependency_dict is None:
        """
        --A
            } -------->D     
        --B              } ---> G
            } --> E -->F
        --C         
        """
        cmd_dict = {
            'A': 'echo "task:A"',
            'B': 'echo "task:B"',
            'C': 'echo "task:C"',
            'D': 'echo "task:D"',
            'E': 'echo "task:E"',
            'F': 'echo "task:F"',
            'G': 'echo "task:G"',
        }
        dependency_dict = {
            'A': '',
            'B': '',
            'C': '',
            'D': 'A,B',
            'E': 'B,C',
            'F': 'E',
            'G': 'D,F',
        }
    wk = WorkFlow(cmd_dict, dependency_dict)
    wk.run()
    print(wk.task_add_order)


def get_logger(output='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    # format log
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # write log
    fh = logging.FileHandler(output, mode='w+')
    fh.setFormatter(formatter)
    logger.setLevel(logging.INFO)
    logger.addHandler(fh)
    # print log
    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(formatter)
    sh.setLevel(logging.WARN)
    logger.addHandler(sh)
    return logger


def run_cmd(cmd):
    return subprocess.check_call(cmd, shell=True)


def make_dirs(path):
    os.system('mkdir -p {}'.format(path))
    return path


def introduce_command(func):
    import argparse
    import inspect
    import json

    def wrapper():
        parser = argparse.ArgumentParser(description=func.__doc__, formatter_class=argparse.RawTextHelpFormatter)
        func_args = inspect.getargspec(func)
        arg_names = func_args.args
        arg_defaults = func_args.defaults
        arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
        for arg, value in zip(arg_names, arg_defaults):
            if value == 'None':
                parser.add_argument('-'+arg, required=True, metavar=arg)
            elif type(value) == bool:
                if value:
                    parser.add_argument('--'+arg, action="store_false", help='default: True')
                else:
                    parser.add_argument('--'+arg, action="store_true", help='default: False')
            elif value is None:
                parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
            else:
                parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
        if func_args.varargs is not None:
            print("warning: *varargs is not supported, and will be neglected! ")
        if func_args.keywords is not None:
            print("warning: **keywords args is not supported, and will be neglected! ")
        args = parser.parse_args().__dict__
        with open("Argument_detail_for_{}.json".format(os.path.basename(__file__).split(".")[0]), 'w') as f:
            json.dump(args, f, indent=2, sort_keys=True)
        start = time.time()
        func(**args)
        print("total time: {}s".format(time.time() - start))
    return wrapper


if __name__ == '__main__':
    test_workflow()

