# coding=utf-8
import os
import subprocess
import glob
import time
import re


def run_cmd(cmd):
    return subprocess.call(cmd, shell=True)


def run_pipeline(fastq_dir, batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE',
                 out_dir='pipeline_result', exclude_steps=None, is_continue=False):
    sample_id_list = [os.path.basename(x).split('_')[0] for x in glob.glob(fastq_dir+'/*_R1_*.fastq*')]
    if not sample_id_list:
        raise Exception("find no fastq file!")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    cfg = '{}/simulate_project.cfg'.format(out_dir)
    if not is_continue:
        with open(cfg, 'w') as fw:
            fw.write('<project>\n')
            fw.write('fastq_dir = {}\n'.format(fastq_dir))
            fw.write('sample_id = {}\n'.format(','.join(sample_id_list)))
            fw.write('sample_name = {}\n'.format(','.join(sample_id_list)))
            fw.write('sample_cfg = {}-{}\n'.format(batch_id, panel))
            fw.write('sample_type = {}\n'.format(sample_type))
            fw.write('experiment_id = simulate_project\n')
            fw.write('</project>\n')
    # cmd = 'python /share/home/test.rd/app/rdscripts_release/pipeline/run_pipeline.py '
    cmd = 'python /share/home/test.rd/dev/Test.scripts/pipeline/run_pipeline.py '
    # cmd = 'python /share/home/test.rd/dev/Cscripts/pipeline/run_pipeline.py '
    # cmd += '-c {} '.format('/share/home/test.rd/app/rdscripts_release/configure_file/cfg.template.rd ')
    cmd += '-c /share/home/test.rd/dev/Test.scripts/configure_file/cfg.template '
    if exclude_steps:
        cmd += '-e {} '.format(exclude_steps)
    if is_continue:
        cmd += "--iscontinue "
    cmd += '-p {} '.format(cfg)
    cmd += '-o {} '.format(out_dir)
    run_cmd(cmd)
    # return sample_id_list
    return out_dir


def after_simulate_analysis(rawdata_dir_pattern="*_af*[0-9]",
                            batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE',
                            exclude_steps="qc,fusion,sv,cnv,msi,extend",
                            is_continue=False,
                            sleep=10000,
                            submit_size=5,
                            output_sample_dir="/share/home/deqing.gu/devdata/pipeline_output_Sample",):
    """
     批量分析模拟数据, 作为一键化模拟分析的补充脚本，提供对模拟数据进行续跑的分析选择和选择跳过一些步骤的选择。
    :param rawdata_dir_pattern: 字符串匹配模式，用于匹配模拟数据fastq所在文件夹名称，默认'*_af*[0-9]'
    :param batch_id: batch_id， 默认LK291， ctDNA对应LK292
    :param panel: panel，默认'Oncoscreen.520.v3'
    :param sample_type: 默认FFPE，ctDNA对应'ctDNA'
    :param exclude_steps: 'qc,fusion,sv,cnv,msi,extend'
    :param is_continue: 默认False
    :param sleep: 每次投递一批任务的间隔时间，默认10000秒，如果数据量大，请考虑设置更大的时间间隔以免占用太多服务器资源
    :param submit_size: 每批任务包含多少次pipeline分析，默认为5，如果数据量大，请考虑设置小一点
    :param output_sample_dir: pipeline分析结果的模板文件夹，因为跳过一些步骤后，为了最终report不报错，需把目录补齐
    :return:
    """
    out_dir_set = glob.glob(rawdata_dir_pattern)
    if not out_dir_set:
        raise Exception("Find No matched rawdata dir")

    for each in out_dir_set:
        fq_list = glob.glob(each + '/*_out_reads/*.fastq*')
        if not fq_list:
            fq_list = glob.glob(each+'/*/*.fastq*')
        if not os.path.exists(each+'/simulated_fastqs'):
            os.mkdir(each+'/simulated_fastqs')
            for fq in fq_list:
                os.system("mv {} {}/simulated_fastqs/".format(fq, each))
    out_dir_list = list(out_dir_set)
    for i in range(0, len(out_dir_list), submit_size):
        # 每次批量投递5个项目的分析
        batch = out_dir_list[i:i+submit_size]
        for each in batch:
            tmp_sample_dir = each+'/simulated_fastqs'
            samples = [os.path.basename(x).split('_')[0] for x in glob.glob(tmp_sample_dir + '/*_R1_*.fastq*')]
            if not is_continue:
                for sample in samples:
                    tmp_path = each+'_runPipeline'+'/'+sample
                    os.system('mkdir -p {}'.format(tmp_path))
                    for dir_name in ['CNV', 'QC', 'SV', 'fusion', 'MSI']:
                        os.system('cp -r {}/{} {}'.format(output_sample_dir, dir_name, tmp_path))
                        os.system('rename simName {} {}/{}/*'.format(sample, tmp_path, dir_name))

            run_pipeline(tmp_sample_dir, batch_id=batch_id, panel=panel, sample_type=sample_type,
                         out_dir=each+'_runPipeline', exclude_steps=exclude_steps, is_continue=is_continue)
        # 一定时间后再投递另外一批分析
        if not (i+submit_size >= len(out_dir_list)):
            time.sleep(sleep)
        else:
            print("全部任务投递完成！")


def introduce_command(func):
    import argparse
    import inspect
    import json
    import time
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
    with open("Argument_detail.json", 'w') as f:
        json.dump(args, f, indent=2, sort_keys=True)
    start = time.time()
    func(**args)
    print("total time: {}s".format(time.time() - start))


if __name__ == '__main__':
    introduce_command(after_simulate_analysis)
