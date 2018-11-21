# coding=utf-8
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pyflow import WorkflowRunner
import glob
import time
import logging
import random
import linecache
import re
import argparse

# 对于一种类型样本（如tissue），模拟要求如下：
# 第一步：基于源测序数据，按照1000x的要求抽出5个样本数据，seed分别为1，2，3，4，5
# 第二步：对第一步得到的5个样本数据，完成质控分析和比对分析，得到5个bam文件
# 数据模拟关键参数：（1）bam文件 （2）AF值 （3）bed文件（包含突变信息）（4）seed值
# 第三步：突变模拟，涉及组合：5bam5seed * 4AF * 20bed = 400 次模拟，其中20bed = 10snv + 10indel，4AF=0.01,0.05,0.10,0.20
# 总结：5次downsample，5次比对分析，400次模拟。


class WorkFlow(WorkflowRunner):
    def __init__(self, cmd_dict, dependency_dict, cmd_mem_dict=None, cmd_cpu_dict=None):
        """
        :param cmd_dict: A dict such as {"task_name2": cmd2, "task_nam3": cmd3， "task_name4": cmd4}
        :param dependency_dict: A dict likes {"task_name2": "task_name3, task_name4", task_name3:"", "task_name4": ""},
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
            cmd_mem_dict = {x: 2048 for x in self.cmd_dict}
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


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w+')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger


def run_cmd(cmd):
    return subprocess.call(cmd, shell=True)


def sampling(sample_name, raw_dir, qc_result, depth=3000, seeds=(1,3,5,7,9), out_dir=None):
    "sampling a sample with `len(seeds)` replicates"
    cmd = "python /share/home/test.rd/app/rdscripts_release/tools/downsample2.py "
    cmd += '--rawdata {} '.format(raw_dir)
    cmd += '--qc_result {} '.format(qc_result)
    cmd += '--number {} '.format(depth)
    cmd_list = list()
    out_dirs = list()
    for seed in seeds:
        rep_sample_name = sample_name + '_Seed{}'.format(seed)
        out_dirs.append(rep_sample_name)
        if not os.path.exists(rep_sample_name):
            os.mkdir(rep_sample_name)
        new_cmd = cmd + '--outdir {} '.format(rep_sample_name)
        new_cmd += '--seed {} '.format(seed)
        cmd_list.append(new_cmd)
    with ThreadPoolExecutor(6) as pool:
        pool.map(run_cmd, cmd_list)
    if not out_dir:
        if not os.path.exists("downsampled_fastqs"):
            os.mkdir("downsampled_fastqs")
        out_dir = os.path.abspath("downsampled_fastqs")
    for each in out_dirs:
        fastqs = glob.glob('{}/*fastq*'.format(each))
        for fq in fastqs:
            sn = each.split('_')[-1] + os.path.basename(fq)
            status = os.system("mv {} {}/{}".format(fq, out_dir, sn))
            if int(status) != 0:
                raise Exception("failed to mv {}".format(fq))
    fq_gz = glob.glob(out_dir+'/*gz')
    for each in fq_gz:
        os.system("gzip -d {} ".format(each))
    return os.path.abspath(out_dir)


def run_pipeline(raw_dir, batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE',
                 out_dir='pipeline_result', exclude_steps=None):
    sample_id_list = [os.path.basename(x).split('_')[0] for x in glob.glob(raw_dir+'/*_R1_*.fastq*')]
    if not sample_id_list:
        raise Exception("find no fastq file!")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    cfg = '{}/simulate_project.cfg'.format(out_dir)
    with open(cfg, 'w') as fw:
        fw.write('<project>\n')
        fw.write('fastq_dir = {}\n'.format(raw_dir))
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
        cmd += '-e {} '.format(','.join(exclude_steps))
    cmd += '-p {} '.format(cfg)
    cmd += '-o {} '.format(out_dir)
    run_cmd(cmd)
    # return sample_id_list
    return out_dir


def wait_bam_born(bam_dir):
    time.sleep(3)
    while 1:
        time.sleep(10)
        bam_bai_file = glob.glob(bam_dir+'/*RAW.*bam.bai')
        if bam_bai_file:
            bam_bai_file = bam_bai_file[0]
            bai_size = os.path.getsize(bam_bai_file)
            time.sleep(99)
            bai_size_latest = os.path.getsize(bam_bai_file)
            if bai_size_latest == bai_size:
                break
        else:
            continue


def filter_vcf(raw_vcfs, exclude_vcfs, out_vcf=None, sep='\t'):
    """
     filter vcf based on provided vcf
    :param raw_vcf: raw vcf
    :param exclude_vcf: to be filtered vcf
    :param out_vcf:  out file
    :return:
    """
    logger = set_logger('filter_vcf.log', 'filter_vcf')
    if type(exclude_vcfs) == str:
        exclude_vcfs = [exclude_vcfs]
    else:
        assert type(exclude_vcfs) == list
    if type(raw_vcfs) == str:
        raw_vcfs = [raw_vcfs]
    else:
        assert type(raw_vcfs) == list

    raw_vcf_set = set()
    for each in raw_vcfs:
        with open(each) as fr:
            lines = (x.strip().split(sep)[0:5] for x in fr if not x.startswith('#') and x.strip())
            for line in lines:
                # chr_, pos, name, ref, alt = line
                line.pop(2)
                raw_vcf_set.add('_'.join(line).lower())
    if not raw_vcf_set:
        logger.error("No mutation was found in the provided vcf file !")
        raise Exception("No mutation was found in the provided vcf file !")

    for each in exclude_vcfs:
        with open(each) as fr:
            lines = (x.strip().split(sep)[0:5] for x in fr if not x.startswith('#') and x.strip())
            for line in lines:
                # chr_, pos, name, ref, alt = line
                line.pop(2)
                mutation = '_'.join(line).lower()
                if mutation.split('chr', 1)[1] in raw_vcf_set or mutation in raw_vcf_set:
                    logger.info('mutation: ' + mutation + 'is excluded for it is included in {}'.format(each))
                    raw_vcf_set.remove(mutation)
    if not raw_vcf_set:
        logger.error("All mutations were excluded by referring input exclude vcf files")
        raise Exception("All mutations were excluded by referring input exclude vcf files")
    if out_vcf is None:
        if not os.path.exists('filtered_mutation'):
            os.mkdir('filtered_mutation')
        out_vcf = 'filtered_mutation/filtered.vcf'
    logger.info("{} mutations were selected".format(len(raw_vcf_set)))
    fw = open(out_vcf, 'w')
    for each in raw_vcfs:
        with open(each) as fr:
            lines = (x for x in fr if not x.startswith('#') and x.strip())
            for line in lines:
                line_list = line.strip().split(sep)[0:5]
                line_list.pop(2)
                if '_'.join(line_list).lower() in raw_vcf_set:
                    fw.write(line)
    fw.close()
    return os.path.abspath(out_vcf)


def random_read(in_file):
    # 随机读取文件的行，且保证每一行只读一次
    count = len(open(in_file, 'rU').readlines())  # 返回文件行数
    s = range(1, count + 1)  # 生成列表，从1到文件行长度
    # 洗牌，次数为文件长度，每次交换两个位置
    for i in s:
        s1 = random.randint(0, count - 1)
        s2 = random.randint(0, count - 1)
        s[s1], s[s2] = s[s2], s[s1]  # 根据生成的随机数交换位置
    # 读取文件内容并显示
    for i in s:
        yield i, linecache.getline(in_file, i) # 从指定文件读取指定行


def split_target_mutation_vcf(target_vcf, size=100, drop_smaller_size=True, random_mode=True, need_file_num=10):
    logger = set_logger('split_target_mutation_vcf.log', 'x1')
    split_result = list()
    fw = open(target_vcf+'_'+'0', 'w')
    eff_num = -1
    if random_mode:
        fr = random_read(target_vcf)
    else:
        fr = enumerate(open(target_vcf))
    for line_num, line in fr:
        chr_, pos, name, ref, alt = line.strip().split("\t")[:5]
        if 'n' in ref.lower() or 'n' in alt.lower():
            logger.warning('line:{} has invalid base "N/n" in ref/alt in {}'.format(line_num, target_vcf))
            continue
        else:
            eff_num += 1
        if eff_num%int(size) == 0:
            fw.close()
            if len(split_result) == need_file_num:
                eff_num -= 1
                break
            file_name = target_vcf+'_'+str(eff_num//int(size))
            fw = open(file_name, 'w')
            split_result.append(os.path.abspath(file_name))
        fw.write(line)
    fw.close()
    if not random_mode:
        fr.close()
    if drop_smaller_size and (eff_num+1)%(int(size)) != 0:
        logger.warning("{} was dropped for not enough mutation number!".format(split_result.pop()))
    if len(split_result) < need_file_num:
        logger.info('split file number: {}'.format(len(split_result)))
        logger.info("split file number is smaller than needed")
    return split_result


def simulate_mutation_cmd(vcf_path, bam_path, fq1, fq2, af=0.2, mu_type='snv', seed=1):
    if mu_type.lower() == 'snv':
        cmd = "python /share/home/deqing.gu/simulation/Simulation/snv_scripts/run_snv_simu.py "
    elif mu_type.lower() == 'indel':
        cmd = "python /share/home/deqing.gu/simulation/Simulation/indel_scripts/run_indel_simu.py "
    elif mu_type.lower() == 'fusion':
        cmd = "python /share/home/deqing.gu/simulation/Simulation/fusion_scripts/run_fusion_simu_mp.py "
    elif mu_type.lower() == 'cnv':
        cmd = "python /share/home/test.rd/dev/SimuVar/cnv_scripts/run_cnv_simu_v4.py "
        out_dir = mu_type + vcf_path.split('_')[-1] + '_af' + str(af).replace('.', '')
        if not os.path.exists(out_dir):
            os.system("mkdir -p {} ".format(out_dir))
        new_sample_name = os.path.basename(bam_path)[:-8]
        new_sample_name += mu_type.capitalize()
        new_sample_name += os.path.basename(vcf_path).split('_')[-1].replace('.', '_')
        new_sample_name += 'af' + str(af).replace('.', '') + 'seed' + str(seed)
        cmd += "-i {} ".format(vcf_path)
        cmd += "-b {} ".format(bam_path)
        cmd += "--sam {} ".format(bam_path+'.2.sam')
        cmd += "-r {} ".format("/share/c400/kai.feng/test/bam_surgeon/ref_genome/human_g1k_v37.fasta")
        cmd += "-s {} ".format(seed)
        cmd += "-o {} ".format(out_dir)
        cmd += "-p {} ".format(new_sample_name)
        cmd += "-f1 {} ".format(fq1)
        cmd += "-f2 {} ".format(fq2)
        return cmd, out_dir
    else:
        cmd = "？"
    out_dir = mu_type+vcf_path.split('_')[-1]+'_af'+str(af).replace('.', '')
    # new_sample_name = os.path.basename(bam_path)[:-8]+'_'+mu_type+vcf_path.split('_')[-1]+'_af'+str(int(af*100))+'_seed'+str(seed)
    new_sample_name = os.path.basename(bam_path)[:-8]
    new_sample_name += mu_type.capitalize()
    new_sample_name += os.path.basename(vcf_path).split('_')[-1].replace('.', '_')
    new_sample_name += 'af'+str(af).replace('.', '')+'seed'+str(seed)
    cmd += "-file {} ".format(vcf_path)
    cmd += "-bam {} ".format(bam_path)
    cmd += "-af {} ".format(af)
    cmd += "-seed {} ".format(seed)
    cmd += "-out {} ".format(out_dir)
    cmd += "-pre {} ".format(new_sample_name)
    cmd += "-f1 {} ".format(fq1)
    cmd += "-f2 {} ".format(fq2)
    return cmd, out_dir


def simulate_mutation_cmd2(vcf_path, bam_path, fq1, fq2, af=0.2, mu_type='snv', seed=1,
                           expect_mutation_number=100, depth_threshold=600):
    if mu_type.lower() == 'snv' or mu_type.lower() == 'indel':
        cmd = "python /share/home/deqing.gu/pycharm_projects/mutation_simulation/simulate_mutation_v2.py "
        new_sample_name = os.path.basename(bam_path)[:-8]
        new_sample_name += mu_type.capitalize()
        new_sample_name += os.path.basename(vcf_path).split('_')[-1].replace('.', '_')
        new_sample_name += 'af' + str(af).replace('.', '') + 'seed' + str(seed)
        out_dir = mu_type + vcf_path.split('_')[-1] + '_af' + str(af).replace('.', '')
        if not os.path.exists(out_dir):
            os.system("mkdir -p {} ".format(out_dir))
        cmd += "-mutation_file {} ".format(vcf_path)
        cmd += "-bam_file {} ".format(bam_path)
        cmd += "-af {} ".format(af)
        cmd += '-output_dir {} '.format(out_dir+'/'+new_sample_name)
        cmd += "-out_fastq_prefix {} ".format(new_sample_name)
        cmd += "-fq1 {} ".format(fq1)
        cmd += "-fq2 {} ".format(fq2)
        cmd += "-expect_mutation_number {} ".format(expect_mutation_number)
        cmd += "-depth_threshold {} ".format(depth_threshold)
        return cmd, out_dir

    if mu_type.lower() == 'fusion':
        cmd = "python /share/home/deqing.gu/simulation/Simulation/fusion_scripts/run_fusion_simu_mp.py "
        out_dir = mu_type + vcf_path.split('_')[-1] + '_af' + str(af).replace('.', '')
        if not os.path.exists(out_dir):
            os.system("mkdir -p {} ".format(out_dir))
        new_sample_name = os.path.basename(bam_path)[:-8]
        new_sample_name += mu_type.capitalize()
        new_sample_name += os.path.basename(vcf_path).split('_')[-1].replace('.', '_')
        new_sample_name += 'af' + str(af).replace('.', '') + 'seed' + str(seed)
        cmd += "-file {} ".format(vcf_path)
        cmd += "-bam {} ".format(bam_path)
        cmd += "-af {} ".format(af)
        cmd += "-seed {} ".format(seed)
        cmd += "-out {} ".format(out_dir)
        cmd += "-pre {} ".format(new_sample_name)
        cmd += "-f1 {} ".format(fq1)
        cmd += "-f2 {} ".format(fq2)
        return cmd, out_dir

    if mu_type.lower() == 'cnv':
        cmd = "python /share/home/test.rd/dev/SimuVar/cnv_scripts/run_cnv_simu_v4.py "
        out_dir = mu_type + vcf_path.split('_')[-1] + '_af' + str(af).replace('.', '')
        if not os.path.exists(out_dir):
            os.system("mkdir -p {} ".format(out_dir))
        new_sample_name = os.path.basename(bam_path)[:-8]
        new_sample_name += mu_type.capitalize()
        new_sample_name += os.path.basename(vcf_path).split('_')[-1].replace('.', '')
        new_sample_name += 'af' + str(af).replace('.', '') + 'seed' + str(seed)
        cmd += "-i {} ".format(vcf_path)
        cmd += "-b {} ".format(bam_path)
        cmd += "--sam {} ".format(bam_path+'.2.sam')
        cmd += "-r {} ".format("/share/c400/kai.feng/test/bam_surgeon/ref_genome/human_g1k_v37.fasta")
        cmd += "-s {} ".format(seed)
        cmd += "-o {} ".format(out_dir)
        cmd += "-p {} ".format(new_sample_name)
        cmd += "-f1 {} ".format(fq1)
        cmd += "-f2 {} ".format(fq2)
        return cmd, out_dir


def batch_simulate(sampled_dir, pipeline_out='pipeline_result', vcf_list=None, pyflow=True,
                   mu_type='snv', af_list=(0.01, 0.02, 0.05, 0.1, 0.2),
                   expect_mutation_number=100, depth_threshold=600):
    time_mark = str(time.time()).replace('.', '_')
    logger = set_logger('batch_simulation_{}.log'.format(time_mark), 'x2')
    sample_id_list = [os.path.basename(x).split('_')[0] for x in glob.glob(sampled_dir + '/*_R1_*.fastq*')]
    # vcf_list = split_target_mutation_vcf(target_vcf, size=chunk_vcf_size)
    cmd_list = list()
    task_ids = list()
    out_dir_set = set()
    for seed, sample in enumerate(sample_id_list):
        bam_path = pipeline_out + '/{sample}/bam/{sample}.RAW.bam'.format(sample=sample)
        if mu_type == 'cnv':
            os.system('samtools view -@ 4 -O SAM -o {out} {bam}'.format(out=bam_path+'.2.sam', bam=bam_path))
        fq1 = glob.glob(sampled_dir+'/{}*_R1_*fastq*'.format(sample))[0]
        fq2 = glob.glob(sampled_dir+'/{}*_R2_*fastq*'.format(sample))[0]
        for vcf in vcf_list:
            for af in af_list:
                try:
                    # 这样可以保证sample的seed和simulate的seed用的一样
                    seed = int(re.match('Seed(\d+).*', sample).groups()[0])
                except Exception as e:
                    logger.warning('Failed to grab seed from sample name, but this is not important!')
                    seed = seed*2 + 1
                task_ids.append(sample+mu_type.capitalize()+os.path.basename(vcf).split('_')[-1].replace('.', '')+'af'+str(af).replace('.', '')+'seed'+str(seed))
                cmd, out_dir = simulate_mutation_cmd2(vcf, bam_path, fq1, fq2, af=af, mu_type=mu_type, seed=seed,
                                                      expect_mutation_number=expect_mutation_number,
                                                      depth_threshold=depth_threshold)
                logger.warning(cmd)
                cmd_list.append(cmd)
                out_dir_set.add(out_dir)
    if not pyflow:
        with ThreadPoolExecutor(20) as pool:
            pool.map(run_cmd, cmd_list)
        logger.warning('simulation finished')
    else:
        step_len = 50
        for i in range(0, len(cmd_list), step_len):
            tmp_task_ids = task_ids[i: i+step_len]
            tmp_cmd_list = cmd_list[i: i+step_len]
            cmd_dict = dict(zip(tmp_task_ids, tmp_cmd_list))
            dependency_dict = {x: '' for x in tmp_task_ids}
            time_mark = str(time.time()).replace('.', '_')
            status = WorkFlow(cmd_dict, dependency_dict).run(mode='sge', dataDirRoot='run_simulate_{}.pyflow'.format(time_mark))
            if status == 0:
                logger.warning("simulation finished")
            else:
                raise Exception("Failed to run batch simulation")
    return out_dir_set


def after_simulate_analysis(out_dir_set=None,  batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE',
                            output_sample_dir="/share/home/deqing.gu/devdata/pipeline_output_Sample"):
    if out_dir_set is None:
        out_dir_set = glob.glob('*[0-9]_af*[0-9]')
    for each in out_dir_set:
        fq_list = glob.glob(each + '/*_out_reads/*.fastq*')
        if not fq_list:
            fq_list = glob.glob(each+'/*/*.fastq*')
        if not os.path.exists(each+'/simulated_fastqs'):
            os.mkdir(each+'/simulated_fastqs')
            for fq in fq_list:
                os.system("mv {} {}/simulated_fastqs/".format(fq, each))
    out_dir_list = list(out_dir_set)
    submit_size = 5
    if sample_type.lower() == 'ctdna':
        submit_size = 3
    for i in range(0, len(out_dir_list), submit_size):
        # 每次批量投递5个项目的分析
        batch = out_dir_list[i:i+submit_size]
        for each in batch:
            tmp_sample_dir = each+'/simulated_fastqs'
            samples = [os.path.basename(x).split('_')[0] for x in glob.glob(tmp_sample_dir + '/*_R1_*.fastq*')]
            for sample in samples:
                tmp_path = each+'_runPipeline'+'/'+sample
                os.system('mkdir -p {}'.format(tmp_path))
                for dir_name in ['CNV', 'QC', 'SV', 'fusion', 'MSI']:
                    os.system('cp -r {}/{} {}'.format(output_sample_dir, dir_name, tmp_path))
                    os.system('rename simName {} {}/{}/*'.format(sample, tmp_path, dir_name))

            run_pipeline(tmp_sample_dir, batch_id=batch_id, panel=panel, sample_type=sample_type,
                         out_dir=each+'_runPipeline',
                         exclude_steps=['qc', 'fusion', 'sv', 'cnv', 'msi', 'extend'])
        # 一定时间后再投递另外一批分析
        time.sleep(8888+1111)
        if sample_type.lower() == 'ctdna':
            time.sleep(9999+6666)


def main(sample_name=None, raw_dir=None, qc_result=None, depth=1000, seeds="1,3,5,7,9",
         batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE',
         downsample_dir=None,
         downsample_pipeline_out=None,
         target_vcfs=None, chunk_vcf_size=120, chunk_number=10, exclude_vcfs=None,
         mu_type='snv', af_list="0.01,0.02,0.05,0.1,0.2", expect_mutation_number=100, depth_threshold=600,
         pyflow=True, analysis_simulated=True,
         output_sample_dir="/share/home/deqing.gu/devdata/pipeline_output_Sample"):
    """
    downsample+模拟+分析
    :param sample_name: 被downsample的样本id，与raw_dir目录中的信息应该一致，如果不提供，则需提供downsample_dir
    :param raw_dir: 被downsample的样本的fastq所在目录
    :param qc_result: 被downsample的样本的QC分析结果
    :param depth: 要模拟的测序深度，默认1000
    :param seeds: downsample时的seed参数，默认"1,3,5,7,9"，表示要进行5次downsample
    :param batch_id: 一键化分析需要的参数，用于生成config文件，默认'LK291'
    :param panel: 一键化分析需要的参数，用于生成config文件，默认'Oncoscreen.520.v3'
    :param sample_type:  一键化分析需要的参数，用于生成config文件，默认'FFPE'
    :param downsample_dir: downsample的结果目录
    :param downsample_pipeline_out: downsample数据的一键化分析结果目录，如果提供，则可以跳过该一键化分析步骤。
    :param target_vcfs: 待模拟的vcf，包含要模拟的位点信息，请提供比预期还要多的模拟位点（因为还要过滤），如不提供，则不进行模拟分析。
    :param chunk_vcf_size: 对待模拟的vcf文件进行切割的大小，如100，表示每份包含100个突变，
    :param chunk_number: 对待模拟的vcf文件进行切割的次数，如10，表示切割10份，这10份将分别用于模拟。
    :param exclude_vcfs: 一般输入已知的germline的突变信息，包含在这里的突变位点将不被模拟。
    :param mu_type: 模拟的突变类型，为'snv'或'indel'
    :param af_list: 模拟的频率，默认"0.01,0.02,0.05,0.1,0.2"，最终模拟的样本数=seeds长度*chunk_number*af_list长度
    :param expect_mutation_number:  决定引入的突变个数，即期望从每份vcf中提取的模拟位点个数，多余的位点信息将不被考虑
    :param depth_threshold: 挑选模拟位点时用到的覆盖度阈值
    :param pyflow: 是否使用pyflow管理模拟进程，默认使用，如设置则不使用pyflow，且该脚本需投递执行。
    :param analysis_simulated: 是否对模拟数据进行分析，默认为True，如设置，则不分析
    :param output_sample_dir: pipeline分析结果的模板文件夹，默认跳过一些步骤，为了最终report不报错，需把目录补齐
    :return:
    使用示例：
    python simulate.py
    -sample_name RD1804206FFP
    -raw_dir /share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/FASTQ
    -qc_result /share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/QC/RD1804206FFP.QC.xls
    -seeds 1,2
    -af_list 0.02,0.05
    -target_vcfs /share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/snv.vcf
    -chunk_vcf_size 130
    -chunk_number 1
    -exclude_vcfs /share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/NA12878.vcf
    -mu_type snv
    -expect_mutation_number 100
    -depth_threshold 600
    """
    if sample_name:
        if not raw_dir:
            print("Please provide rawdata directory for down-sampling fastq!")
            exit(0)
        if not qc_result:
            print('Please provide qc result file for down-sampling fastq')
            exit(0)
        seeds = [int(x) for x in  seeds.strip().split(",")]
        downsample_dir = sampling(sample_name, raw_dir, qc_result, depth=depth, seeds=seeds)
    else:
        print("你没有提供样本名，所以跳过downsample阶段")
        if not os.path.exists(downsample_dir) or (downsample_dir is None):
            print("downsample result directory is not existed!")
            exit(0)

    if downsample_pipeline_out:
        if not os.path.exists(downsample_pipeline_out):
            print("你提供的downsample数据一键化分析结果的路径不对")
        else:
            print("你提供了downsample数据一键化分析结果的路径，故跳过对downsample数据的分析,直接进入模拟阶段")
    else:
        downsample_pipeline_out = 'pipeline_result'
        run_pipeline(downsample_dir, batch_id=batch_id, panel=panel, sample_type=sample_type,
                     out_dir=downsample_pipeline_out)
        time.sleep(3666)  # wait for run_pipeline
        sample_id_list = [os.path.basename(x).split('_')[0] for x in glob.glob(downsample_dir + '/*_R1_*.fastq*')]
        bam_dir_list = (downsample_pipeline_out + '/{}/bam/'.format(x) for x in sample_id_list )
        map(wait_bam_born, bam_dir_list)  # wait for bam.bai

    if not target_vcfs:
        print("你没有提供待模拟vcf文件，因此不做模拟分析，但仍然可以完成downsample工作和对downsample的一键化分析工作")
        return

    if not exclude_vcfs:
        print("你没有提供类似记录germline突变的vcf文件，因此不对你提供的待模拟vcf文件进行过滤")
    else:
        target_vcfs = filter_vcf(target_vcfs, exclude_vcfs)
    vcf_list = split_target_mutation_vcf(target_vcfs, size=chunk_vcf_size,
                                         drop_smaller_size=True, need_file_num=chunk_number)

    # simulation
    af_list = [float(x) for x in af_list.strip().split(',')]
    out_dir_set = batch_simulate(downsample_dir, pipeline_out=downsample_pipeline_out, vcf_list=vcf_list,
                                 mu_type=mu_type, af_list=af_list, pyflow=pyflow,
                                 expect_mutation_number=expect_mutation_number,
                                 depth_threshold=depth_threshold)
    if not analysis_simulated:
        print("你选择暂不对模拟数据进行一键化分析，后续你可以选择用analysis_simulated_data.py进行分析")
        return out_dir_set
    after_simulate_analysis(out_dir_set, batch_id=batch_id, panel=panel, sample_type=sample_type,
                            output_sample_dir=output_sample_dir)


def main2(sampled_dir=None, pipeline_out=None, batch_id='LK291', panel='Oncoscreen.520.v3',
         sample_type='FFPE', target_vcfs=None, chunk_vcf_size=100, chunk_number=10,
         exclude_vcfs=None, mu_type='snv', af_list=(0.01, 0.02, 0.05, 0.1, 0.2), pyflow=True, final_analysis=True):
    # sampled_dir = sampling(sample_name, raw_dir, qc_result, depth=depth, seeds=seeds)
    # sampled_dir = test_sampling()
    if not sampled_dir or not os.path.exists(sampled_dir):
        raise Exception('Please provide valid fastq dir')
    sample_id_list = [os.path.basename(x).split('_')[0] for x in glob.glob(sampled_dir + '/*_R1_*.fastq*')]
    if not pipeline_out or not os.path.exists(pipeline_out):
        pipeline_out  = "pipeline_result"
        run_pipeline(sampled_dir, batch_id=batch_id, panel=panel, sample_type=sample_type, out_dir=pipeline_out)
        time.sleep(4666)  # wait for run_pipeline
        bam_dir_list = (pipeline_out + '/{}/bam/'.format(x) for x in sample_id_list )
        map(wait_bam_born, bam_dir_list)  # wait for bam.bai
    target_vcf = filter_vcf(target_vcfs, exclude_vcfs)
    vcf_list = split_target_mutation_vcf(target_vcf, size=chunk_vcf_size, drop_smaller_size=True, need_file_num=chunk_number)
    out_dir_set = batch_simulate(sampled_dir, pipeline_out=pipeline_out, vcf_list=vcf_list,
                                 mu_type=mu_type, af_list=af_list, pyflow=pyflow)
    if not final_analysis:
        return out_dir_set
    after_simulate_analysis(out_dir_set, batch_id=batch_id, panel=panel, sample_type=sample_type)


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
    with open("Argument_detail_for_{}.json".format(os.path.basename(__file__).split(".")[0]), 'w') as f:
        json.dump(args, f, indent=2, sort_keys=True)
    start = time.time()
    func(**args)
    print("total time: {}s".format(time.time() - start))




# -----------------------下面是各种测试用例-----------------------------
# for FFPE simulate
def test_sampling():
    sample_name = "RD1804206FFP"
    raw_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/FASTQ"
    qc_result = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/QC/RD1804206FFP.QC.xls"
    depth = 1000
    seeds = [1,3,5]
    sampled_dir = sampling(sample_name, raw_dir, qc_result, depth=depth, seeds=seeds)
    sample_name = "RD1804207FFP"
    raw_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/FASTQ"
    qc_result = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/QC/RD1804207FFP.QC.xls"
    depth = 1000
    seeds = [11, 7, 9]
    sampled_dir = sampling(sample_name, raw_dir, qc_result, depth=depth, seeds=seeds, out_dir=sampled_dir)
    return sampled_dir


def test_snv_main3():
    # sample_name = "RD1804206FFP"
    # raw_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/FASTQ"
    # qc_result = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/QC/RD1804206FFP.QC.xls"
    exclude_vcfs = ['/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/NA12878.vcf', ]
    target_vcfs = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/snv.vcf"
    mu_type = 'snv'
    af_list = (0.01, 0.02, 0.05, 0.1, 0.2)
    sampled_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/downsampled_fastqs"
    pipeline_out = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/pipeline_result"
    os.system('ln -s {} .'.format(pipeline_out))
    main2(sampled_dir=sampled_dir, pipeline_out=pipeline_out,
         batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE',
         target_vcfs=target_vcfs, chunk_vcf_size=130, exclude_vcfs=exclude_vcfs,  # > split target vcf
         mu_type=mu_type, af_list=af_list,  final_analysis=True,  # > simulation
    )


def test_indel_main4():
    # sample_name = "RD1804206FFP"
    # raw_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/FASTQ"
    # qc_result = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/QC/RD1804206FFP.QC.xls"
    exclude_vcfs = ['/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/NA12878.vcf', ]
    target_vcfs = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/indel.vcf"
    mu_type = 'indel'
    af_list = (0.01, 0.02, 0.05, 0.1, 0.2)
    sampled_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/downsampled_fastqs"
    pipeline_out = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/pipeline_result"
    os.system('ln -s {} .'.format(pipeline_out))
    main2(sampled_dir=sampled_dir, pipeline_out=pipeline_out,
         batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE',
         target_vcfs=target_vcfs, chunk_vcf_size=130, exclude_vcfs=exclude_vcfs,  # > split target vcf
         mu_type=mu_type, af_list=af_list,  # > simulation
    )


# for ctDNA
def test_sampling_ctdna():
    sample_name = "RD1803971PLA"
    raw_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_ctDNA/fastqs"
    qc_result = "/share/home/deqing.gu/devdata/simulate_rawdata_for_ctDNA/RD1803971PLA.QC.xls"
    depth = 10000
    seeds = [1,3,5]
    sampled_dir = sampling(sample_name, raw_dir, qc_result, depth=depth, seeds=seeds)
    sample_name = "RD1803972PLA"
    raw_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_ctDNA/fastqs"
    qc_result = "/share/home/deqing.gu/devdata/simulate_rawdata_for_ctDNA/RD1803972PLA.QC.xls"
    depth = 10000
    seeds = [11, 7, 9]
    sampled_dir = sampling(sample_name, raw_dir, qc_result, depth=depth, seeds=seeds, out_dir=sampled_dir)
    return sampled_dir


def test_indel_main_ctdna():
    exclude_vcfs = ['/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/NA12878.vcf', ]
    target_vcfs = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/indel.vcf"
    mu_type = 'indel'
    af_list = (0.002, 0.005, 0.01, 0.02, 0.05)
    sampled_dir = "/share/home/deqing.gu/devdata/simulate_rawdata_for_ctDNA/downSample/downsampled_fastqs"
    main2(sampled_dir=sampled_dir, pipeline_out=None,
         batch_id='LK292', panel='Oncoscreen.520.v3', sample_type='ctDNA',
         target_vcfs=target_vcfs, chunk_vcf_size=130, exclude_vcfs=exclude_vcfs,  # > split target vcf
         mu_type=mu_type, af_list=af_list,  pyflow=True
    )


def test_snv_main_ctdna():
    exclude_vcfs = ['/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/NA12878.vcf', ]
    target_vcfs = "/share/home/deqing.gu/devdata/simulate_rawdata_for_FFPE/mutate_candidates/snv.vcf"
    mu_type = 'snv'
    af_list = (0.002, 0.005, 0.01, 0.02, 0.05)
    sampled_dir = '/share/home/deqing.gu/devdata/simulate_rawdata_for_ctDNA/downSample/downsampled_fastqs'
    pipeline_out = 'pipeline_result'
    main2(sampled_dir=sampled_dir, pipeline_out=pipeline_out,
         batch_id='LK292', panel='Oncoscreen.520.v3', sample_type='ctDNA',
         target_vcfs=target_vcfs, chunk_vcf_size=130, exclude_vcfs=exclude_vcfs,  # > split target vcf
         mu_type=mu_type, af_list=af_list,  pyflow=True
    )


if __name__ == '__main__':
    # ffpe simulate
    # test_sampling()
    # test_indel_main4()
    # test_snv_main3()

    # ctDNA simulate
    # test_sampling_ctdna()
    # test_indel_main_ctdna()
    # test_snv_main_ctdna

    # out_dir_set = None
    # time.sleep(3600*48)
    # after_simulate_analysis(out_dir_set=out_dir_set, batch_id='LK292', panel='Oncoscreen.520.v3', sample_type='ctDNA', output_sample_dir="/share/home/deqing.gu/devdata/pipeline_output_Sample")
    # after_simulate_analysis(out_dir_set=out_dir_set, batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE', output_sample_dir="/share/home/deqing.gu/devdata/pipeline_output_Sample")
    introduce_command(main)
