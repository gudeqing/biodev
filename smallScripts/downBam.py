import pandas as pd
import os
from concurrent.futures import ThreadPoolExecutor
import subprocess


def run_cmd(cmd):
    return subprocess.call(cmd, shell=True)


def downBam2Fq(data):
    """
    samtools sort + view + fastq -> down sample bam to fastq
    :param data: 第一列样本名，第二列为bam路径，其他列为百分比，第二列的header需指定为"path"
    :return:
    """
    data = pd.read_csv(data, header=0, index_col=0, sep=None, engine='python')
    cmd_list = list()
    for sample in data.index:
        bam = data.loc[sample, 'path']
        for name in data.columns[2:]:
            ratio = data.loc[sample, name]
            if not os.path.exists(name):
                os.mkdir(name)
            fq = f'{name}/{sample}.{name}.R1.fq.gz'
            fq2 = f'{name}/{sample}.{name}.R2.fq.gz'
            cmd = f'samtools sort -n {bam}| samtools view -s {ratio} -O BAM - | samtools fastq -1 {fq} -2 {fq2} - '
            cmd_list.append(cmd)
    with ThreadPoolExecutor(6) as pool:
        pool.map(run_cmd, cmd_list)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['downBam2Fq'])
