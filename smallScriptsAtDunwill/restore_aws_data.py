import re
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor as Pool


def run_cmd(cmd):
    return subprocess.call(cmd, shell=True)


def get_target_file_path(target_dir, out='target_file.list', full_match_exp='.*\.fastq.gz'):
    """
    获取
    :param target_dir: example "s3://epionengs/80011001_HCC_Metastase/Methylation/"
    :param out: 路径的输出文件
    :param full_match_exp: 使用re.fullmatch进行匹配的表达式, 只保留匹配上的路径
    :return: 路径列表
    """
    os.system(f'aws s3 ls --recursive {target_dir} > all.files')
    paths = [re.split('\s+', x, 3)[3] for x in open('all.files')]
    if full_match_exp:
        paths = [x for x in paths if re.fullmatch(full_match_exp, x.strip())]
    with open(out, 'w') as f:
        for each in paths:
                f.write(each)
    return paths


def tagging_files(path_file, bucket='epionengs', key='type', value='fastq'):
    """
    给目标文件打标签, 方便生命周期策略的使用
    :param path_file: 目标文件的路径信息文件, 通常是get_target_file_path的输出文件
    :param bucket: 存储桶的名称
    :param key: tag的键，先用这个键对文件进行分大类, 然后用tag的值分为小类, 使用生命周期策略时可以用'tag键值对'筛选出目标文件进行管理.
    :param value: tag的值
    :return: None
    """
    with open(path_file) as f:
        path_lst = [x.strip() for x in f]
    for each in path_lst:
        cmd = "aws s3api put-object-tagging --bucket {} ".format(bucket)
        cmd += '--key {} '.format(each)
        cmd += "--tagging 'TagSet={{Key={}, Value={}}}'".format(key, value)
        print(cmd)
        subprocess.check_call(cmd, shell=True)


def restore_files(path_file, bucket='epionengs'):
    """
    提取一份目标文件并使用标准存储, 且只保留2天
    :param path_file: 目标文件的路径信息文件, 通常是get_target_file_path的输出文件
    :param bucket: 存储桶的名称
    :return: None
    """
    with open(path_file) as f:
        path_lst = [x.strip() for x in f]
    for each in path_lst:
        cmd = 'aws s3api restore-object '
        cmd += '--bucket {} '.format(bucket)
        cmd += '--key "{}" '.format(each)
        cmd += """--restore-request '{"Days":2,"GlacierJobParameters":{"Tier":"Standard"}}' """
        print(cmd)
        subprocess.check_call(cmd, shell=True)


def restore_data(target_path, threads=0):
    """
    对指定的路径中的所有文件进行还原提取
    :param target_path: such as s3://epionengs/80011001_HCC_Metastase/Methylation/
    :param threads: thread number, 为0时不使用多线程，默认不使用，因为命令执行很快完成，aws会自动并发还原，无需等待
    :return:
    """

    print('还原该目录: ')
    os.system(f'aws s3 ls --recursive {target_path} > all.files')
    paths = [re.split('\s+', x, 3)[3] for x in open('all.files')]

    with open('all.files.path', 'w') as f:
        _ = [f.write(x) for x in paths]

    with open('restore.cmd.list', 'w') as f:
        cmds = []
        for path in paths:
            cmd = 'aws s3api restore-object '
            cmd += '--bucket epionengs '
            cmd += '--key "{}" '.format(path.strip())
            cmd += """--restore-request '{"Days":2,"GlacierJobParameters":{"Tier":"Standard"}}' """
            cmds.append(cmd)
            f.write(cmd + '\n')

    if threads:
        with Pool(threads) as pool:
            pool.map(run_cmd, cmds)
    else:
        for each in cmds:
            subprocess.check_call(each, shell=True)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())


