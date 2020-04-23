import re
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor as Pool

_help = '''
从aws下载深度存储的原始数据, 还原数据需要收费，还原后的数据也要收费，谨慎还原！

1. 从aws下载数据，如果数据已经深度归档，首先要还原，然后才能下载。
    我这边提供了脚本辅助批量还原：
    （1）首先使用get_target_file_path获得某个目录下所有文件路径信息，提供正则匹配参数筛选目标文件。
    （2）使用restore_files 对目标文件进行还原，还原通常需要12小时。

2. 辅助脚本路径：
    /data/users/dqgu/PycharmProjects/biodev/smallScriptsAtDunwill/restore_aws_data.py 

3. 该脚本包含若干字命令，列出子命令：
```
$ python restore_aws_data.py 
The tool has the following sub-commands: 
Pool
run_cmd
get_target_file_path（*本次用到）
tagging_files（对目标文件打标签，便于进行生命周期策略管理）
restore_files（*本次用到，基于提供的目标文件路径逐一还原）
restore_data （对目标目录的所有文件进行还原）
```

4. 获得子命令使用说明示例：
```
$ python restore_aws_data.py get_target_file_path -h
usage: get_target_file_path [-h] -target_dir target_dir
                            [-out Default:target_file.list]
                            [-full_match_exp Default:.*\.fastq.gz]

optional arguments:
  -h, --help            show this help message and exit
  -target_dir target_dir
                         example "s3://epionengs/80011001_HCC_Metastase/Methylation/"
  -out Default:target_file.list
                         路径的输出文件
  -full_match_exp Default:.*\.fastq.gz
                         使用re.fullmatch进行匹配的表达式, 只保留匹配上的路径
```

5. 还原举例：
（1）python restore_aws_data.py get_target_file_path -t s3://epionengs/80011001_HCC_Metastase/WES -full_match_exp '.*\.fastq.gz' -out target_file.list
（2）python restore_aws_data.py restore_files -path  target_file.list

6. 下载已经还原的数据，注意，还原数据是比较慢的，通常隔天才能完成，除非你愿意加钱。
（1）举例：aws s3 sync --only-show-errors s3://epionengs/80028_sc_RNA local_path > download.log &  
 (2) 如果有多个文件需要下载，可以写shell循环完成，暂未提供批量下载的脚本

'''


def readme():
    print(_help)


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
    return [x.strip() for x in paths]


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
    if type(path_file) == str:
        with open(path_file) as f:
            path_lst = [x.strip() for x in f]
    else:
        assert type(path_file) == list
        path_lst = path_file
    for each in path_lst:
        cmd = 'aws s3api restore-object '
        cmd += '--bucket {} '.format(bucket)
        cmd += '--key "{}" '.format(each)
        cmd += """--restore-request '{"Days":2,"GlacierJobParameters":{"Tier":"Standard"}}' """
        print(cmd)
        subprocess.check_call(cmd, shell=True)
        print('please wait about 12 hours and then download your target file by using other command')
        print('such as: aws s3 sync s3://epionengs/80011024_CRISPR/190/ 190/')
        print('Restored data will be deleted after two days!')



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


def download(files, outdir=os.getcwd(), bucket='epionengs', threads=3):
    cmds = list()
    for each in files:
        cmd = 'aws s3 cp s3://{bucket}/{each} {outdir}'.format(
            bucket=bucket, each=each, outdir=outdir
        )
        cmds.append(cmd)
    with Pool(threads) as pool:
        pool.map(run_cmd, cmds)


def find_restore_download(target_dir, outdir, match='.*\.fastq.gz', restore=False, threads=3, do_not_download=False):
    """

    :param target_dir: example "s3://epionengs/80011001_HCC_Metastase/Methylation/"
    :param outdir: 下载结果的输出路径
    :param match: 匹配文件的表达式
    :param restore: 如果提供该参数，则下载前需要对文件进行还原
    :param threads: 下载的并发数
    :param do_not_download: 如果提供该参数，则不下载
    :return:
    """
    out = os.path.join(outdir, 'target_file.list')
    bucket = target_dir.split('/')[2]
    target_files = get_target_file_path(target_dir, out, match)
    if restore:
        if target_files:
            restore_files(target_files, bucket=bucket)
        else:
            print('Nothing matched!')
            return
        import time

    if not do_not_download:
        if restore:
            import time
            time.sleep(3600*12)
        download(target_files, outdir, bucket, threads)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())


