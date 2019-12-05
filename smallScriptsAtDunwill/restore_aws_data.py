import re
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor as Pool


def run_cmd(cmd):
    return subprocess.call(cmd, shell=True)


def restore_data(target_path, threads=0):
    """
    还原aws上深度存储的数据，方便下载
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
    xcmds.xcmds(locals(), include=['restore_data'])


