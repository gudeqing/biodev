import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor as Pool


def run_cmd(cmd):
    return subprocess.call(cmd, shell=True)

file = sys.argv[1]
targets = [x.strip() for x in open(file)]
cmds = []
for each in targets:
    cmd = 'wget {}'.format(each)
    cmds.append(cmd)

with Pool(5) as pool:
    pool.map(run_cmd, cmds)

