import re
import os

os.system('aws s3 ls --recursive s3://epionengs/80011001_HCC_Metastase/Methylation/ > all.files')
obj = open('all.files')
paths = [re.split('\s+', x, 3)[3] for x in open('all.files')]
obj.close()

with open('all.files.path', 'w') as f:
    _ = [f.write(x) for x in paths]
cmds = []
for path in paths:
    cmd = 'aws s3api restore-object '
    cmd += '--bucket epionengs '
    cmd += '--key "{}" '.format(path.strip())
    cmd += """--restore-request '{"Days":25,"GlacierJobParameters":{"Tier":"Standard"}}' """
    cmds.append(cmd)
    os.system(cmd)
