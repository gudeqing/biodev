import os
import subprocess


def find_md5_files(path):
    result = list()
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.endswith('.md5'):
                result.append(os.path.join(root, name))
                # print(result[-1])
    return result


def check_md5(path):
    for each in find_md5_files(path):
        print('start to check: {}'.format(each[:-4]))
        check_value = subprocess.check_output('md5sum {}'.format(each[:-4]), shell=True)
        check_value = check_value.split()[0].decode()
        known_value = open(each).read().strip().split()[0]
        if check_value != known_value:
            print('Failed to check md5 Value of {}'.format(each[:-4]))
            print('check_value: {}'.format(check_value))
            print('known_value: {}'.format(known_value))
        else:
            print('Success to check md5 Value of {}'.format(each[:-4]))


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['check_md5'])
