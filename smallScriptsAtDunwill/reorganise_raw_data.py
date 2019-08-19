#! /data/users/dqgu/anaconda3/bin/python
# coding=utf-8
import pandas as pd
import hashlib as hash
import easygui as eg
import os
import shutil
import sys
import glob
import gzip
import logging


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w+')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger

def query_type(query):
    query_dict = dict(
        TA01='AVENIO_Targeted_ctDNA',
        TA02='AVENIO_Expanded_ctDNA',
        TA03='AVENIO_Surveillance_ctDNA',
        TA11='AVENIO_Targeted_Tissue',
        TA12='AVENIO_Expanded_Tissue',
        TA13='AVENIO_Surveillance_Tissue',
        TU03='CHPv2',
        TU01='OCPv3',
        TB01='Lung_cfDNA',
        TB03='Lung_cfTNA',
        TB02='Breast_cfDNA_v1',
        TB04='Breast_cfDNA_v2',
        TB05='Colon_cfDNA',
        TB06='TCR_Beta_LR',
        TU05='TCR_Beta_SR_DNA',
        TU07='TCR_Beta_SR_RNA',
        TU04='Immune_Response',
        TE01='WES',
        TT01='WTS',
        TG03='WGS',
        Meth='TM01',
        TF01='Epione_800_Tissue',
        TF02='Epione_800_ctDNA',
        TC01='Epione_800_Cell',
        TC02='scWGS',
        TC03='scWES',
        TU06='Oncomine_TML',
        TB07='Oncomine_BRCA',
        TT02='CMS_RNA',
        TU08='Oncomine_Focus'
    )
    query = query.strip()
    if query in query_dict:
        return query_dict[query]

    else:
        return None


def check_integrity(file_dir, blocksize=2 ** 20):
    veri_failed = []
    veri_success = []
    hasher = hash.md5()
    for file in os.listdir(file_dir):
        if file.endswith(".gz"):
            with gzip.open(os.path.join(file_dir, file), mode='rb') as f:
                while True:
                    buf = f.read(blocksize)
                    if not buf:
                        break
                    hasher.update(buf)
                checksum = hasher.hexdigest()
            md5_file = os.path.join(file_dir, file + '.md5')
            with open(md5_file) as f:
                original_hashkey = f.readline().split()[0]
                # GUI tell if verification failed or successful
            if checksum != original_hashkey:
                veri_failed.append(file)

            else:
                veri_success.append(file)

    return veri_failed, veri_success


def copy_file_to_target_dir(source, target_part, target_full, veri_success, target_directory):
    file_exists = []
    file_des = []

    if os.path.exists(target_full):
        pass
    else:
        os.mkdir(target_part)
        os.mkdir(target_full)
    for item in os.listdir(source):
        item_path = os.path.join(source, item)
        # print(item_path, target_directory)
        if item_path == target_directory:
            if item.endswith('xlsx'):
                target_final = os.path.join(target_full, item)
                shutil.copy(item_path, target_final)
                continue
            if not os.path.isdir(item_path):
                continue
            for files in os.listdir(item_path):
                if files in veri_success:
                    file_path = os.path.join(source, item, files)
                    target_almost = os.path.join(target_full, item)
                    if not os.path.exists(target_almost):
                        os.mkdir(target_almost)
                    # print(target_almost, "THIS IS TARGET_ALMOST")
                    target_final = os.path.join(target_almost, files)

                    # Tell user that there is a duplicate directory (GUI)
                    if os.path.exists(target_final):
                        ph = "{}".format(target_final)
                        file_exists.append(ph)

                    else:
                        print('Copying  {} --> {}'.format(file_path, target_final))
                        shutil.copy2(file_path, target_final)
                        ph = "{}".format(target_final)
                        file_des.append(ph)
                        shutil.copy(file_path + '.md5', target_final + '.md5')

    return (file_exists, file_des)


def run(original_path, sample_info_excel, target_path, check=False):
    """
    重新组织原始数据
    :param original_path: 原下机数据路径
    :param sample_info_excel: 样本信息表
    :param target_path: 新的数据存放路径
    :param check: 是否检查数据完整性，默认不检查
    :return:
    """
    logger = set_logger(os.path.join(target_path, 'logger.txt'))
    missing_file = []
    veri_failed = []
    file_exists = []
    file_des = []
    veri_success = []
    integ = 0 if not check else 1

    info_pd = pd.read_excel(sample_info_excel).iloc[:, [1, 2, 3, 4, 5, 6, 7, 8]]
    for row in info_pd.iterrows():
        _, cols = row
        file_name_pattern = '*-{}'.format(cols[2].strip())
        file_name_pattern_path = os.path.join(original_path, file_name_pattern)
        find_result = glob.glob(file_name_pattern_path)
        if not find_result or len(find_result) > 1:
            # 如果用第四列信息找不到文件或者找到多个符合条件的文件，就改用第8列作为搜素条件
            file_name_pattern = '*-{}'.format(cols[6].strip().replace('_', '-'))
            file_name_pattern_path = os.path.join(original_path, file_name_pattern)
            find_result = glob.glob(file_name_pattern_path)

        if find_result:
            target_directory = find_result[0]
            x = cols[1].split('_')[1]
            x = query_type(x)
            target_part = os.path.join(target_path, str(cols[0]))
            target_full = os.path.join(target_path, str(cols[0]), x)

            if integ:
                check_failed, check_success = check_integrity(target_directory)
                veri_failed += check_failed
                veri_success += check_success
            else:
                for item in os.listdir(original_path):
                    item_path = os.path.join(original_path, item)
                    if not os.path.isdir(item_path):
                        continue
                    for files in os.listdir(item_path):
                        if files.endswith(".gz"):
                            veri_success.append(files)

            existed, copied = copy_file_to_target_dir(original_path, target_part, target_full, veri_success,
                                                            target_directory)
            file_exists += existed
            file_des += copied

        else:
            missing_file.append(cols[2])

    # GUI print failed verify
    if integ:
        logger.info("下列文件没有通过完整性检验, 已跳过这些文件:")
        msg = '\n'.join(veri_failed)
        logger.info(msg)
        title = "下列文件没有通过完整性检验，已跳过这些文件"
        eg.msgbox(msg, title=title)
    else:
        pass

    # GUI print missing files
    if missing_file:
        logger.info('找不到下列样本的文件:')
        msg = '\n'.join(missing_file)
        logger.info(msg)
        title = "找不到下列样本的文件，已跳过这些样本"
        eg.msgbox(msg, title=title)

    # GUI print already file at destination
    if file_exists:
        logger.info("下列文件已经在目的地存在, 请核对:")
        logger.info("\n".join(file_exists))
        title = "下列文件已经在目的地存在, 已跳过再次复制这些文件, 请核对!"
        msg = "\n".join(file_exists)
        eg.msgbox(msg, title=title)
    else:
        pass

    # GUI print final destination
    if file_des:
        logger.info("成功复制下列文件到目的地:")
        logger.info("\n".join(file_des))
        title = "成功复制下列文件到目的地"
        msg = "\n".join(file_des)
        eg.msgbox(msg, title=title)
    else:
        eg.msgbox("没有复制任何文件！！")
        logger.info("请知悉：没有复制任何文件！！！")


def gui():
    original_path = eg.diropenbox(u"选择下机数据目录")
    if not original_path:
        sys.exit(0)
    print(original_path)

    target_path = eg.diropenbox(u"选择拆分数据存放目录")
    if not target_path:
        sys.exit(0)
    print(target_path)

    excel = eg.fileopenbox(u"选择项目或样本信息表(*.xlsx)")
    if not excel:
        sys.exit(0)
    print(excel)

    msg = u"是否检查数据的完整性？"
    check = False
    if eg.ynbox(msg, choices=['是(较耗时)', '否']):
        check = True

    run(original_path, excel, target_path, check=check)


# gui()
if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['run'])
