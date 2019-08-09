# coding=utf-8
import pandas as pd
import hashlib as hash
import easygui as eg
import os
import shutil
import sys
import glob
import gzip


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
                        print('copying  {}'.format(file_path))
                        shutil.copy2(file_path, target_final)
                        ph = "{}".format(target_final)
                        file_des.append(ph)
                        shutil.copy(file_path + '.md5', target_final + '.md5')

    return (file_exists, file_des)


def run(original_path, sample_info_excel, target_path):
    missing_file = []
    veri_failed = []
    file_exists = []
    file_des = []
    veri_success = []

    msg = u"是否检查数据的完整性？"
    title = u"请确认是否检查数据的完整性"
    integ = 0
    if eg.ynbox(msg, title, ('Yes', 'No')):
        integ = 1
    else:
        pass

    info_pd = pd.read_excel(sample_info_excel).iloc[:, [1, 2, 3, 4, 5, 6, 7, 8]]
    for row in info_pd.iterrows():
        _, cols = row
        file_name_pattern = '*-{}'.format(cols[2].strip())
        # print(file_name_pattern, "THIS IS THE NAME PATTERN")
        file_name_pattern_path = os.path.join(original_path, file_name_pattern)
        find_result = glob.glob(file_name_pattern_path)
        if not find_result or len(find_result) > 1:
            file_name_pattern = '*-{}'.format(cols[6].strip().replace('_', '-'))
            file_name_pattern_path = os.path.join(original_path, file_name_pattern)
            find_result = glob.glob(file_name_pattern_path)
        target_directory = ""
        if find_result:
            target_directory = find_result[0]

            x = cols[1].split('_')[1]
            x = query_type(x)
            target_part = os.path.join(target_path, str(cols[0]))
            target_full = os.path.join(target_path, str(cols[0]), x)

            if integ:
                veri_failed, veri_success = check_integrity(target_directory)
            else:
                for item in os.listdir(original_path):
                    item_path = os.path.join(original_path, item)
                    if not os.path.isdir(item_path):
                        continue
                    for files in os.listdir(item_path):
                        if files.endswith(".gz"):
                            veri_success.append(files)

            file_exists, file_des = copy_file_to_target_dir(original_path, target_part, target_full, veri_success,
                                                            target_directory)

        else:
            missing_file.append(cols[2])

    # GUI print failed verify
    if integ:
        title = "Please Confirm Error(s)"
        msg = "Couldn't Verify these files."
        veri = eg.choicebox(msg, title, veri_failed)
        if not veri:
            sys.exit(0)
    else:
        pass

    # GUI print missing files
    if missing_file:
        title = "Please Confirm Error(s)"
        msg = "Couldn't find these files."
        veri = eg.choicebox(msg, title, missing_file)
        if not veri:
            sys.exit(0)

    else:
        pass

    # GUI print already file at destination
    if file_exists:
        title = "Please Confirm Error(s)"
        msg = "The following file(s) already exists at destination."
        veri = eg.choicebox(msg, title, file_exists)
        if not veri:
            sys.exit(0)
    else:
        pass

    # GUI print final destination
    if file_des:
        title = "Please Confirm Result"
        msg = "These files were copied"
        veri = eg.choicebox(msg, title, file_des)
        if not veri:
            sys.exit(0)
    else:
        pass


def gui():
    # Select original_path
    eg.msgbox(u"选择下机数据目录", title=u"选择下机数据目录")
    original_path = eg.diropenbox(u"选择下机数据目录", title=u"选择下机数据目录")
    if not original_path:
        sys.exit(0)
    print(original_path)
    # Select final_path
    eg.msgbox(u"选择拆分数据存放目录", title=u"选择拆分数据存放目录")
    target_path = eg.diropenbox(u"选择拆分数据存放目录", title=u"选择拆分数据存放目录")
    if not target_path:
        sys.exit(0)
    print(target_path)
    # Select input file
    eg.msgbox(u"选择拆分项目数据所需的样本项目信息文件(.xlsx)", title=u"选择拆分项目数据所需的样本项目信息文件(.xlsx)")
    excel = eg.fileopenbox(u"选择拆分项目数据所需的样本项目信息文件(.xlsx)", title=u"选择拆分项目数据所需的样本项目信息文件(.xlsx)")
    if not excel:
        sys.exit(0)
    print(excel)
    ##    title = "Files to be copied"
    ##    msg ="These are the files that will be copied"
    ##    choices = input.iloc[:,3].tolist()
    ##    eg.choicebox(msg, title, choices)
    run(original_path, excel, target_path)


gui()
