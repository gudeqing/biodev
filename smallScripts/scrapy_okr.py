import itertools
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as  EC
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException
from selenium .webdriver.support.select import Select
# from bs4 import BeautifulSoup
import sys
import time


# 设置下载选项
profile = webdriver.FirefoxProfile()
profile.set_preference('browser.download.dir', 'd:\\')
profile.set_preference('browser.download.folderList', 2)
profile.set_preference('browser.download.manager.showWhenStarting', False)
# 下面的类型根据https://www.w3school.com.cn/media/media_mimeref.asp查询
profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'application/pdf')
profile.set_preference("pdfjs.disabled", True) #这一行代码如果注释掉，会导致弹框依然出现
# 打开浏览器
options = Options()
options.add_argument('-headless')  # 无头参数
executable_path = 'D:\geckodriver-v0.26.0-win64\geckodriver.exe'
browser = webdriver.Firefox(executable_path=executable_path, firefox_profile=profile)
# wait = WebDriverWait(browser, 10)
# browser.implicitly_wait(10)

# 登录网站
login_site = "http://10.62.2.16:8088/#/login"
browser.get(login_site)
time.sleep(5)  # 等待加载完成，否则后续可能失败
browser.find_element_by_id('username').send_keys('gudeqing')
browser.find_element_by_id('password').send_keys('dunwill1220')
browser.find_element_by_id('submit').click()
time.sleep(2)


def get_gene_description():
    browser.get('http://10.62.2.16:8088/#/narrative')
    time.sleep(3)
    # key = browser.find_element_by_id('narrativeKey')
    # new_class_name = 'ui-select-container select2 select2-container ng-valid select2-container-active select2-dropdown-open open'
    # browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", key, 'class', new_class_name)
    # key.click()
    # time.sleep(1)
    # ul = key.find_element_by_class_name('select2-result-single')
    # lis = ul.find_elements_by_tag_name('li')
    # print('Available gene number is', len(lis))
    i = -1
    gene_num = 2
    desc_file = open('all.113.okr.Biomarker.description.txt', 'wb')
    while True:
        i += 1
        if i >= gene_num:
            break
        key = browser.find_element_by_id('narrativeKey')
        new_class_name = 'ui-select-container select2 select2-container ng-valid select2-container-active select2-dropdown-open open'
        browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", key, 'class', new_class_name)
        key.click()
        time.sleep(1)
        ul = key.find_element_by_class_name('select2-result-single')
        lis = ul.find_elements_by_tag_name('li')
        gene_num = len(lis)
        gene_obj = lis[i]
        gene_name = gene_obj.text.encode(encoding='UTF-8')
        gene_obj.click()
        time.sleep(1)
        # rows = key.find_elements_by_xpath('//div[@class="row ng-scope"]')
        rows = key.find_elements_by_xpath('//div[@ng-bind="section.text" and @class="ng-binding"]')
        print(gene_name)
        desc_file.write(gene_name+b'\n')
        desc_file.write(b'Background: '+ rows[0].text.encode(encoding='UTF-8') + b'\n')
        desc_file.write(b'Alterations and prevalence: '+ rows[1].text.encode(encoding='UTF-8') + b'\n')
        desc_file.write(b'Potential clinical relevance: '+ rows[2].text.encode(encoding='UTF-8') + b'\n')
        time.sleep(1)

    desc_file.close()


def get_report(mutations:list, cancer, filter_preset, template, report_name='report'):
    """
    :param mutations: 突变列表，如['ERBB3 mutation', 'TP53 mutation', 'Tumor Mutational Burden']
    :param cancer:
    :param filter_preset:
    :param template:
    :param report_name:
    :return:
    """
    # 切换到目标页面
    browsite = 'http://10.62.2.16:8088/#/browse'
    browser.get(browsite)
    time.sleep(3)

    # --------输入突变-------
    mut_txts = []
    for mut in mutations:
        print('selecting', mut)
        # 找到下面这个元素并点击，触发下拉框
        classes = browser.find_element_by_id('classes')
        # classes.click()
        # time.sleep(2)
        # 限定选择范围
        choices = classes.find_element_by_class_name('select2-search-field').find_element_by_tag_name('input')
        choices.send_keys(mut)
        time.sleep(1)
        # 触发下拉框后，可以定位所有下拉框的选项，他们存储在<ul>中
        ul = classes.find_element_by_class_name('select2-result-single')
        mutation_lst = ul.find_elements_by_tag_name('li')
        if len(mutation_lst) <= 0:
            raise Exception('可选突变列表为空了！')
        else:
            # print('Mutation number is', len(mutation_lst))
            pass
        for x in mutation_lst:
            if x.text == mut:
                x.click()
                mut_txts.append(mut)
                break
        time.sleep(1)

    # --------选择癌症类型-------
    # 找到下面这个元素并点击，触发下拉框
    indication = browser.find_element_by_id('indication')
    new_class_name = 'ui-select-container select2 select2-container ng-invalid ng-invalid-required ng-touched select2-container-active select2-dropdown-open open direction-up'
    browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", indication, 'class',
                           new_class_name)
    indication.click()
    time.sleep(1)

    # 触发下拉框后，可以定位所有下拉框的选项，他们存储在<ul>中
    class_name = '.ui-select-choices ui-select-choices-content select2-results ng-scope'.replace(' ', '.')
    # class_name = '.select2-result-single select2-result-sub'.replace(' ', '.') 这个名字已经失效，使用上面的
    uls = indication.find_elements_by_css_selector(class_name)
    cancer_lst = []
    for ul in uls:
        cancer_lst += ul.find_elements_by_tag_name('li')
    # print(cancer_lst)
    # print([x.text.encode() for x in cancer_lst])
    # with open('cancer_type.txt', 'wb') as f:
    #     for each in cancer_lst:
    #         f.write(each.text.encode(encoding='UTF-8') + b'\n')
    print('Cancer type number is', len(cancer_lst))
    # 在所有癌症类型中选择目标类型
    if type(cancer) == int:
        cancer = cancer_lst[cancer]
        cancer_txt = cancer.text
        cancer.click()
    else:
        for x in cancer_lst:
            if x.text == cancer:
                x.click()
                cancer_txt = cancer
                break
    time.sleep(2)

    # --------filterPreset-------
    # 找到下面这个元素并点击，触发下拉框
    preset = browser.find_element_by_id('filterPreset')
    preset.click()
    time.sleep(1)
    # 触发下拉框后，可以定位所有下拉框的选项，他们存储在<ul>中
    ul = preset.find_element_by_class_name('select2-result-single')
    filter_lst = ul.find_elements_by_tag_name('li')
    # print('Filter number', len(filter_lst))
    # print([x.text for x in filter_lst])

    if type(filter_preset) == int:
        filt = filter_lst[filter_preset]
        filt_txt = filt.text
        filt.click()
    else:
        for x in filter_lst:
            if x.text == filter_preset:
                x.click()
                filt_txt = filter_preset
                break
    time.sleep(2)

    # --------提交选项-------
    # 提交选项
    browser.find_element_by_css_selector('.btn.btn-primary').click()
    time.sleep(3)

    # ---进入报告页面并点击‘generate report’---
    browser.find_element_by_css_selector('.btn.btn-primary').click()
    time.sleep(4)

    # ---选择模板---
    temp = browser.find_element_by_id('reportTemplate')
    new_class_name = 'ui-select-container select2 select2-container ng-invalid ng-invalid-required ng-touched select2-container-active select2-dropdown-open open'
    # 修改temp的class才能触发下拉框
    browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", temp, 'class', new_class_name)
    temp.click()
    time.sleep(2)
    ul = temp.find_element_by_class_name('select2-result-single')
    tmp_list = ul.find_elements_by_tag_name('li')
    # print('report template number is', len(lis))
    # print([x.text for x in lis])
    if type(template) == int:
        tmp = tmp_list[template]
        tmp_txt = tmp.text
        tmp.click()
    else:
        for x in tmp_list:
            if x.text == template:
                tmp_txt = template
                x.click()
                break
    time.sleep(1)

    # summary
    print('query:', mut_txts)
    print('cancer:', cancer_txt)
    print('filterPreset:', filt_txt)
    print('reportTemplate:', tmp_txt)
    print('report_name:', report_name)

    # ---设置报告名称---
    report_name_obj = browser.find_element_by_name('reportFilename')
    report_name_obj.clear()
    report_name_obj.send_keys(report_name)
    time.sleep(1)

    # --------- download pdf -----------
    txt_button = browser.find_element_by_xpath("//input[@value='{}']".format('PDF'))
    new_class_name = 'ng-valid ng-valid-required ng-dirty ng-touched ng-valid-parse'
    # 修改txt_button才能勾选Txt选项
    browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", txt_button, 'class',
                           new_class_name)
    txt_button.click()
    time.sleep(1)

    # download report
    download_button = browser.find_element_by_css_selector('.btn.btn-primary')
    download_button.click()
    time.sleep(5)

    # # 模拟鼠标点击，由于已可以自动下载，无需
    # import mouse
    # # print(mouse.get_position())
    # time.sleep(2)
    # mouse.move(718, 459)
    # mouse.click('left')
    # time.sleep(2)

    # ----------download TXT---------
    txt_button = browser.find_element_by_xpath("//input[@value='{}']".format('TXT'))
    new_class_name = 'ng-valid ng-valid-required ng-dirty ng-touched ng-valid-parse'
    # 修改txt_button才能勾选Txt选项
    browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", txt_button, 'class',
                           new_class_name)
    txt_button.click()
    time.sleep(1)

    # download report
    download_button = browser.find_element_by_css_selector('.btn.btn-primary')
    download_button.click()
    time.sleep(3)
    #
    # import mouse
    # # print(mouse.get_position())
    # time.sleep(2)
    # mouse.move(718, 459)
    # mouse.click('left')
    # time.sleep(1)


if __name__ == '__main__':
    # mutations = ['FBXW7 deleterious mutation', 'KRAS G12 mutation', 'ERBB3 mutation', 'TP53 mutation']
    # mutations = ['KRAS G12 mutation', 'ERBB3 mutation', 'TP53 mutation']
    mutations = ['KRAS G12 mutation', 'ERBB3 mutation']
    # mutations.append('Tumor Mutational Burden'),
    # mutations += ['Microsatellite stable']
    # mutations += ['Microsatellite instability-High']
    # mutations += ['Microsatellite instability-Low']
    get_report(
        mutations,
        # cancer='Non-Small Cell Lung Cancer',
        cancer='Rectal Cancer',
        filter_preset='Epi800_Test',
        template='CL_800_Panel',
    )

# 串行步骤：
# （1）借助python的paramiko包把命令投递到服务器完成annovar等注释
# （2）window OKR爬虫注释，借助paramiko把结果上传到服务器
# （3）前面两步准备好输入文件后，完成报告生成

# if __name__ == '__main__':
#     from xcmds import xcmds
#     xcmds.xcmds(locals(), include=['get_report'])
