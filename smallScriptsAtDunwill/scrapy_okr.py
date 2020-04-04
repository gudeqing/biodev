import itertools
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
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
profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/plain, application/pdf, application/octet-stream, text/plain;charset=UTF-8, text/css')

# 打开浏览器
options = Options()
options.add_argument('-headless')  # 无头参数
executable_path = 'D:\\firefoxdriver\\geckodriver.exe'
browser = webdriver.Firefox(executable_path=executable_path, firefox_profile=profile)
# wait = WebDriverWait(browser, 10)
# browser.implicitly_wait(10)

# 登录网站
login_site = "http://10.62.2.16:8088/#/login"
browser.get(login_site)
time.sleep(5)  # 等待加载完成，否则后续可能失败
browser.find_element_by_id('username').send_keys('dqgu')
browser.find_element_by_id('password').send_keys('dunwill1220')
browser.find_element_by_id('submit').click()
time.sleep(2)


def get_report():
    mutation_index = range(1058)
    cancer_index = range(67)
    filter_index = range(6)
    number = 0 
    for m_idx, c_idx, f_idx in itertools.product(mutation_index, cancer_index, filter_index):
        number += 1
        if number > 3:
            break

        # 切换到目标页面
        browsite = 'http://10.62.2.16:8088/#/browse'
        browser.get(browsite)
        time.sleep(3)

        for mu_id in [m_idx, m_idx+3]:
            # 找到下面这个元素并点击，触发下拉框
            classes = browser.find_element_by_id('classes')
            classes.click()
            time.sleep(2)
            # 触发下拉框后，可以定位所有下拉框的选项，他们存储在<ul>中
            ul = classes.find_element_by_class_name('select2-result-single')
            mutation_lst = ul.find_elements_by_tag_name('li')
            print('Mutation number is', len(mutation_lst))
            mutation = mutation_lst[m_idx]
            mutation_txt = mutation.text
            mutation.click()
            time.sleep(1)

        # 找到下面这个元素并点击，触发下拉框
        indication = browser.find_element_by_id('indication')
        new_class_name = 'ui-select-container select2 select2-container ng-invalid ng-invalid-required ng-touched select2-container-active select2-dropdown-open open direction-up'
        browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", indication, 'class', new_class_name)
        indication.click()
        time.sleep(1)
        
        # 触发下拉框后，可以定位所有下拉框的选项，他们存储在<ul>中
        class_name = '.select2-result-single select2-result-sub'.replace(' ', '.')
        uls = indication.find_elements_by_css_selector(class_name)
        cancer_lst = []
        for ul in uls:
            cancer_lst += ul.find_elements_by_tag_name('li')
        print([x.text.encode() for x in cancer_lst])
        with open('cancer_type.txt', 'wb') as f:
            for each in cancer_lst:
                f.write(each.text.encode(encoding='UTF-8')+b'\n')

        print('Cancer type number is', len(cancer_lst))
        cancer = cancer_lst[c_idx]
        cancer_txt = cancer.text
        cancer.click()
        time.sleep(2)

        # 找到下面这个元素并点击，触发下拉框
        preset = browser.find_element_by_id('filterPreset')
        preset.click()
        time.sleep(1)
        # 触发下拉框后，可以定位所有下拉框的选项，他们存储在<ul>中
        ul = preset.find_element_by_class_name('select2-result-single')
        filter_lst = ul.find_elements_by_tag_name('li')
        print('Filter number', len(filter_lst))
        print([x.text for x in filter_lst])
        filt = filter_lst[f_idx]
        filt_txt = filt.text
        filt.click()
        time.sleep(2)

        # 提交选项
        browser.find_element_by_css_selector('.btn.btn-primary').click()
        time.sleep(3)

        # 进入报告页面并点击‘generate report’
        browser.find_element_by_css_selector('.btn.btn-primary').click()
        time.sleep(4)

        # 选择模板
        # select report template
        temp = browser.find_element_by_id('reportTemplate')
        new_class_name = 'ui-select-container select2 select2-container ng-invalid ng-invalid-required ng-touched select2-container-active select2-dropdown-open open'
        # 修改temp的class才能触发下拉框
        browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", temp, 'class', new_class_name)
        temp.click()
        time.sleep(2)
        ul = temp.find_element_by_class_name('select2-result-single')
        lis = ul.find_elements_by_tag_name('li')
        print('report template number is', len(lis))
        print([x.text for x in lis])
        lis[2].click()
        time.sleep(1)

        # 设置报告名称
        my_report_name = '('+mutation_txt+')'+'('+cancer_txt+')'+'('+filt_txt+')'
        report_name_obj = browser.find_element_by_name('reportFilename')
        report_name_obj.clear()
        report_name_obj.send_keys(my_report_name)
        time.sleep(1)

        # format selection
        txt_button = browser.find_element_by_xpath("//input[@value='PDF']")
        new_class_name = 'ng-valid ng-valid-required ng-dirty ng-touched ng-valid-parse'
        # 修改txt_button才能勾选Txt选项
        browser.execute_script("arguments[0].setAttribute(arguments[1],arguments[2])", txt_button, 'class', new_class_name)
        txt_button.click()
        time.sleep(1)

        # download report
        browser.find_element_by_css_selector('.btn.btn-primary').click()
        time.sleep(3)
        # 切换到目标页面
        browsite = 'http://10.62.2.16:8088/#/browse'
        browser.get(browsite)
        time.sleep(3)

# browser.quit()

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
        

get_gene_description()


with open('known.mutations.xls') as f, open('final.hotspots.vcf') as f2:
    query = dict()
    for line in f2:
        hgvs = line.strip().split('=', 1)[1]
        key = hgvs.split(':')[0] + ':' + hgvs.split(':')[3]
        query[key] = hgvs
    with open('new.known.mutations.xls', 'w') as f3:
        f3.write(f.readline())
        for line in f:
            sample, hgvs, af = line.strip().split()
            key = hgvs.split(':')[0] + ':' + hgvs.split(':')[3]
            if key in query:
                f3.write(f'{sample}\t{query[key]}\t{af}\n')
            else:
                print('failed:', line)
