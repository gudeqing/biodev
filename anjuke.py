import os
import itertools
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.common.exceptions import NoSuchElementException
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

import json
import os
import random
import re
import sys
import time
import traceback
from datetime import datetime
from typing import List, Union, Iterable
from urllib.parse import urljoin

from lxml import etree
from lxml.etree import Element, _Element
from lxml.html import tostring

from selenium import webdriver
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait

def random_sleep():
    time.sleep(random.random() * (random.random() * 2))


def random_wait():
    time.sleep(random.random() * (random.random() * 2) + random.randint(2, 4))

# 阻塞【等待所有数据加载完成】
TIME_OUT = 15
TOTAL_PAGE = 10
function_type = type(lambda: None)

def wait_wrapper(func: function_type, browser):
    wait = WebDriverWait(browser, TIME_OUT)
    is_success = False
    for i in range(3):
        try:
            wait.until(func)
            is_success = True
            break
        except TimeoutException:
            browser.refresh()
            random_wait()
    if not is_success:
        raise Exception(f'点击跳转失败！{browser.current_url}')


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

# 进入安居客城市选项页面, 获得所有城市的网址列表
city_info_file = 'anjuke_city_web_sites.txt'
if not os.path.exists(city_info_file):
    browser.get('https://www.anjuke.com/sy-city.html')
    time.sleep(2)
    city_web_sites = dict()

    choices = browser.find_elements_by_class_name('city_list')
    for choice in choices:
        for each in choice.find_elements_by_tag_name('a'):
            city = each.text
            city_web_site = each.get_attribute('href')
            city_web_sites[city] = city_web_site
    print(f'there are {len(city_web_sites)} cities!')
    with open(city_info_file, 'w') as f:
        for k, v in city_web_sites.items():
            f.write(f'{k}\t{v}\n')
else:
    city_web_sites = dict(x.strip().split() for x in open(city_info_file))

# 进入某个城市的页面, 获得租、售网址
target_city = list(city_web_sites.keys())[0]
city_site = city_web_sites[target_city]

browser.get(city_site)
zu_shou_sites = set()
for choice in browser.find_elements_by_class_name('third_navlist'):
    for each in choice.find_elements_by_tag_name('a'):
        sell_type_site = each.get_attribute('href')
        if sell_type_site.endswith('/sp-zu/'):
            zu_site = sell_type_site
            zu_shou_sites.add(zu_site)
        elif sell_type_site.endswith('/sp-shou/'):
            sh_site = sell_type_site
            zu_shou_sites.add(sh_site)
        if len(zu_shou_sites) == 2:
            break

zu_shou_sites = sorted(zu_shou_sites)

# 进入商铺出售的页面
browser.get(zu_shou_sites[0])
time.sleep(2)
while True:
    for choice in browser.find_elements_by_class_name('list-item'):
        # 提取detail_url，如果想进入具体页面抓取信息，则可能需要进行验证
        detail_url = choice.find_element_by_tag_name('a').get_attribute('href')
        # 下面直接在当前页面提取信息，唯一不能提取的信息是经纬度地址
        # 提取title，地址信息
        item_info = choice.find_element_by_class_name('item-info')
        title = item_info.find_element_by_class_name('item-title').find_element_by_tag_name('span').text
        print(title)
        info_lst = []
        for desc in item_info.find_elements_by_class_name('item-descript'):
            text = ''
            for each in desc.find_elements_by_tag_name('span'):
                text += each.text
            info_lst.append(text)
        address, floor_info, source_info = info_lst
        floor_info = floor_info.split()[0]
        print('Record', title, info_lst)

        # 提取面积信息
        area_info = choice.find_element_by_class_name('item-area')
        area_size = area_info.find_element_by_class_name('area').find_element_by_tag_name('span').text

        # 提取单价信息
        try:
            price_info = choice.find_element_by_class_name('item-price').find_element_by_class_name('price-daily')
            price_detail = price_info.find_elements_by_tag_name('span')
            _, price, unit = [x.text for x in price_detail]
            if '万' in unit:
                price = str(float(price)*10000)
                unit = "元/m²"
        except NoSuchElementException as e:
            price = 'unknown'
            unit = 'unknown'

        print(area_size, price, unit)

        # 点击进入详情页获取信息
        detail_url.click()
        random_wait()
        random_sleep()
        # 获取经纬度



    # 点击下一页
    try:
        next_page = browser.find_element_by_class_name('aNxt')
        next_page.click()
        random_sleep()
        random_wait()
        wait_wrapper(EC.visibility_of_element_located((By.CLASS_NAME, 'aNxt')), browser)
    except Exception as e:
        print(e)
        break






# 在上面的页面找到“https://sh.sydc.anjuke.com/sp-zu/”

# 进入出租页面



