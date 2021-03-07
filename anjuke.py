# coding=utf-8
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
from pprint import pprint

# 获取IP：http://h.zhimaruanjian.com/getapi/

target_info_dict = dict(
    province = '省/直辖市',
    city = '市',
    district = '区',
    url = '详情页URL',
    type = '商铺类型[出租/出售]',
    title = '发布信息Title',
    price = '价格[单位：元]',
    unit = '价格单位[元 / 月 ， 出售]',
    payment_method = '支付方式[例如：押一付三]',
    area = '面积[平方米]',
    width = '宽度[米]',
    height = '高度[米]',
    level = '楼层[高，中，低]',
    houses_name = '楼盘名称',
    houses_address = '楼盘地址',
    pos_type = '经纬度获取方式[raw: 爬取，translate：百度转换]',
    longitude = '经度[百度]',
    latitude = '纬度[百度]',
    floor_number = '楼层数',
    published_datetime = '发布日期',
)


def set_chrome(proxy):
    options = webdriver.ChromeOptions()
    options.add_argument("ignore-certificate-errors")
    options.add_argument("--ignore-ssl-errors")
    # 设置代理
    if proxy:
        ip, port = proxy['ip'], proxy['port']
        print(f'***  使用的IP：{ip}, PORT：{port}')
        options.add_argument(f"--proxy-server={ip}:{port}")
    else:
        print('未使用代理！')
    options.add_argument('lang=zh-CN,zh,zh-TW,en-US,en')
    options.add_argument('disable-infobars')
    options.add_experimental_option('useAutomationExtension', False)
    options.add_experimental_option('excludeSwitches', ['enable-automation'])
    preferences = {
        "profile.managed_default_content_settings.images": 2,
        "webrtc.ip_handling_policy": "disable_non_proxied_udp",
        "webrtc.multiple_routes_enabled": False,
        "webrtc.nonproxied_udp_enabled": False,
    }
    options.add_experimental_option("prefs", preferences)
    # 无头模式，不开启浏览器
    # options.add_argument('--headless')

    browser = webdriver.Chrome(options=options, executable_path='D:\Drivers\chromedriver.exe')
    # 反屏蔽[检测是否selenium自动爬取，其大多数情况下，检测基本原理是检测当前浏览器窗口下的 window.navigator 对象是否包含 webdriver 这个属性
    # ，这个属性是 undefined来解决]
    browser.execute_cdp_cmd(
        'Page.addScriptToEvaluateOnNewDocument',
        {'source': 'Object.defineProperty(navigator, "webdriver", {get: () => undefined})'}
    )
    return browser


def random_sleep():
    time.sleep(random.random() * (random.random() * 2))


def random_wait():
    time.sleep(random.random() * (random.random() * 2) + random.randint(2, 4))


def wait_wrapper(func, browser):
    # 阻塞【等待所有数据加载完成】
    TIME_OUT = 15
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


def set_firefox(proxy:dict=None):
    profile = webdriver.FirefoxProfile()
    profile.set_preference('browser.download.dir', 'd:\\')
    profile.set_preference('browser.download.folderList', 2)
    profile.set_preference('browser.download.manager.showWhenStarting', False)
    # 下面的类型根据https://www.w3school.com.cn/media/media_mimeref.asp查询
    profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'application/pdf')
    profile.set_preference("pdfjs.disabled", True)
    # 这一行代码如果注释掉，会导致弹框依然出现
    # 打开浏览器
    options = Options()
    options.add_argument('-headless')  # 无头参数
    options.add_argument("ignore-certificate-errors")
    options.add_argument("--ignore-ssl-errors")
    # 设置代理
    if proxy:
        ip, port = proxy['ip'], proxy['port']
        print(f'***  使用的IP：{ip}, PORT：{port}')
        options.add_argument(f"--proxy-server={ip}:{port}")
    else:
        print('未使用代理！')
    executable_path = 'D:\Drivers\geckodriver.exe'
    browser = webdriver.Firefox(executable_path=executable_path, firefox_profile=profile)
    # wait = WebDriverWait(browser, 10)
    # browser.implicitly_wait(10)
    return browser


def get_city_urls(start_url='https://www.anjuke.com/sy-city.html'):
    browser.get(start_url)
    time.sleep(2)
    wait_wrapper(EC.visibility_of_element_located((By.CLASS_NAME, 'city_list')), browser)
    city_web_sites = dict()
    choices = browser.find_elements_by_class_name('city_list')
    for choice in choices:
        for each in choice.find_elements_by_tag_name('a'):
            city = each.text
            city_web_site = each.get_attribute('href')
            city_web_sites[city] = city_web_site
    print(f'there are {len(city_web_sites)} cities!')
    return city_web_sites


def get_zu_shou_url(city_url):
    browser.get(city_url)
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
    return sorted(zu_shou_sites)


def get_detail_url(start_url):
    browser.get(start_url)
    time.sleep(2)
    urls = []
    while True:
        # 存在页面失效问题，需要先把详情页面的网址先爬取完成，再逐一解析
        for choice in browser.find_elements_by_class_name('list-item'):
            # 提取detail_url，如果想进入具体页面抓取信息，则可能需要进行验证
            detail_url = choice.find_element_by_tag_name('a').get_attribute('href')
            urls.append(detail_url)
        # 点击下一页
        try:
            next_page = browser.find_element_by_class_name('aNxt')
            next_page.click()
            random_sleep()
            wait_wrapper(EC.visibility_of_element_located((By.CLASS_NAME, 'aNxt')), browser)
        except Exception as e:
            print(e)
            break
    print(f'find {len(urls)} detail urls')
    return urls


def get_detail_info(url, city):
    browser.get(url)
    random_wait()
    # 提取title
    try:
        title = browser.find_element_by_tag_name('head').find_element_by_xpath('//meta[@name="keywords"]')
        title = title.get_attribute('content')
    except NoSuchElementException:
        print('可能需要验证了, 请三十秒的内通过验证！')
        time.sleep(30)
        try:
            browser.get(url)
            random_wait()
            title = browser.find_element_by_tag_name('head').find_element_by_xpath('//meta[@name="keywords"]')
            title = title.get_attribute('content')
        except Exception as e:
            print(e)
            return 'change_IP'
    except Exception as e:
        print(e)
        return 'change_IP'

    # 获取更新日期
    try:
        pup_date = browser.find_element_by_class_name('site-item-date').text
        # pup_date = re.search('(\d+月\d+)日', pup_date).groups()[0]
    except NoSuchElementException as e:
        print(e)
        pup_date = 'unknown'
    # 提取基础信息
    basic_info = browser.find_element_by_class_name('basic-info-wrapper').find_elements_by_class_name('item')
    basic_info_dict = dict()
    for each in basic_info:
        tmp = []
        for i in each.find_elements_by_tag_name('span'):
            tmp.append(i.text)
        if not tmp[1].strip():
            tmp[1] = 'unknown'
        basic_info_dict[tmp[0]] = tmp[1]

    # 获取经纬度 baiduMap['lat'] 和 baidMap['lng']
    page = browser.page_source
    try:
        lat = re.search('lat: (.*),?', page).groups()[0]
        lng = re.search('lng: (.*),?', page).groups()[0]
    except Exception as e:
        print(e)
        lat = 'unknown'
        lng = 'unknown'
    basic_info_dict['latitude'] = lat
    basic_info_dict['longitude'] = lng
    basic_info_dict['published_datetime'] = pup_date
    # 加工信息
    district = basic_info_dict['地址'].split()[0]
    basic_info_dict['district'] = district
    if '总价' in basic_info_dict:
        basic_info_dict['type'] = 1
        if "面" in basic_info_dict['总价']:
            price = '面议'
        else:
            if '万' in basic_info_dict['总价']:
                price = float(basic_info_dict['总价'].split('万')[0])*10000/float(basic_info_dict['建筑面积'].split('m')[0])
            elif '亿' in basic_info_dict['总价']:
                price = float(basic_info_dict['总价'].split('亿')[0]) * 10**9 / float(
                    basic_info_dict['建筑面积'].split('m')[0])
            else:
                price = float(basic_info_dict['总价'].split('元')[0])*10000/float(basic_info_dict['建筑面积'].split('m')[0])
        unit = '元/m2'
    else:
        basic_info_dict['type'] = 0
        yuezu = basic_info_dict['月租'].split('/')
        unit = '元/月/m2'
        if '万' in yuezu:
            price = float(yuezu.split('万')[0])*10000/float(basic_info_dict['建筑面积'].split('m')[0])
        elif '亿' in yuezu:
            price = float(yuezu.split('亿')[0])*10**9/float(basic_info_dict['建筑面积'].split('m')[0])
        elif '/' in yuezu:
            price = float(yuezu.split('/')[0])/float(basic_info_dict['建筑面积'].split('m')[0])
        else:
            price = '面议'
    basic_info_dict['price'] = price
    basic_info_dict['unit'] = unit
    width, height, level = re.findall('(\d+)m', basic_info_dict['规格'])
    basic_info_dict['width'] = width
    basic_info_dict['height'] = height
    basic_info_dict['level'] = level
    basic_info_dict['url'] = url
    basic_info_dict['title'] = title
    basic_info_dict['city'] = city
    basic_info_dict['pos_type'] = 'raw'
    # rename
    basic_info_dict['payment_method'] = basic_info_dict.get('押付', 'unknown')
    basic_info_dict['area'] = basic_info_dict['建筑面积']
    basic_info_dict['houses_name'] = basic_info_dict['物业']
    basic_info_dict['houses_address'] = basic_info_dict['地址']
    basic_info_dict['floor_number'] = basic_info_dict['楼层']
    return basic_info_dict


# 获取出售信息
def pipeline():
    proxy_lst = [
        {"ip": "122.232.230.213", "port": 4278, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "223.215.177.242", "port": 4245, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "36.62.210.254", "port": 4216, "expire_time": "2021-03-07 19:12:11"},
         {"ip": "42.56.239.116", "port": 4278, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "27.190.81.193", "port": 4282, "expire_time": "2021-03-07 20:51:56"},
         {"ip": "49.85.111.120", "port": 4257, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "175.154.202.156", "port": 4258, "expire_time": "2021-03-07 20:52:40"},
         {"ip": "117.70.41.177", "port": 4245, "expire_time": "2021-03-07 21:31:28"},
         {"ip": "27.152.193.52", "port": 4273, "expire_time": "2021-03-07 19:07:20"},
         {"ip": "59.54.125.30", "port": 4275, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "153.99.10.54", "port": 4207, "expire_time": "2021-03-07 20:39:15"},
         {"ip": "106.57.168.90", "port": 4282, "expire_time": "2021-03-07 19:36:00"},
         {"ip": "49.67.92.148", "port": 4264, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "220.161.33.254", "port": 4245, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "117.28.60.189", "port": 4252, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "114.104.139.166", "port": 4214, "expire_time": "2021-03-07 19:34:00"},
         {"ip": "114.106.170.127", "port": 4245, "expire_time": "2021-03-07 20:08:01"},
         {"ip": "221.234.31.198", "port": 4245, "expire_time": "2021-03-07 19:30:22"},
         {"ip": "125.87.86.51", "port": 4278, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "117.64.254.209", "port": 4251, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "27.44.218.67", "port": 4245, "expire_time": "2021-03-07 20:09:27"},
         {"ip": "106.57.168.239", "port": 4280, "expire_time": "2021-03-07 20:56:03"},
         {"ip": "116.54.250.216", "port": 4283, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "113.141.223.37", "port": 4236, "expire_time": "2021-03-07 20:33:40"},
         {"ip": "27.190.82.153", "port": 4278, "expire_time": "2021-03-07 20:28:20"},
         {"ip": "125.78.218.73", "port": 4237, "expire_time": "2021-03-07 20:35:23"},
         {"ip": "60.184.199.115", "port": 4223, "expire_time": "2021-03-07 21:00:18"},
         {"ip": "183.7.113.72", "port": 4230, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "60.166.181.141", "port": 4231, "expire_time": "2021-03-07 21:08:54"},
         {"ip": "125.78.217.236", "port": 4237, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "60.185.32.204", "port": 4234, "expire_time": "2021-03-07 19:04:02"},
         {"ip": "183.7.140.235", "port": 4230, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "220.164.105.105", "port": 4235, "expire_time": "2021-03-07 21:38:56"},
         {"ip": "122.136.164.246", "port": 4278, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "113.100.9.62", "port": 4245, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "125.87.87.5", "port": 4278, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "119.54.46.145", "port": 4212, "expire_time": "2021-03-07 21:01:57"},
         {"ip": "106.40.145.249", "port": 4283, "expire_time": "2021-03-07 19:47:27"},
         {"ip": "110.90.220.187", "port": 4210, "expire_time": "2021-03-07 19:52:46"},
         {"ip": "123.181.148.224", "port": 4278, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "222.219.252.56", "port": 4256, "expire_time": "2021-03-07 20:45:39"},
         {"ip": "110.90.222.128", "port": 4260, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "36.57.68.99", "port": 4227, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "182.240.228.141", "port": 4245, "expire_time": "2021-03-07 19:09:08"},
         {"ip": "49.65.165.20", "port": 4232, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "111.227.42.3", "port": 4278, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "180.119.19.209", "port": 4245, "expire_time": "2021-03-07 19:03:22"},
         {"ip": "36.6.146.194", "port": 4225, "expire_time": "2021-03-07 19:13:53"},
         {"ip": "58.21.243.144", "port": 4278, "expire_time": "2021-03-07 21:20:41"},
         {"ip": "60.161.152.114", "port": 4246, "expire_time": "2021-03-07 20:39:24"},
         {"ip": "220.201.85.200", "port": 4260, "expire_time": "2021-03-07 20:54:14"},
         {"ip": "60.173.34.155", "port": 4276, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "112.123.40.112", "port": 4254, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "116.248.172.91", "port": 4221, "expire_time": "2021-03-07 20:58:43"},
         {"ip": "36.102.175.175", "port": 4285, "expire_time": "2021-03-07 20:03:50"},
         {"ip": "220.164.105.63", "port": 4235, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "36.102.169.219", "port": 4235, "expire_time": "2021-03-07 19:12:13"},
         {"ip": "49.85.85.2", "port": 4257, "expire_time": "2021-03-07 21:18:03"},
         {"ip": "121.205.214.33", "port": 4245, "expire_time": "2021-03-07 20:58:11"},
         {"ip": "117.69.201.179", "port": 4262, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "112.248.12.235", "port": 4282, "expire_time": "2021-03-07 19:27:46"},
         {"ip": "112.113.194.167", "port": 4242, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "42.59.117.14", "port": 4260, "expire_time": "2021-03-07 21:27:00"},
         {"ip": "106.111.157.5", "port": 4217, "expire_time": "2021-03-07 19:56:01"},
         {"ip": "183.165.10.244", "port": 4272, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "117.70.47.32", "port": 4210, "expire_time": "2021-03-07 18:56:44"},
         {"ip": "183.147.30.227", "port": 4274, "expire_time": "2021-03-07 19:08:02"},
         {"ip": "183.165.237.232", "port": 4210, "expire_time": "2021-03-07 21:50:50"},
         {"ip": "120.35.202.127", "port": 4278, "expire_time": "2021-03-07 20:46:12"},
         {"ip": "27.44.222.232", "port": 4278, "expire_time": "2021-03-07 21:50:50"}
    ]
    global browser
    # browser = set_firefox(proxy_lst.pop())
    browser = set_firefox()

    if os.path.exists('city_urls.json'):
        city_urls = json.load(open('city_urls.json'))
    else:
        city_urls = get_city_urls()
        with open('city_urls.json', 'w', encoding='utf-8') as f:
            json.dump(city_urls, f)

    target_city_lst = ['深圳', '成都', '湛江', '曲靖', '泸州', '连云港', '通化', '潍坊']
    if os.path.exists('success.list'):
        success = {x.strip() for x in open('success.list')}
    else:
        success = set()
    success_log_file = open('success.list', 'w')
    fw = open('result.txt', 'a+', encoding='utf-8')
    for city, city_url in city_urls.items():
        if city.strip('市') in target_city_lst:
            shou_url, zu_rul = get_zu_shou_url(city_url)
            if not os.path.exists(f'{city}.target.urls.txt'):
                shou_urls = get_detail_url(shou_url)
                zu_urls = get_detail_url(zu_rul)
                print(city, 'shou:', len(shou_urls))
                print(city, 'zu:', len(zu_urls))
                all_urls = shou_urls + zu_urls
                with open(f'{city}.target.urls.txt', 'w') as f:
                    _ = [f.write(x+'\n') for x in all_urls]
            else:
                all_urls = [x.strip() for x in open(f'{city}.target.urls.txt')]
                print(f'all_url_of_{city}:{len(all_urls)}')

            for url in all_urls:
                if url in success:
                    continue
                info = get_detail_info(url, city)
                if info == 'change_IP':
                    browser = set_firefox(proxy_lst.pop())
                    continue
                else:
                    success.add(url)
                    success_log_file.write(url+'\n')
                    fw.write(f'{info}\n')
                    # pprint(info)
    fw.close()
    success_log_file.close()

# run
pipeline()
