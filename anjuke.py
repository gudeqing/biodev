# coding=utf-8
import os
import itertools
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.remote.webelement import WebElement
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


def set_chrome(proxy=None):
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
        # "profile.managed_default_content_settings.images": 2,
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
    print('google')
    browser.get("http://httpbin.org/ip")
    print(browser.page_source)
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
    # 打开浏览器
    options = webdriver.FirefoxOptions()
    # options.add_argument('-headless')  # 无头参数
    options.add_argument("ignore-certificate-errors")
    options.add_argument("--ignore-ssl-errors")
    executable_path = 'D:\Drivers\geckodriver.exe'
    profile = webdriver.FirefoxProfile()
    profile.set_preference("network.proxy.type", 1)
    # 设置代理
    if proxy:
        ip, port = proxy['ip'], proxy['port']
        print(f'使用代理IP：{ip}, PORT：{port}')
        profile.set_preference("network.proxy.http", f"{ip}")
        profile.set_preference("network.proxy.http_port", f'{port}')
    else:
        print('未使用代理！')
    profile.update_preferences()
    browser = webdriver.Firefox(firefox_profile=profile, executable_path=executable_path, options=options)
    browser.get("http://httpbin.org/ip")
    print(browser.page_source)
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
    random_wait()
    zu_shou_sites = set()
    # for choice in browser.find_elements_by_class_name('third_navlist'):
    for each in browser.find_elements_by_tag_name('a'):
        sell_type_site = each.get_attribute('href')
        if sell_type_site.endswith('.anjuke.com/sp-zu/'):
            zu_shou_sites.add(sell_type_site)
        elif sell_type_site.endswith('.anjuke.com/sp-shou/'):
            zu_shou_sites.add(sell_type_site)
        if len(zu_shou_sites) == 2:
            break
    print(city_url)
    print('?', zu_shou_sites)
    return sorted(zu_shou_sites)


def get_detail_info(url, city=None):
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
        # print(e)
        try:
            pup_date = browser.find_element_by_class_name('time').text
        except Exception as e:
            print(e)
            pup_date = 'unknown'
    # 提取基础信息
    try:
        basic_info = browser.find_element_by_class_name('basic-info-wrapper').find_elements_by_class_name('item')
    except Exception as e:
        print(e)
        return 'change_IP'

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
    # 获取standardHouseId
    # standardHouseId = re.search('standardHouseId: (.*),?', page).groups()[0]
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
    tmp = re.findall('(\d+)m', basic_info_dict['规格'])
    if not tmp:
        tmp = ['unknown', 'unknown', 'unknown']
    elif len(tmp) == 2:
        tmp.append('unknown')
    elif len(tmp) == 1:
        tmp += ['unknown', 'unknown']
    else:
        tmp = ['unknown', 'unknown', 'unknown']
    width, height, level = tmp
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


def get_all_detail(browser, start_url, city, proxy_lst, outfile='result.txt'):
    browser.get(start_url)
    time.sleep(3)
    success_set = set()
    if os.path.exists(outfile):
        for line in open(outfile, encoding='utf-8'):
            tmp_dict = eval(line.strip())
            if tmp_dict.get('url') and len(tmp_dict['url']) > 10:
                success_set.add(tmp_dict['url'].split('?', 1)[0])
    fw = open(outfile, 'a+', encoding='utf-8')
    while True:
        indexes = list(range(len(browser.find_elements_by_class_name('list-item'))))
        for i in indexes:
            choice: WebElement = browser.find_elements_by_class_name('list-item')[i]
            # 提取detail_url，如果想进入具体页面抓取信息，则可能需要进行验证
            detail = choice.find_element_by_tag_name('a')
            try:
                detail.find_element_by_class_name('mianyi')
                continue
            except Exception:
                pass
            detail_url = detail.get_attribute('href')
            if detail_url.split('?', 1)[0] in success_set:
                print('已经爬取，跳过： ', detail_url)
                continue
            # 进入详情页
            detail.click()
            random_sleep()
            # 切换tab
            window_handles = browser.window_handles
            browser.switch_to.window(window_handles[-1])
            wait_wrapper(EC.visibility_of_element_located((By.CLASS_NAME, 'basic-info-wrapper')), browser)
            # 获取数据
            info = get_detail_info(detail_url, city)
            if info == 'change_IP':
                # 会跳过爬取失败的目录
                # browser = set_firefox(proxy_lst.pop())
                browser = set_chrome(proxy_lst.pop())
                indexes.append(i)
                continue
            # print(info)
            fw.write(f'{info}\n')
            # 关闭tab
            browser.close()
            assert len(browser.window_handles) == 1
            random_sleep()
            browser.switch_to.window(window_handles[0])
        # 点击下一页
        try:
            next_page = browser.find_element_by_class_name('aNxt')
            next_page.click()
            random_sleep()
            wait_wrapper(EC.visibility_of_element_located((By.CLASS_NAME, 'aNxt')), browser)
        except Exception as e:
            print(e)
            break

    fw.close()


def pipeline(proxy_and_city):
    global browser
    proxy_lst, target_city_lst, outfile = proxy_and_city
    browser = set_chrome(proxy_lst.pop())
    # browser = set_chrome()
    if os.path.exists('city_urls.json'):
        city_urls = json.load(open('city_urls.json'))
    else:
        city_urls = get_city_urls()
        with open('city_urls.json', 'w', encoding='utf-8') as f:
            json.dump(city_urls, f)

    for city, city_url in city_urls.items():
        if city.strip('市') in target_city_lst:
            shou_url, zu_rul = get_zu_shou_url(city_url)
            get_all_detail(browser, shou_url, city, proxy_lst, outfile)
            print(f'提取完{city}的商铺出售信息')
            get_all_detail(browser, zu_rul, city, proxy_lst, outfile)
            print(f'提取完成{city}的商铺出租信息')
    browser.close()


if __name__ == '__main__':
    # batch run
    from multiprocessing import Pool

    proxy_lst = [{"ip":"183.164.255.223","port":4252,"expire_time":"2021-03-10 02:41:01","city":"安徽省淮北市"},{"ip":"116.16.176.178","port":4245,"expire_time":"2021-03-10 03:47:50","city":"广东省揭阳市"},{"ip":"114.106.170.30","port":4245,"expire_time":"2021-03-10 01:57:07","city":"安徽省池州市"},{"ip":"171.41.129.111","port":4278,"expire_time":"2021-03-10 02:50:47","city":"湖北省荆门市"},{"ip":"59.58.150.136","port":4206,"expire_time":"2021-03-10 03:47:50","city":"福建省莆田市"},{"ip":"113.229.0.52","port":4256,"expire_time":"2021-03-10 03:47:50","city":"辽宁省营口市"},{"ip":"182.241.34.196","port":4240,"expire_time":"2021-03-10 03:29:43","city":"云南省大理白族自治州"},{"ip":"122.143.55.203","port":4278,"expire_time":"2021-03-10 01:53:00","city":"吉林省白城市"},{"ip":"113.117.26.160","port":4245,"expire_time":"2021-03-10 03:47:50","city":"广东省揭阳市"},{"ip":"114.99.12.237","port":4227,"expire_time":"2021-03-10 03:47:50","city":"安徽省铜陵市"},{"ip":"180.122.102.65","port":4216,"expire_time":"2021-03-10 03:47:50","city":"江苏省泰州市"},{"ip":"139.208.108.61","port":4278,"expire_time":"2021-03-10 03:21:56","city":"吉林省延边朝鲜族自治州"},{"ip":"114.104.140.70","port":4247,"expire_time":"2021-03-10 03:47:50","city":"安徽省黄山市"},{"ip":"49.87.124.170","port":4253,"expire_time":"2021-03-10 02:55:16","city":"江苏省淮安市"},{"ip":"218.91.133.204","port":4245,"expire_time":"2021-03-10 02:57:01","city":"江苏省扬州市"},{"ip":"14.134.188.48","port":4272,"expire_time":"2021-03-10 03:47:50","city":"宁夏固原市"},{"ip":"183.166.86.41","port":4251,"expire_time":"2021-03-10 03:47:50","city":"安徽省淮南市"},{"ip":"42.85.237.96","port":4256,"expire_time":"2021-03-10 00:56:24","city":"辽宁省营口市"},{"ip":"182.241.34.78","port":4240,"expire_time":"2021-03-10 03:47:50","city":"云南省大理白族自治州"},{"ip":"221.9.150.157","port":4258,"expire_time":"2021-03-10 03:47:50","city":"吉林省四平市"},{"ip":"120.38.18.105","port":4232,"expire_time":"2021-03-10 01:05:37","city":"福建省漳州市"},{"ip":"120.34.231.74","port":4216,"expire_time":"2021-03-10 01:54:31","city":"福建省南平市"},{"ip":"113.237.2.168","port":4256,"expire_time":"2021-03-10 02:23:05","city":"辽宁省营口市"},{"ip":"125.119.66.7","port":4245,"expire_time":"2021-03-10 03:47:50","city":"浙江省杭州市"},{"ip":"27.190.74.2","port":4278,"expire_time":"2021-03-10 01:34:15","city":"河北省唐山市"},{"ip":"60.175.39.178","port":4235,"expire_time":"2021-03-10 03:47:50","city":"安徽省安庆市"},{"ip":"59.60.153.55","port":4245,"expire_time":"2021-03-10 02:48:37","city":"福建省漳州市"}]

    target_cities_lst = [
        ['连云港', '通化', '潍坊'],
        ['湛江', '曲靖', '泸州'],
        ['深圳', '成都']
    ]
    target_cities_lst = [
        ['通化', '潍坊'],
        ['湛江', '曲靖'],
        ['深圳']
    ]
    proxy_num = int(len(proxy_lst)/len(target_cities_lst))
    print('proxy number for each run:', proxy_num)
    arguments = []
    for ind, cities in enumerate(target_cities_lst):
        proxies = [proxy_lst.pop() for x in range(proxy_num)]
        arguments.append((proxies, cities, f'{ind}.result.txt'))

    pool = Pool(3)
    pool.map(pipeline, arguments)
    pool.close()
    pool.join()

# # run
# pipeline()
# while True:
#     try:
#         pipeline()
#     except Exception as e:
#         # time.sleep(300)
#         print(e)
#         pipeline()

# 17521126310
# 1q2w3e4r5t
