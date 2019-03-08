import requests
from bs4 import BeautifulSoup as bs
import json
import re
import os
from pprint import pprint


def download(url='http://10.62.2.15/report/146/', login_url='http://10.62.2.15/login/'):
    if not url:
        return None
    user_agent = 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'
    headers = {'User-Agenet': user_agent}
    payload = {
        'username': 'ionadmin',
        'password': 'ionadmin'
    }
    # Use 'with' to ensure the session context is closed after use.
    with requests.Session() as s:
        p = s.post(login_url, data=payload, headers=headers)
        # An authorised request.
        r = s.get(url, headers=headers)
        if r.status_code == 200:
            r.encoding = 'utf-8'
            return r.text


def main(url='http://10.62.2.15/report/146/', login_url='http://10.62.2.15/login/'):
    page = download(url, login_url)
    soup = bs(page, "html.parser")
    overall_info = dict()
    sample_detail = dict()

    run_summary = soup.find('div', id="RunSummary")
    report_labels = run_summary.find_all(class_='report-label')
    summary_dict = dict()
    for tag in report_labels:
        key = tag.string
        if tag.string == 'Projects':
            value = tag.next_sibling.next_sibling.contents[1].string
        else:
            value = tag.next_sibling.next_sibling.string
        summary_dict[key] = value
    with open('summary.json', 'w') as f:
        json.dump(summary_dict, f, indent=2)
        overall_info.update(summary_dict)

    bead = soup.find('div', id="beadfind")
    bead_dict = dict()
    bead_dict['Totabl Bases'] = bead.find('small', string="Total Bases").parent.find('h2').string
    bead_dict['Key Signal'] = bead.find('small', string='Key Signal').parent.find('h2').string
    bead_dict['ISP Loading'] = bead.find('div', id="bead-signal")['data-beadloading']+'%'
    # print(bead_dict)

    basecall = soup.find('div', id='basecaller')
    basecall_dict = dict()
    basecall_dict['Total Reads'] = basecall.find('small', string="Total Reads").parent.find('h2').string.replace(',', '')
    basecall_dict['Usable Reads'] = basecall.find('div', id="usable_sequence")['data-usablesequence']+'%'
    # print(basecall_dict)

    read_length = dict()
    read_len_info = soup.find('div', id="readlength")
    for each in ['Mean', 'Median', 'Mode']:
        read_length[each] = read_len_info.find('small', string=each).parent.find('h2').string
    # print(read_length)

    with open('unaligned_reads.json', 'w') as f:
        tmp_dict = dict()
        tmp_dict.update(bead_dict)
        tmp_dict.update(basecall_dict)
        tmp_dict.update(read_length)
        json.dump(tmp_dict, f, indent=2)
        overall_info.update(tmp_dict)

    Chip_well_details = dict()
    table = soup.find('strong', string="Chip well details").next_sibling.next_sibling.get_text()
    table = re.sub('[\n]+', '\n', table).strip().split('\n')
    Chip_well_details[table[0]] = table[1].replace(',', '')
    Chip_well_details.update(dict(zip(table[2:][::3], table[3:][::3])))
    Chip_well_details = {k: int(v.replace(',', '')) for k,v in Chip_well_details.items()}
    # print(Chip_well_details)
    with open('unaligned_reads.Chip_well_details.json', 'w') as f:
        json.dump(Chip_well_details, f, indent=2)
        overall_info.update(Chip_well_details)


    Library_ISP_details = dict()
    table = soup.find('strong', string="Library ISP details").next_sibling.next_sibling.get_text()
    table = re.sub('[\n]+', '\n', table).strip().split('\n')
    Library_ISP_details[table[0]] = table[1].replace(',', '')
    Library_ISP_details.update(dict(zip(table[2:][::3], table[3:][::3])))
    Library_ISP_details = {k: int(v.replace(',', '')) for k,v in Library_ISP_details.items()}
    # print(Library_ISP_details)
    with open('unaligned_reads.Library_ISP_details.json', 'w') as f:
        json.dump(Library_ISP_details, f, indent=2)
        overall_info.update(Library_ISP_details)


    alignment_stat = dict()
    align = soup.find('div', id='alignMap')
    alignment_stat['Total Aligned Bases'] = align.find('small', string="Total Aligned Bases").parent.find('h2').string
    alignment_stat['Reference Coverage'] = align.find('small', string="Reference Coverage").parent.find('h2').string
    table = align.find('table').get_text()
    table = re.sub('[\n]+', '\n', table).strip().split('\n')
    tmp_dict = dict(zip(table[2:][::3], table[3:][::3]))
    tmp_dict = {k: int(v.replace(',', '')) for k, v in tmp_dict.items()}
    alignment_stat.update(tmp_dict)
    # raw aligned
    align = soup.find('div', id='rawAligned')
    alignment_stat['Mean Raw Accuracy 1x'] = align.find('small', string="Mean Raw Accuracy 1x").parent.find('h2').string
    # alignment
    align = soup.find('div', id='alignment')
    alignment_stat['AQ17 Total Bases'] = align.find('small', string="AQ17 Total Bases").parent.find('h2').string
    # Alignment Quality table
    table = align.find('table').get_text()
    table = re.sub('[\n]+', '\n', table).strip().split('\n')
    header = table[:3]
    for ind in range(3, len(table), 4):
        # alignment_stat[table[ind]] = dict(zip(header, (x.strip() for x in table[ind+1:ind+4])))
        for k, v in dict(zip(header, (x.strip() for x in table[ind+1:ind+4]))).items():
            alignment_stat[k+'_'+table[ind].replace(' ', '')] = v
    # pprint(alignment_stat)
    with open('aligned_reads.json', 'w') as f:
        json.dump(alignment_stat, f, indent=2)
        overall_info.update(alignment_stat)

    # var barcodes_json
    json_content = re.findall('var barcodes_json = (.*?);', page)
    if json_content:
        json_content = json.loads(json_content[0])
    # pprint(json_content)
    with open('output_files.json', 'w') as f:
        json.dump(json_content, f, indent=2)
        for each in json_content:
            if 'barcode_name' in each:
                sample_detail.setdefault(each['barcode_name'], dict())
                sample_detail[each['barcode_name']].update(each)


    # detail
    test_fragments = soup.find('div', id="TestFragments")
    table = test_fragments.find('table').get_text()
    table = re.sub('[\n]+', '\n', table).strip().split('\n')
    test_fragments_dict = dict(zip(table[:4], table[6:]))
    test_fragments_dict['Test Fragment'] = test_fragments_dict['Test Fragment'].replace(',', '')
    # print(test_fragments_dict)
    with open('TestFragments.json', 'w') as f:
        json.dump(test_fragments_dict, f, indent=2)
        overall_info.update(test_fragments_dict)

    chef_summary_dict = dict()
    chef_summary = soup.find('div', id="ChefSummary")
    table = chef_summary.find('table')
    rows = table.find_all('tr')
    for row in rows:
        k, v = list(row.children)
        chef_summary_dict[k.string] = v.string
    # pprint(chef_summary_dict)
    with open('ChefSummary.json', 'w') as f:
        json.dump(chef_summary_dict, f, indent=2)
        overall_info.update(chef_summary_dict)

    # S5ConsumableSummary
    S5ConsumableSummary = dict()
    s5 = soup.find('div', id="S5ConsumableSummary")
    ps = s5.find_all('p')
    for p in ps:
        S5ConsumableSummary[p.string.split(':')[0]] = p.string.split(':')[1].strip()
    # print(S5ConsumableSummary)
    table = s5.find('table')
    rows = list(table.find_all('tr'))
    header = [x.string.strip() for x in rows[0].children if x.string.strip()]
    for row in rows[1:]:
        row_data = [x.string.strip() for x in row.children if x.string.strip()]
        S5ConsumableSummary[row_data[0]] = dict(zip(header[1:], row_data[1:]))
    # pprint(S5ConsumableSummary)
    with open('S5ConsumableSummary.json', 'w') as f:
        json.dump(S5ConsumableSummary, f, indent=2)

    #
    AnalysisDetails_dict = dict()
    detail = soup.find('div', id="AnalysisDetails")
    table = detail.find('table')
    rows = table.find_all('tr')
    for row in rows:
        row = bs(str(row).replace('\n', ''), "html.parser")
        k, v = [x.string.strip() for x in row.find_all('td')]
        AnalysisDetails_dict[k] = v
    # pprint(AnalysisDetails_dict)
    with open('AnalysisDetails.json', 'w') as f:
        json.dump(AnalysisDetails_dict, f, indent=2)
        overall_info.update(AnalysisDetails_dict)

    SoftwareVersion_dict = dict()
    version = soup.find('div', id="SoftwareVersion")
    table = version.find('table')
    rows = table.find_all('tr')
    for row in rows:
        row = bs(str(row).replace('\n', ''), "html.parser")
        k, v = [x.string.strip() for x in row.find_all('td')]
        SoftwareVersion_dict[k] = v
    # pprint(SoftwareVersion_dict)
    with open('SoftwareVersion.json', 'w') as f:
        json.dump(SoftwareVersion_dict, f, indent=2)
        overall_info.update(SoftwareVersion_dict)

    # 连接服务器并且找到目标html
    lst = re.findall(r'/output/Home/Auto_user.*?/', page)
    plugin_result_dir = '/results/analysis' + lst[0] + 'plugin_out/'
    import paramiko, os
    ssh = paramiko.SSHClient()  # 创建SSH对象
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())  # 允许连接不在know_hosts文件中的主机
    ssh.connect(hostname='10.62.2.15', port=22, username='ionadmin', password='ionadmin')  # 连接服务器
    stdin, stdout, stderr = ssh.exec_command('ls {}/*/*html'.format(plugin_result_dir))  # 执行命令并获取命令结果
    stdin2, stdout2, stderr2 = ssh.exec_command('ls {}/coverageAnalysis*/*/*.stats.cov.txt'.format(plugin_result_dir))  # 执行命令并获取命令结果
    # stdin为输入的命令
    # stdout为命令返回的结果
    # stderr为命令错误时返回的结果
    res, err = stdout.read(), stderr.read()
    result = res if res else err
    target_html_lst = result.strip().split()
    target_html_lst += stdout2.read().strip().split()
    # 下载目标文件
    transport = paramiko.Transport(('10.62.2.15', 22))
    transport.connect(username='ionadmin', password='ionadmin')
    sftp = paramiko.SFTPClient.from_transport(transport)
    # 目标文件下载到本地 local_path
    if not os.path.exists('plugin.html'):
        os.mkdir('plugin.html')
    for each in target_html_lst:
        sftp.get(each, os.path.join(b'plugin.html', os.path.basename(each)))
    transport.close()

    # coverageAnalysis.html
    target_html_lst = os.listdir('plugin.html')
    if 'coverageAnalysis.html' in target_html_lst:
        page = open(os.path.join('plugin.html', 'coverageAnalysis.html')).read()
        json_content = re.findall('var barcodes_json = (.*?);', page)
        if json_content:
            json_content = json.loads(json_content[0])
            with open('coverageAnalysis.json', 'w') as f:
                json.dump(json_content, f, indent=2)
                for each in json_content:
                    if 'barcode_name' in each:
                        sample_detail.setdefault(each['barcode_name'], dict())
                        sample_detail[each['barcode_name']].update(each)
    from glob import glob
    cov_stat = glob('plugin.html/*.stats.cov.txt')
    json_content = list()
    if cov_stat:
        for each in cov_stat:
            tmp_dict = dict()
            tmp_dict['barcode_name'] = os.path.basename(each).split('_Auto')[0]
            with open(each) as f:
                _ = f.readline()
                for line in f:
                    if line.strip():
                        k, v = [x.strip() for x in line.strip().split(':')]
                        tmp_dict[k] = v
            json_content.append(tmp_dict)
        with open('stats.cov.json', 'w') as f:
            json.dump(json_content, f, indent=2)
            for each in json_content:
                if 'barcode_name' in each:
                    sample_detail.setdefault(each['barcode_name'], dict())
                    sample_detail[each['barcode_name']].update(each)


    # variantCaller.html
    if 'variantCaller.html' in target_html_lst:
        page = open(os.path.join('plugin.html',  'variantCaller.html')).read()
        json_content = re.findall('barcode_json\[.*?\] = .*?;', page)
        if json_content:
            result = list()
            for each in json_content:
                if "barcode_json['index']" in each or "barcode_json['exports']" in each:
                    continue
                if each.startswith("barcode_json['barcode_details']"):
                    barcode_json = dict()
                    result.append(barcode_json)
                exec(each)
            with open('variantCaller.json', 'w') as f:
                json.dump(result, f, indent=2)
                json_content = result
                for each in json_content:
                    if 'barcode_name' in each:
                        sample_detail.setdefault(each['barcode_name'], dict())
                        sample_detail[each['barcode_name']].update(each)

    # write out overall info
    import pandas as pd
    pd.DataFrame(overall_info, index=['all']).to_csv('overall.info.xls', index=False, header=True, sep='\t')

    tmp_dict = dict()
    for x, y in sample_detail.items():
        tmp_y = dict()
        for k, v in y.items():
            if type(v) in [dict, list]:
                continue
            if k.endswith('_link'):
                continue
            if k.endswith('barcode_details'):
                continue
            if k == 'median_num_fam3':
                k = 'Median Mol Cov'
            elif k == 'median_depth':
                k = 'Median Read Cov'
            elif k == 'fm3_pass80':
                k = 'LOD %'
            tmp_y[k] = v
        tmp_dict[x] = tmp_y

    pd.DataFrame(tmp_dict).to_csv('sample_detail.xls', header=True, index=True, sep='\t')


if __name__ == '__main__':
    class Func2Command(object):
        def __init__(self, callable_dict):
            self.call(callable_dict)

        @staticmethod
        def introduce_command(func):
            import argparse
            import inspect
            import json
            import time
            if isinstance(func, type):
                description = func.__init__.__doc__
            else:
                description = func.__doc__
            if description:
                _ = [print(x.strip()) for x in description.split('\n') if x.strip()]
                parser = argparse.ArgumentParser(add_help=False)
            else:
                parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=description)
            func_args = inspect.getfullargspec(func)
            arg_names = func_args.args
            if not arg_names:
                func()
                return
            arg_defaults = func_args.defaults
            if not arg_defaults:
                arg_defaults = []
            arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
            sig = inspect.signature(func)
            for arg, value in zip(arg_names, arg_defaults):
                arg_type = sig.parameters[arg].annotation
                if arg == 'self':
                    continue
                if value == 'None':
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=True, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=True, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=True, metavar=arg)
                elif type(value) == bool:
                    if value:
                        parser.add_argument('--'+arg, action="store_false", help='default: True')
                    else:
                        parser.add_argument('--'+arg, action="store_true", help='default: False')
                elif value is None:
                    parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
                else:
                    if sig.parameters[arg].annotation in [list, tuple]:
                        parser.add_argument('-' + arg, default=value, nargs='+', type=type(value), metavar='Default:' + str(value), )
                    else:
                        parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
            if func_args.varargs is not None:
                print("warning: *varargs is not supported, and will be neglected! ")
            if func_args.varkw is not None:
                print("warning: **keywords args is not supported, and will be neglected! ")
            args = parser.parse_args().__dict__
            try:
                with open("Argument_detail.json", 'w') as f:
                    json.dump(args, f, indent=2, sort_keys=True)
            except IOError:
                print('Current Directory is not writable, thus argument log is not written !')
            start = time.time()
            func(**args)
            print("total time: {}s".format(time.time() - start))

        def call(self, callable_dict):
            import sys
            excludes = ['introduce_command', 'Func2Command']
            _ = [callable_dict.pop(x) for x in excludes if x in callable_dict]
            if len(callable_dict) >= 2:
                if len(sys.argv) <= 1:
                    print("The tool has the following sub-commands: ")
                    _ = [print(x) for x in callable_dict]
                    return
                sub_cmd = sys.argv[1]
                sys.argv.remove(sub_cmd)

                if sub_cmd in callable_dict:
                    self.introduce_command(callable_dict[sub_cmd])
                else:
                    print('sub-command: {} is not defined'.format(sub_cmd))
            elif len(callable_dict) == 1:
                self.introduce_command(callable_dict.pop(list(callable_dict.keys())[0]))

    callable_dict = {x: y for x, y in locals().items() if callable(y)}
    _ = [callable_dict.pop(x) for x in {'Func2Command'} if x in callable_dict]
    Func2Command(callable_dict)
