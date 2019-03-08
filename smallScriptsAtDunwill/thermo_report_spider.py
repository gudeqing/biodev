# coding=utf-8
import requests
from bs4 import BeautifulSoup as bs
import json
import re
import os
from collections import OrderedDict
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
    print(plugin_result_dir)
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

    # 过滤overall info
    import pandas as pd
    final_overall = dict()
    for k, v in overall_info.items():
        if k.endswith('Args'):
            continue
        final_overall[k] = v

    # 过滤sample detail
    sample_detail_dict = dict()
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
        sample_detail_dict[x] = tmp_y

    overall_target_fields = "Run Name, Run Date, Report Name, Report Date, Research Application, Projects, Chip Type, " \
                           "Library Read Length, Flows, Number of Barcodes, Sample Tube Label, Chip Barcode, Chip Lot, " \
                           "Chip Exp., Library Con., Total Bases, Key Signal, ISP Loading, Enrichment, Library ISPs (%)," \
                           " Clonal, Final Library, Test Fragments, Adapter Dimer, Low Quality, Usable Reads, Total Reads," \
                           " Mean Read Len., Media Read Len., Mode Read Len., Total Aligned Bases, Aligned Bases (%), " \
                           "Unaligned Bases (%), Average Coverage Depth of Reference, TotalAlignedReads, AlignedReads, " \
                           "UnalignedReads, AlignedReads (%), UnalignedReads (%), Mean Raw Accuracy 1× (%), AQ17 Bases," \
                           " AQ20 Bases, Prefect Bases, AQ17_MeanLength (bp), AQ17_LongestAlignment (bp), " \
                           "AQ17_MeanCoverageDepth, AQ20_MeanLength (bp), AQ20_LongestAlignment (bp), " \
                           "AQ20_MeanCoverageDepth, Perfect_MeanLength (bp), Perfect_LongestAlignment (bp), " \
                           "Perfect_MeanCoverageDepth, AddressableWells, WithISPs, Live, TestFragment, LibraryISPs, " \
                           "FilteredPolyclonal, FilteredLowQuality, FilteredAdapterDimer, FinalLibraryISPs, " \
                           "FilteredPolyclonal_r, Sample, OperationMode, FlowOrder, LibraryKey, TFKey, ChipCheck, " \
                           "ChipData, ChipWafer, BarcodeSet, AnalysisFlows, runID"
    sample_target_fields = "Report Name, Barcode Name, 建库编号, 测序编号, Library Type, Library kit, Reference, " \
                           "Target regions, Bases, >=Q20 Bases (%), Reads, Mean Read Length (bp), Mapped Reads, " \
                           "Mapped Reads (%), On Target (%), Mean Depth, Mean Depth / Mapped Reads (%), Uniformity, " \
                           "Median Read Cov, Median Mol Cov, Media Mol Cov / Media Read Cov (%), " \
                           "Median Read Cov / Reads (%), Median Mol Cov / Reads (%), LOD %, Variants, Hotspot Variants," \
                           " Calculate Amplicons, Number of amplicons, Percent assigned amplicon reads," \
                           " Average reads per amplicon, Uniformity of amplicon coverage," \
                           " Amplicons with at least 1 read, Amplicons with at least 20 reads," \
                           " Amplicons with at least 100 reads, Amplicons with at least 500 reads," \
                           " Amplicons with no strand bias, Amplicons reading end-to-end," \
                           " Amplicon base composition bias, Bases in target regions, Percent base reads on target," \
                           " Average base coverage depth, Uniformity of base coverage, Target base coverage at 1x," \
                           " Target base coverage at 20x, Target base coverage at 100x, Target base coverage at 500x," \
                           " Target bases with no strand bias, Percent end-to-end reads"
    overall_target_fields = [x.strip() for x in overall_target_fields.split(',')]
    sample_target_fields = [x.strip() for x in sample_target_fields.split(',')]

    # def match_fields(my_fields, ref_fields:list):
    #     from fuzzywuzzy import fuzz
    #     from fuzzywuzzy import process
    #     choices = set(my_fields)
    #     match_dict = OrderedDict()
    #     for each in ref_fields:
    #         match, score = process.extractOne(each, choices)
    #         if score > 80:
    #             match_dict[each] = match
    #         else:
    #             match_dict[each] = 'unknown'
    #     return match_dict

    # pprint(match_fields(final_overall.keys(), overall_target_fields))
    # one_sample = list(sample_detail_dict.keys())[0]
    # pprint(match_fields(sample_detail_dict[one_sample].keys(), sample_target_fields))

    final_overall_match_dict = OrderedDict([
        ('Run Name', final_overall['Run Name']),
        ('Run Date', final_overall['Run Date']),
        ('Report Name', final_overall['Report Name']),
        ('Report Date', final_overall['Report Date']),
        ('Research Application', 'unknown'),
        ('Projects', final_overall['Projects']),
        ('Chip Type', final_overall['Chip Type']),
        ('Library Read Length', 'unknown'),
        ('Flows', 'unknown'),
        ('Number of Barcodes', 'unknown'),
        ('Sample Tube Label', final_overall['Sample Tube Label']),
        ('Chip Barcode', final_overall['Chip Barcode']),
        ('Chip Lot', final_overall['Chip Lot Number']),
        ('Chip Exp.', 'unknown'),
        ('Library Con.', 'unknown'),
        ('Total Bases', final_overall['Totabl Bases']),
        ('Key Signal', final_overall['Key Signal']),
        ('ISP Loading', final_overall['ISP Loading']),
        ('Enrichment', 'unknown'),
        ('Library ISPs (%)', 'unknown'),
        ('Clonal', 1- int(final_overall['Filtered: Polyclonal'])/int(final_overall['Library ISPs'])),
        ('Final Library', final_overall["Final Library ISPs"]),
        ('Test Fragments', final_overall['Test Fragment']),
        ('Adapter Dimer', final_overall['Filtered: Adapter Dimer']),
        ('Low Quality', final_overall['Filtered: Low Quality']),
        ('Usable Reads', final_overall['Usable Reads']),
        ('Total Reads', final_overall['Total Reads']),
        ('Mean Read Len.', final_overall['Mean']),
        ('Media Read Len.', final_overall['Median']),
        ('Mode Read Len.', final_overall['Mode']),
        ('Total Aligned Bases', final_overall['Total Aligned Bases']),
        ('Aligned Bases (%)', 'unknown'),
        ('Unaligned Bases (%)', 'unknown'),
        ('Average Coverage Depth of Reference', 'unknown'),
        ('TotalAlignedReads', 'unknown'),
        ('AlignedReads', final_overall['Aligned Reads']),
        ('UnalignedReads', final_overall['Unaligned Reads']),
        ('AlignedReads (%)', int(final_overall['Aligned Reads'])/int(final_overall['Total Reads'])),
        ('UnalignedReads (%)', 1-int(final_overall['Aligned Reads'])/int(final_overall['Total Reads'])),
        ('Mean Raw Accuracy 1× (%)', final_overall['Mean Raw Accuracy 1x']),
        ('AQ17 Bases', final_overall['AQ17 Total Bases']),
        ('AQ20 Bases', final_overall["AQ20_TotalNumberofBases[bp]"]),
        ('Prefect Bases', final_overall["Perfect_TotalNumberofBases[bp]"]),
        ('AQ17_MeanLength (bp)', final_overall['AQ17_MeanLength[bp]']),
        ('AQ17_LongestAlignment (bp)', final_overall['AQ17_LongestAlignment[bp]']),
        ('AQ17_MeanCoverageDepth', final_overall['AQ17_MeanCoverageDepth[x]']),
        ('AQ20_MeanLength (bp)', final_overall['AQ20_MeanLength[bp]']),
        ('AQ20_LongestAlignment (bp)', final_overall['AQ20_LongestAlignment[bp]']),
        ('AQ20_MeanCoverageDepth', final_overall['AQ20_MeanCoverageDepth[x]']),
        ('Perfect_MeanLength (bp)', final_overall['Perfect_MeanLength[bp]']),
        ('Perfect_LongestAlignment (bp)', final_overall['Perfect_LongestAlignment[bp]']),
        ('Perfect_MeanCoverageDepth', final_overall['Perfect_MeanCoverageDepth[x]']),
        ('AddressableWells', final_overall['Addressable Wells']),
        ('WithISPs', final_overall['With ISPs']),
        ('Live', final_overall['Live']),
        ('TestFragment', final_overall['Test Fragment']),
        ('LibraryISPs', final_overall['Library ISPs']),
        ('FilteredPolyclonal', final_overall['Filtered: Polyclonal']),
        ('FilteredLowQuality', final_overall["Filtered: Low Quality"]),
        ('FilteredAdapterDimer', final_overall["Filtered: Adapter Dimer"]),
        ('FinalLibraryISPs', final_overall["Final Library ISPs"]),
        ('FilteredPolyclonal_r', 'unknown'),
        ('Sample', 'unknown'),
        ('OperationMode', final_overall['Sequencer Operation Mode']),
        ('FlowOrder', final_overall['Flow Order']),
        ('LibraryKey', final_overall['Library Key']),
        ('TFKey', final_overall['TF Key']),
        ('ChipCheck', final_overall['Chip Check']),
        ('ChipData', final_overall['Chip Data']),
        ('ChipWafer', final_overall['Chip Wafer']),
        ('BarcodeSet', final_overall['Barcode Set']),
        ('AnalysisFlows', final_overall['Analysis Flows']),
        ('runID', final_overall['runID'])
    ])

    pd.DataFrame(final_overall_match_dict, index=['all']).to_csv('overall.info.xls', index=False, header=True, sep='\t')

    final_sample_detail = list()
    for x, y in sample_detail_dict.items():
        y = OrderedDict([
            ('Report Name', final_overall['Report Name']),
            ('Barcode Name', y['barcode_name']),
            (u'建库编号', 'unknown'),
            (u'测序编号', 'unknown'),
            ('Library Type', 'unknown'),
            ('Library kit', 'unknown'),
            ('Reference', y['reference']),
            ('Target regions', y['Target Regions']),
            ('Bases', y["total_bases"]),
            ('>=Q20 Bases (%)', '{}%'.format(round(int(y['Q20_bases'])/int(y["total_bases"])*100, 2))),
            ('Reads', 'unknown'),
            ('Mean Read Length (bp)', y["mean_read_length"]),
            ('Mapped Reads', y['mapped_reads']),
            ('Mapped Reads (%)', 'unknown'),
            ('On Target (%)', y['Percent base reads on target']),
            ('Mean Depth', y['mean_depth']),
            ('Mean Depth / Mapped Reads (%)', 'unknown'),
            ('Uniformity', y['uniformity']),
            ('Median Read Cov', y['Median Read Cov']),
            ('Median Mol Cov', y['Median Mol Cov']),
            ('Media Mol Cov / Media Read Cov (%)', 'unknown'),
            ('Median Read Cov / Reads (%)', 'unknown'),
            ('Median Mol Cov / Reads (%)', 'unknown'),
            ('LOD %', y['LOD %']),
            ('Variants', y['variants']),
            ('Hotspot Variants', y['hotspot_variants']),
            ('Calculate Amplicons', 'unknown'),
            ('Number of amplicons', y['Number of amplicons']),
            ('Percent assigned amplicon reads', y['Percent assigned amplicon reads']),
            ('Average reads per amplicon', y['Average reads per amplicon']),
            ('Uniformity of amplicon coverage', y['Uniformity of amplicon coverage']),
            ('Amplicons with at least 1 read', y['Amplicons with at least 1 read']),
            ('Amplicons with at least 20 reads', y['Amplicons with at least 20 reads']),
            ('Amplicons with at least 100 reads', y['Amplicons with at least 100 reads']),
            ('Amplicons with at least 500 reads', y['Amplicons with at least 500 reads']),
            ('Amplicons with no strand bias', y['Amplicons with no strand bias']),
            ('Amplicons reading end-to-end', y['Amplicons reading end-to-end']),
            ('Amplicon base composition bias', y['Amplicon base composition bias']),
            ('Bases in target regions', y['Bases in target regions']),
            ('Percent base reads on target', y['Percent base reads on target']),
            ('Average base coverage depth', y['Average base coverage depth']),
            ('Uniformity of base coverage', y['Uniformity of base coverage']),
            ('Target base coverage at 1x', y['Target base coverage at 1x']),
            ('Target base coverage at 20x', y['Target base coverage at 20x']),
            ('Target base coverage at 100x', y['Target base coverage at 100x']),
            ('Target base coverage at 500x', y['Target base coverage at 500x']),
            ('Target bases with no strand bias', y['Target bases with no strand bias']),
            ('Percent end-to-end reads', y['Percent end-to-end reads'])
        ])
        final_sample_detail.append(y)
    pd.DataFrame(final_sample_detail).to_csv('sample_detail.xls', header=True, index=False, sep='\t')


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
