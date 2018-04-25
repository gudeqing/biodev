# coding=utf-8
import argparse
import subprocess
import pandas as pd
import re


def run_hmmscan(args):
    other_args = "".join(args.other_args)
    # format cmd and run cmd
    cmd = "./hmmscan -o {} --tblout {} --domtblout {} --seed 0 --notextw ".format(
        args.o, args.tblout, args.domtblout
    )
    cmd += "--cpu {} ".format(args.cpu)
    if "-T " not in args.other_args:
        cmd += "-E {} ".format(args.E)
        cmd += "--domE {} ".format(args.domE)
    if args.other_args:
        cmd += other_args
    cmd += " {hmm_db} {fasta}".format(hmm_db=args.hmmdb, fasta=args.seqfile)
    subprocess.check_call(cmd, shell=True)

    # reformat domtblout
    with open(args.domtblout) as fr, open("domain_predict.txt", 'w') as fw:
        for line in fr:
            if line.startswith("# target name"):
                col_names = [
                    "domain_family",
                    "pfam_id",
                    "domain_len",
                    "query_id",
                    "query_accession",
                    "query_len",
                    "e_value",
                    "score",
                    "bias",
                    "domain_order",
                    "domain_num",
                    "c_evalue",
                    "i_evalue",
                    "score",
                    "bias",
                    "hmm_start",
                    "hmm_end",
                    "ali_start",
                    "ali_end",
                    "env_start",
                    "env_end",
                    "acc",
                    "domain_description",
                ]
                fw.write("\t".join(col_names))
            elif line.startswith("#"):
                continue
            else:
                fw.write(line.replace(" ", '\t', 22))


def get_domain_stat(domain_predict):
    domain_predict_pd = pd.read_table(domain_predict, sep='\t', header=0)
    domain2query = domain_predict_pd.loc[:, ['pfam_id', 'query_id']]
    domain2query.columns = ['query_id', 'pfam_id']
    domain_stat_dict = dict()
    group_domain = domain2query.groupby("query_id").groups
    for key in group_domain:
        tmp_pd = domain2query.iloc[group_domain[key], :]
        tmp_count = tmp_pd.groupby("pfam_id").count().to_dict()["query_id"]
        new_tmp_count = dict()
        for k, v in tmp_count.items():
            new_tmp_count[k.split(".")[0]] = v
        domain_stat_dict[key] = new_tmp_count
    return domain_stat_dict


def judge_plant_tf(domain_stat_dict, rule_dict_list):
    """
    :param domain_stat_dict:
    :param rule_dict_list:
    :return: {query_id:
                    {
                        family_name1: {d1,d2},
                        family_name2: {d2,d3},
                    },
            ......
            }
    对于没有参与判定转录因子家族的domain将不予记录
    对于没有被判定为转录因子的query将不予记录
    """
    query2family = dict()
    for query_id, domain_dict in domain_stat_dict.items():
        predict_domains = set(domain_dict.keys())
        query2family[query_id] = dict()
        for each_tf in rule_dict_list:
            judge_rule = each_tf["judge_rule"]
            binding_domains = set(judge_rule['binding'].keys())
            auxiliary_domains = set(judge_rule["auxiliary"])
            forbidden_domains = set(judge_rule["forbidden"])
            if len(predict_domains & forbidden_domains) != 0:
                continue
            if len(binding_domains - predict_domains) != 0:
                continue
            range_match = [1 for x in binding_domains if domain_dict[x] in judge_rule['binding'][x]]
            if sum(range_match) == len(binding_domains):
                if len(auxiliary_domains & predict_domains) >= 1:
                    if each_tf['Family'] in query2family[query_id]:
                        query2family[query_id].pop(each_tf['Family'])
                    based_domains = binding_domains | (auxiliary_domains & predict_domains)
                    query2family[query_id][each_tf['SubFamily']] = based_domains
                else:
                    if each_tf['SubFamily'] not in query2family[query_id]:
                        query2family[query_id][each_tf['Family']] = binding_domains
            else:
                continue
        else:
            if not query2family[query_id]:
                _ = "{} was not judged as any kind of TF".format(query_id)
                query2family.pop(query_id)
    return query2family


def judge_animal_tf(domain_stat_dict, rule_dict):
    query2family = dict()
    for query_id, domain_dict in domain_stat_dict.items():
        query2family[query_id] = dict()
        for domain in domain_dict:
            family = rule_dict.get(domain)
            if family:
                query2family[query_id][family] = {domain, }
        else:
            if not query2family[query_id]:
                _ = "{} was not judged as any kind of TF".format(query_id)
                query2family.pop(query_id)
    else:
        return query2family


def get

def run_blast(**kwargs):
    """
    比对的目的是为了能够分配一个已知的TF_id，后续靶基因预测可以根据这个已知的TF_id获得转录因子靶向Motif.
    :param kwargs:
    :return:
    """


# tf_factor_rule
plant_tf_family_assignment_rules = """
Family	SubFamily	DNA-binding domain	Auxiliary domain	Forbidden domain
AP2/ERF	AP2	AP2 (>=2) (PF00847)	-	-
AP2/ERF	ERF	AP2 (1) (PF00847)	-	-
AP2/ERF	RAV	AP2 (PF00847) and B3 (PF02362)	-	-
B3 superfamily	ARF	B3 (PF02362)	Auxin_resp (PF06507)	-
B3 superfamily	B3	B3 (PF02362)	-	-
BBR-BPC	BBR-BPC	GAGA_bind (PF06217)	-	-
BES1	BES1	DUF822 (PF05687)	-	-
bHLH	bHLH	HLH (PF00010)	-	-
bZIP	bZIP	bZIP_1 (PF00170)	-	-
C2C2	CO-like	zf-B_box (PF00643)	CCT (PF06203)	-
C2C2	Dof	Zf-Dof (PF02701)	-	-
C2C2	GATA	GATA-zf (PF00320)	-	-
C2C2	LSD	Zf-LSD1 (PF06943)	-	Peptidase_C14 (PF00656)
C2C2	YABBY	YABBY (PF04690)	-	-
C2H2	C2H2	zf-C2H2 (PF00096)	-	RNase_T (PF00929)
C3H	C3H	Zf-CCCH (PF00642)	-	RRM_1 (PF00076) or Helicase_C (PF00271)
CAMTA	CAMTA	CG1 (PF03859)	-	-
CPP	CPP	TCR (PF03638)	-	-
DBB	DBB	zf-B_box (>=2) (PF00643)	-	-
E2F/DP	E2F/DP	E2F_TDP (PF02319)	-	-
EIL	EIL	EIN3 (PF04873)	-	-
FAR1	FAR1	FAR1 (PF03101)	-	-
GARP	ARR-B	G2-like (self-build)	Response_reg (PF00072)	-
GARP	G2-like	G2-like (self-build)	-	-
GeBP	GeBP	DUF573 (PF04504)	-	-
GRAS	GRAS	GRAS (PF03514)	-	-
GRF	GRF	WRC (PF08879)	QLQ (PF08880)	-
HB	HD-ZIP	Homeobox (PF00046)	HD-ZIP_I/II (self-build) or SMART (PF01852)	-
HB	TALE	Homeobox (PF00046)	BELL (self-build)or ELK (PF03789)	-
HB	WOX	homeobox (PF00046)	Wus type homeobox (self-build)	-
HB	HB-PHD	homeobox (PF00046)	PHD (PF00628)	-
HB	HB-other	homeobox (PF00046)	-	-
HRT-like	HRT-like	HRT-like (self-build)	-	-
HSF	HSF	HSF_dna_bind (PF00447)	-	-
LBD (AS2/LOB)	LBD (AS2/LOB)	DUF260 (PF03195)	-	-
LFY	LFY	FLO_LFY (PF01698)	-	-
MADS	M_type	SRF-TF (PF00319)	-	-
MADS	MIKC	SRF-TF (PF00319)	K-box (PF01486)	-
MYB superfamily	MYB	Myb_dna_bind (>=2) (PF00249)	-	SWIRM (PF04433)
MYB superfamily	MYB_related	Myb_dna_bind (1) (PF00249)	-	SWIRM (PF04433)
NAC	NAC	NAM (PF02365)	-	-
NF-X1	NF-X1	Zf-NF-X1 (PF01422)	-	-
NF-Y	NF-YA	CBFB_NFYA (PF02045)	-	-
NF-Y	NF-YB	NF-YB (self-build)	-	-
NF-Y	NF-YC	NF-YC (self-build)	-	-
Nin-like	Nin-like	RWP-RK (PF02042)	-	-
NZZ/SPL	NZZ/SPL	NOZZLE (PF08744)	-	-
S1Fa-like	S1Fa-like	S1FA (PF04689)	-	-
SAP	SAP	SAP (self-build)	-	-
SBP	SBP	SBP (PF03110)	-	-
SRS	SRS	DUF702 (PF05142)	-	-
STAT	STAT	STAT (self-build)	-	-
TCP	TCP	TCP (PF03634)	-	-
Trihelix	Trihelix	Trihelix (self-build)	-	-
VOZ	VOZ	VOZ (self-build)	-	-
Whirly	Whirly	Whirly (PF08536)	-	-
WRKY	WRKY	WRKY (PF03106)	-	-
ZF-HD	ZF-HD	ZF-HD_dimer (PF04770)	-	-
"""

animal_tf_family_assignment_rules = """
Family	DNA-binding domain	Pfam ID
AF-4	AF-4	PF05110
ARID	ARID	PF01388
bHLH	HLH	PF00010
CBF	CBF_alpha	PF02312
CEP-1	CEP1-DNA_bind	PF09287
CSL	BTD	PF09270
NF-YA	CBFB_NFYA	PF02045
CG-1	CG-1	PF03859
CP2	CP2	PF04516
CSD	CSD	PF00313
E2F	E2F_TDP	PF02319
ETS	Ets	PF00178
Fork head	Fork_head	PF00250
GCM	GCM	PF03615
GTF2I	GTF2I	PF02946
HMG	HMG_box	PF00505
HSF	HSF_DNA-bind	PF00447
HTH	HTH_psq	PF05225
IRF	IRF	PF00605
MYB	Myb_DNA-bd	PF00249
MBD	MBD	PF01429
NCU-G1	NCU-G1	PF15065
NDT80/PhoG	NDT80_PhoG	PF05224
Nrf1	Nrf1_DNA-bind	PF10491
PC4	PC4	PF02229
P53	P53	PF00870
PAX	PAX	PF00292
HPD	HPD	PF05044
RFX	RFX	PF02257
RHD	RHD	PF00554
Runt	Runt	PF00853
SAND	SAND	PF01342
SRF	SRF	PF00319
STAT	STAT_bind	PF02864
T-box	T-box	PF00907
TEA	TEA	PF01285
TSC22	TSC22	PF01166
Tub	Tub	PF01167
CTF/NFI	MH1	PF00859
MH1	MH1	PF03165
Homeobox	Homeobox	PF00046
Pou	Homeobox, Pou	PF00157
CUT	Homeobox, CUT	PF02376
TF_Otx	Homeobox, TF_Otx	PF03529
zf-C2HC	zf-C2HC	PF01530
zf-GAGA	zf-GAGA	PF09237
zf-BED	zf-BED	PF02892
ZBTB	zf-C2H2	PF00651
zf-C2H2	zf-C2H2	PF00096
DM	DM	PF00751
zf-GATA	zf-GATA	PF00320
zf-LITAF-like	zf-LITAF-like	PF10601
zf-MIZ	zf-MIZ	PF02891
zf-NF-X1	zf-NF-X1	PF01422
THAP	THAP	PF05485
"""


def parse_plant_tf_judge_rules(raw_rules):
    """
    最后依据judge_rule_dict的转录因子判定规则:
        1. key包含的所有pfam_id都要被蛋白包含，且每个pfam_id的个数在指定的列表范围，满足此条件后才可进入2，3判读
        2. 对于auxiliary列表包含的pfam_id，如某蛋白含有其中任何一个pfam_id且满足1条件，判定具体到subfamily，否则为母类。
        3. 对于forbidden列表包含的pfam_id， 如某蛋白含有其中任何一个pfam_id，则不能认定是转录因子，即使1，2都符合。
    :param raw_rules:
    :return:
    """
    rule_txt = raw_rules
    rule_list = rule_txt.strip().split('\n')
    header = rule_list[0].split('\t')
    domain_num_limit = 10
    rule_dict_list = list()
    # Family	SubFamily	DNA-binding domain	Auxiliary domain	Forbidden domain
    for each in rule_list[1:]:
        tmp_dict = dict(zip(header, each.split('\t')))
        # sub_family = tmp_dict['SubFamily']
        binding_domain = tmp_dict["DNA-binding domain"]
        auxiliary_domain = tmp_dict['Auxiliary domain']
        forbidden_domain = tmp_dict['Forbidden domain']
        judge_rule_dict = {
            "binding": dict(),
            "auxiliary": list(),
            "forbidden": list(),
        }
        # 这里针对binding_domain 添加判定规则, 只考虑'>'和'<', 只考虑and关系，不考虑or关系，因为目前不存在这种状况
        match_result = re.findall(r"[^()]+?\s+\(([<>=]*)(\d+)\)\s+\((P[^()]+?)\)", binding_domain)
        if match_result:
            for sign, num, pf_id in match_result:
                if not sign:
                    judge_rule_dict['binding'][pf_id] = [int(num)]
                elif ">" in sign:
                    judge_rule_dict['binding'][pf_id] = range(int(num), domain_num_limit)
                elif "<" in sign:
                    judge_rule_dict['binding'][pf_id] = range(1, int(num)+1)
                else:
                    pass
        match_result = re.findall(r"[^()]+?\s+\((P[^()]+?)\)", binding_domain)
        if match_result:
            for pf_id in match_result:
                judge_rule_dict['binding'][pf_id] = range(1, domain_num_limit)
        #  这里针对auxiliary_domain, 只考虑or关系， 不考虑domain数量的要求情况，因为目前其他情况不存在
        match_result = re.findall(r"[^()]+?\s+\((P[^()]+?)\)", auxiliary_domain)
        if match_result:
            for pf_id in match_result:
                judge_rule_dict['auxiliary'].append(pf_id)
        #  这里针对forbidden_domain, 只考虑or关系， 不考虑domain数量的要求情况，因为目前其他情况不存在
        match_result = re.findall(r"[^()]+?\s+\((P[^()]+?)\)", forbidden_domain)
        if match_result:
            for pf_id in match_result:
                judge_rule_dict['forbidden'].append(pf_id)
        # save rule dict
        tmp_dict['judge_rule'] = judge_rule_dict
        rule_dict_list.append(tmp_dict)
    else:
        return rule_dict_list


def parse_animal_tf_judge_rules(raw_rules):
    rule_list = raw_rules.strip().split("\n")
    header = rule_list[0]
    family2pf_id = {x.split('\t')[2]: x.split('\t')[0] for x in rule_list}
    return family2pf_id



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
        A simple Wrapper for hmmscan.
        hmmscan is used to search protein sequences against collections of pro-
        tein  profiles. For each sequence in <seqfile>, use that query sequence
        to search the target database of profiles in <hmmdb>, and output ranked
        lists  of  the  profiles  with  the  most  significant  matches  to the
        sequence.""", )
    parser.add_argument("-hmmdb", required=True, help="""
        The <hmmdb> needs to be  press'ed  using  hmmpress  before  it  can  be
        searched  with  hmmscan.   This  creates  four  binary  files, suffixed
        .h3{fimp}.
        """)
    parser.add_argument("-seqfile", required=True, help="""
         The <seqfile> may contain more than one query sequence. It  can  be  in
           FASTA  format,  or several other common sequence file formats (genbank,
           embl, and uniprot, among others), or in alignment file formats  (stock-
           holm,  aligned  fasta, and others). See the --qformat option for a com-
           plete list.
        """)
    parser.add_argument("-o", default="raw_result.txt", help="""
        Direct the main human-readable output to a file <f>  instead  of
        the default stdout.
        """)
    parser.add_argument("-tblout", default="tblout.txt", help="""
        Save  a  simple  tabular  (space-delimited) file summarizing the
        per-target output, with one  data  line  per  homologous  target
        model found.
        """)
    parser.add_argument("-domtblout", default="domtblout.tmp.txt", help="""
        Save  a  simple  tabular  (space-delimited) file summarizing the
        per-domain output, with one  data  line  per  homologous  domain
        detected in a query sequence for each homologous model.
        """)
    parser.add_argument("-E", default=10, type=float, help="""
        In the per-target output, report target profiles with an E-value
        of  <= <x>.  The default is 10.0, meaning that on average, about
        10 false positives will be reported per query, so  you  can  see
        the  top  of  the  noise  and decide for yourself if it's really
        noise.
        """)
    parser.add_argument("-domE", default=10, type=float, help="""
        In  the per-domain output, for target profiles that have already
        satisfied the per-profile reporting threshold, report individual
        domains  with  a  conditional E-value of <= <x>.  The default is
        10.0.  A conditional E-value means the expected number of  addi-
        tional  false  positive  domains  in the smaller search space of
        those comparisons that already satisfied the per-profile report-
        ing threshold (and thus must have at least one homologous domain
        already).
        """)
    parser.add_argument("-cpu", default=12, type=int, help="""
        Set  the  number of parallel worker threads to <n>.  By default,
        HMMER sets this to the number of CPU cores it  detects  in  your
        machine  -  that is, it tries to maximize the use of your avail-
        able processor cores. Setting <n>  higher  than  the  number  of
        available  cores  is of little if any value, but you may want to
        set it to something less. You can also control  this  number  by
        setting an environment variable, HMMER_NCPU.
          This  option  is only available if HMMER was compiled with POSIX
        threads support. This is the  default,  but  it  may  have  been
        turned off for your site or machine for some reason.
        """)
    parser.add_argument("-other_args", nargs="*", default='', help="""
        All other optional arguments of hmmscan are supported. But, you have to quote""")

    args = parser.parse_args()













