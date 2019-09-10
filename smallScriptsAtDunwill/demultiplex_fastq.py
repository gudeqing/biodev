import os
import gzip
import pandas as pd


def read_sample_info(in_file):
    table = pd.read_csv(in_file, sep=None, header=0, index_col=0, engine='python')
    table.columns = ['p7index', 'p5index']
    table.index.name = 'sample'
    pair_index_dict = dict()
    for sample in table.index:
        p7index = table.loc[sample, 'p7index']
        p5index = table.loc[sample, 'p5index']
        indexes = tuple(sorted([p7index, p5index]))
        if indexes in pair_index_dict:
            raise Exception(f'indexes of {sample} is same as '
                            f'{pair_index_dict[indexes]}')
        pair_index_dict[indexes] = sample
    return pair_index_dict


def judge(seq_info, pair_index_dict, write_obj_dict):
    # 决定当前read是否是想要的及写到哪一个文件
    first_line = seq_info[0]
    header, detail = first_line.strip().decode().split()
    lane = '{:0>3d}'.format(int(header.split(':')[3]))
    read_num, filtered, ctrl_num, indexes = detail.rsplit(':')
    indexes = tuple(sorted(indexes.split('+')))
    sample = pair_index_dict.get(indexes)
    if sample:
        w_id = sample + f'_R{read_num}'
        if w_id not in write_obj_dict:
            write_obj = gzip.open(f'{sample}_S1_L{lane}_R{read_num}_001.fastq.gz', 'w')
            write_obj_dict[w_id] = write_obj
        else:
            write_obj = write_obj_dict[w_id]
        return write_obj


async def fastq_parser(fastq, pair_index_dict, write_obj_dict, size=3000):
    seq_list = list()
    with gzip.open(fastq, 'rb') as f:
        seq_info = list()
        seq_info.append(f.readline())
        seq_info.append(f.readline())
        seq_info.append(f.readline())
        seq_info.append(f.readline())
        for next_line in f:
            if next_line.startswith(b'@'):
                # 处理上条read
                write_obj = judge(seq_info, pair_index_dict, write_obj_dict)
                if write_obj:
                    seq_info.append(write_obj)
                    seq_list.append(seq_info)
                    if len(seq_list) >= size:
                        yield seq_list
                        seq_list = []
                # init seq info
                seq_info = [next_line]
            else:
                seq_info.append(next_line)
        else:
            # 处理最后一条read
            write_obj = judge(seq_info, pair_index_dict, write_obj_dict)
            if write_obj:
                seq_info.append(write_obj)
                seq_list.append(seq_info)
            yield seq_list


async def write_seq(seq_list):
    print(len(seq_list), flush=True)
    for seq_info in seq_list:
        write_obj = seq_info[-1]
        write_obj.write(seq_info[0])
        write_obj.write(seq_info[1])
        write_obj.write(seq_info[2])
        write_obj.write(seq_info[3])


async def splitter(fastq, sample_info, size=3000):
    pair_index_dict = read_sample_info(sample_info)
    write_obj_dict = dict()
    async for seq_list in fastq_parser(fastq, pair_index_dict, write_obj_dict, size=size):
        await write_seq(seq_list)

    # close writing obj
    for _, obj in write_obj_dict.items():
            obj.close()

def run_splitter(fastq, sample_info, size=3000):
    """
    拆分fastq
    :param fastq: undetermined fastq file
    :param sample_info: 第一行为header，第一列为样本名，第二列为p7index，第三列为p5index
    :param size: 指定一次读多少条read后再输出到文件
    :return:
    """
    # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
    # <read>:<is filtered>:<control number>:<index sequence>
    import asyncio
    asyncio.run(splitter(fastq, sample_info, size=size))


def split_fastq(fastq, sample_info):
    """
    拆分fastq
    :param fastq: undetermined fastq file
    :param sample_info: 第一行为header，第一列为样本名，第二列为p7index，第三列为p5index
    :return:
    """
    # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
    # <read>:<is filtered>:<control number>:<index sequence>
    pair_index_dict = read_sample_info(sample_info)
    # print(pair_index_dict)
    write_obj_dict = dict()
    with gzip.open(fastq, 'rb') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith(b'@'):
                header, detail = line.strip().decode().split()
                lane = '{:0>3d}'.format(int(header.split(':')[3]))
                read_num, filtered, ctrl_num, indexes = detail.split(':')
                indexes = tuple(sorted(indexes.split('+')))
                sample = pair_index_dict.get(indexes)
                if sample:
                    w_id = sample+f'_R{read_num}'
                    if w_id not in write_obj_dict:
                        write_obj = gzip.open(f'{sample}_S1_L{lane}_R{read_num}_001.fastq.gz', 'w')
                        write_obj_dict[w_id] = write_obj
                    else:
                        write_obj = write_obj_dict[w_id]
                    write_obj.write(line)
                    write_obj.write(f.readline())
                    write_obj.write(f.readline())
                    write_obj.write(f.readline())
    # close writing obj
    for _, obj in write_obj_dict.items():
        obj.close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['split_fastq', 'run_splitter'])
