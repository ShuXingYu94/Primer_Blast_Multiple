import os
from Bio import SeqIO
from Bio import Entrez
import primer3
from pandas import DataFrame
from pandas import read_csv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML


def readFASTA(FASTA):
    Index = []
    Seq = []
    n = 0
    for i in FASTA:
        if '>' in i and n == 0:
            txt = ''
            Index.append(i.replace('>', '').replace('\n', ''))  # replace('\n','')
            n += 1
        elif '>' in i and n > 0:
            Seq.append(txt)
            txt = ''
            Index.append(i.replace('>', '').replace('\n', ''))
        else:
            txt += i.strip()
    Seq.append(txt)
    return Seq, Index


def extract_pairs(out_address):
    primer3_result_df = read_csv(out_address, header=0, index_col=0)
    left_primer = list(primer3_result_df.loc['PRIMER_LEFT_SEQUENCE'].values)
    right_primer = list(primer3_result_df.loc['PRIMER_RIGHT_SEQUENCE'].values)
    primer_list = []
    while left_primer:
        f = left_primer.pop(0)
        r = right_primer.pop(0)
        primer_list.append([f, r])
    return primer_list


def get_target(email: str, gene_id: dict, gene_name=''):
    Entrez.email = email
    if not gene_name:
        filename = gene_id['id']
    else:
        filename = gene_name
    if not os.path.isfile(filename):
        net_handle = Entrez.efetch(db=gene_id['db'], id=gene_id['id'], rettype=gene_id['rettype'],
                                   retmode=gene_id['retmode'])
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print(filename + ' Data Saved.')
    else:
        print('Data found in given local address.')
    print('Parsing...')
    record = SeqIO.read(filename, 'gb')
    seq = str(record.seq)
    return record, seq


def design_primer(seq_args, global_args, out_address):
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
    print('Succeeded. A total of ', primer3_result['PRIMER_PAIR_NUM_RETURNED'], 'primer pair(s) are designed.',
          '\nFor forward primer: ', primer3_result['PRIMER_LEFT_EXPLAIN'], '\nFor reverse primer: ',
          primer3_result['PRIMER_RIGHT_EXPLAIN'], '\nFor primer pairs: ', primer3_result['PRIMER_PAIR_EXPLAIN'])
    primer3_result_table_dict = {}
    for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        primer_id = str(i)
        for key in primer3_result:
            if '_' + primer_id + '_' in key:
                info_tag = key.replace("_" + primer_id, "")
                try:
                    primer3_result_table_dict[info_tag]
                except:
                    primer3_result_table_dict[info_tag] = []
                finally:
                    primer3_result_table_dict[info_tag].append(primer3_result[key])
    index = []
    for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        index.append("PRIMER_PAIR_" + str(i + 1))
    primer3_result_df = DataFrame(primer3_result_table_dict, index=index)
    primer3_result_df = primer3_result_df.T
    primer3_result_df.to_csv(out_address)
    print('Primer pair(s) Designed: \n', primer3_result_df)
    return primer3_result_df


def extract_pairs(out_address):
    primer3_result_df = read_csv(out_address, header=0, index_col=0)
    left_primer = list(primer3_result_df.loc['PRIMER_LEFT_SEQUENCE'].values)
    right_primer = list(primer3_result_df.loc['PRIMER_RIGHT_SEQUENCE'].values)
    primer_list = []
    while left_primer:
        f = left_primer.pop(0)
        r = right_primer.pop(0)
        primer_list.append([f, r])
    return primer_list


def blastn(query_address: str, db_address: str, out_address1: str = '', evalue=0.001, identity=18, task='blastn',
           dust='yes'):
    blastn_cline = NcbiblastnCommandline(query=query_address, db=db_address, evalue=evalue, outfmt=5, out=out_address1,
                                         task=task, dust=dust)
    stout, stderr = blastn_cline()
    result_handle = open(out_address1)
    blast_record = NCBIXML.read(result_handle)
    e_value_thresh = evalue  # set E_value or other parameter and judge if exist
    identities = identity  # set identity for alignments,for primer design:length of primer-2 is recommended
    count = 0  # count number of blast hits
    name_list = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect <= e_value_thresh and hsp.identities >= identities:
                count += 1
                name_list.append(alignment.title)
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('identity:', hsp.identities)
                print('e value:', hsp.expect)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')
    print(count, ' similar sequence found.')
    return count, name_list
