import pickle
from Bio import SeqIO
from tqdm import tqdm
import re

# gff_path = r'/mnt/d/tinyHMM/data/GCF_000001405.40/chr1.gff'
# genome_path = r'/mnt/d/tinyHMM/data/GCF_000001405.40/chr1.fna'
gff_path = r'/mnt/d/tinyHMM/data/GCF_000001405.40/genomic.gff'
genome_path = r'/mnt/d/tinyHMM/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna'

def pickle_genome(genome_path, pickled_genome_path):
    seq_dic = SeqIO.to_dict(SeqIO.parse(genome_path, 'fasta'))
    seq_dic = dict(sorted(seq_dic.items(), key = lambda item : len(item[1].seq), reverse = True))
    for id, SeqRecord in seq_dic.items():
        seq_dic[id] = str(SeqRecord.seq)
    pickle.dump(seq_dic, open(pickled_genome_path, 'wb'))
    
def pickle_gff(gff_path, pickled_gff_path):
    with open(gff_path, 'r') as f:
        gff_lis = f.readlines()
    pickle.dump(gff_lis, open(pickled_gff_path, 'wb'))



# seq_dic = pickle.load(open('/mnt/d/tinyHMM/data/pickle/seq_dic.pkl', 'rb'))
# gff_lis = pickle.load(open('/mnt/d/tinyHMM/data/pickle/gff_lis.pkl', 'rb'))    

# A stat is a dict with key as the column name and value as the number of times it appears
src_stat = {}
feature_stat = {}
attribute_keys_stat = {}



def get_stat(line, index, stat_dic):
    if line[index] in stat_dic:
        stat_dic[line[index]] += 1
    else:
        stat_dic[line[index]] = 1   

def get_attribute_stat():
    pass
        
        
gene_stat = {} # need to rename

# TODO: get exon number of genes
def get_gene_structure(gff_lis):
    for line in gff_lis:
        if line.startswith('#'):
            continue
        elif line[2] == 'gene':
            gene_stat[(line[0], line[3], line[4])] = 0
            while line[2] != 'gene':
                if line[2] == 'exon':
                    gene_stat[(line[0], line[3], line[4])] += 1

def resolve_attributes(attributes, only_key = False):
    attributes = attributes.split(';')
    attributes = [i.split('=') for i in attributes]
    attributes = dict(attributes)
    
    if only_key:
        return attributes.keys()
    return attributes

def get_interval(line):
    return (int(line[3]), int(line[4])) 

# for line in gff_lis:
#     if line.startswith('#'):
#         continue
#     else:
#         line = line.strip().split('\t')
#         get_stat(line, 1, src_stat)
#         get_stat(line, 2, feature_stat)
#         for attr in resolve_attributes(line[8], only_key=True):
#             if attr in attribute_keys_stat:
#                 attribute_keys_stat[attr] += 1
#             else:
#                 attribute_keys_stat[attr] = 1
# src_stat = sorted(src_stat.items(), key=lambda x: x[1], reverse=True)     
# feature_stat = sorted(feature_stat.items(), key=lambda x: x[1], reverse=True)
    
# print(len(gff_lis))
# print(src_stat)
# print(feature_stat)

# print(resolve_attributes(gff_lis[1000].split('\t')[8], only_key=True))
# print(attribute_keys_stat)

def generate_hmmer_model_from_pickled_gff():
    pass

def generate_hmm_model_from_pickled_gff(pickled_gff, pickled_genome, hmm_model_path):
    gff_lis = pickle.load(open(pickled_gff, 'rb'))
    genome_dic = pickle.load(open(pickled_genome, 'rb'))
    
    coding_dic = {}
    contig = str()
    for line in gff_lis:
        if line.startswith('##sequence-region'):
            line = line.strip().split(' ')
            contig = line[1]
            coding_dic[contig] = []
        elif line.startswith('#'):
            continue
        else:
            line = line.strip().split('\t')
            feature = line[2]
            if feature == 'gene':
                coding_dic[contig].append(get_interval(line))
    del(gff_lis)     
                
    model_dic = {}
    for contig, seq in genome_dic.items():
        model_dic[contig] = 'N' * len(seq)
       
       
    for contig, intervals in coding_dic.items():
        if intervals != []:
            intervals = sorted(intervals, key=lambda x: x[0])
            bar_intervals = tqdm(intervals)
            for interval in bar_intervals:  
                bar_intervals.set_description(f'Processing {contig}')
                # model_dic[contig] = model_dic[contig][:interval[0]-1] + 'C' * (interval[1] - interval[0] + 1) + model_dic[contig][interval[1]:]
                model_dic[contig] = ''.join([model_dic[contig][:interval[0]-1], 'C' * (interval[1] - interval[0] + 1), model_dic[contig][interval[1]:]])
    
    pickle.dump(model_dic, open(hmm_model_path, 'wb'))
     
    # with open(r'./model_dic_test', 'w') as out:
    #     for contig, coding in model_dic.items():
    #         out.writelines(contig + ': ' + coding)
    
    # with open(r'./model_dic_test_read', 'w') as out:
    #     for contig, coding in model_dic.items():
    #         coding_fna = '\n'.join([coding[i:i+50] for i in range(0, len(coding), 50)])
    #         out.writelines(contig + ': \n' + coding_fna)

                

                
