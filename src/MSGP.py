import pickle
import numpy as np
import math
from tqdm import tqdm
from numba import jit
import gff_resolver

# {Token: State}
# anno_dic = {'AGTCTACGCTAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGTGTAGCTACGTGTCAGTTCAGTCGATGCAACATG': 
#             'NNNNNNNNNNNNCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNCNNCNCNNCNCNNCCCNNCNCN'}

def generate_anno_dic(genome_dic, hmm_model):
    anno_dic = {}
    for token, state in hmm_model.items():
        anno_dic[genome_dic[token]] = state
    return anno_dic

@jit
def training(anno_dic):
    trans = {'NN': 0.0, 'NC': 0.0, 'CN': 0.0, 'CC': 0.0}
    emiss = {'nA':0.0, 'nC':0.0, 'nG':0.0, 'nT':0.0, 'cA':0.0, 'cC':0.0, 'cG':0.0, 'cT':0.0}


    for token, state in anno_dic.items():
        print(f'Processing {len(token)/1000000}M sequence')
        # training emission probability
        # bar_range_len_token = tqdm(range(len(token)))
        for i in range(len(token)):
            # bar_range_len_token.set_description(f'Training emission probability: {i}/{len(token)}')
            if token[i] == 'A' and state[i] == 'N':
                emiss['nA'] += 1
            elif token[i] == 'C' and state[i] == 'N':
                emiss['nC'] += 1
            elif token[i] == 'G' and state[i] == 'N':
                emiss['nG'] += 1
            elif token[i] == 'T' and state[i] == 'N':
                emiss['nT'] += 1
            elif token[i] == 'A' and state[i] == 'C':
                emiss['cA'] += 1
            elif token[i] == 'C' and state[i] == 'C':
                emiss['cC'] += 1
            elif token[i] == 'G' and state[i] == 'C':
                emiss['cG'] += 1
            elif token[i] == 'T' and state[i] == 'C':
                emiss['cT'] += 1   

            # if i % 1000000 == 0:
            #     print(f'{i/1000000}/{len(token)/1000000}M')
        
        # training transition probability
        prev = state[0]
        # bar_state = tqdm(state)
        for i in state:
            # bar_state.set_description(f'Training transition probability: {i}/{len(state)}')
            if i == 'N' and prev == 'N': # N -> N
                trans['NN'] += 1
            elif i == 'C' and prev == 'N': # N -> C
                trans['CN'] += 1
            elif i == 'N' and prev == 'C': # C -> N
                trans['NC'] += 1
            elif i == 'C' and prev == 'C': # C -> C
                trans['CC'] += 1
            prev = i
            
            
    (trans['NN'], trans['CN']) = (trans['NN'] / (trans['NN'] + trans['CN']), trans['CN'] / (trans['NN'] + trans['CN']))
    (trans['NC'], trans['CC']) = (trans['NC'] / (trans['NC'] + trans['CC']), trans['CC'] / (trans['NC'] + trans['CC']))
    total_n = emiss['nA'] + emiss['nC'] + emiss['nG'] + emiss['nT']
    total_c = emiss['cA'] + emiss['cC'] + emiss['cG'] + emiss['cT']
    (emiss['nA'], emiss['nC'], emiss['nG'], emiss['nT']) = (emiss['nA'] / total_n, emiss['nC'] / total_n, emiss['nG'] / total_n, emiss['nT'] / total_n)
    (emiss['cA'], emiss['cC'], emiss['cG'], emiss['cT']) = (emiss['cA'] / total_c, emiss['cC'] / total_c, emiss['cG'] / total_c, emiss['cT'] / total_c)

    for key, value in trans.items():
        trans[key] = math.log(value)
    for key, value in emiss.items():
        emiss[key] = math.log(value)
    
    return trans, emiss

# def predicting(token, trans: dict[str, float], emiss: dict[str, float]):
def predicting(token, trans: dict, emiss: dict):
    c_path, n_path = [], []
    pc_prev, pn_prev = trans['NC'], trans['NN']
    for i in token:
        pc = emiss[('c' + i)] + max(pc_prev + trans['CC'], pn_prev + trans['NC'])
        if pc_prev + trans['CC'] > pn_prev + trans['NC']:
            c_path.append('C')
        else:
            c_path.append('N')
        
        pn = emiss[('n' + i)] + max(pc_prev + trans['CN'], pn_prev + trans['NN'])
        if pc_prev + trans['CN'] > pn_prev + trans['NN']:
            n_path.append('C')
        else:
            n_path.append('N')
        
        pc_prev, pn_prev = pc, pn
        print(pc, pn)
        
    return ''.join(c_path), ''.join(n_path)

# genome_dic = pickle.load(open('/mnt/d/TinyHMM/data/pickle/seq_dic.pkl', 'rb'))
# hmm_model = pickle.load(open('/mnt/d/TinyHMM/data/pickle/model_dic.pkl', 'rb'))
# trans, emiss = training(generate_anno_dic(genome_dic, hmm_model))

genome_path = '/mnt/d/tinyHMM/data/MPVR/GCF_014621545.1/GCF_014621545.1_ASM1462154v1_genomic.fna'
gff_path = '/mnt/d/tinyHMM/data/MPVR/GCF_014621545.1/genomic.gff'
pickled_genome_path = '/mnt/d/tinyHMM/data/pickle/MPVR/seq_dic.pkl'
pickled_gff_path = '/mnt/d/tinyHMM/data/pickle/MPVR/gff_dic.pkl'
hmm_model_path = '/mnt/d/tinyHMM/data/pickle/MPVR/model_dic.pkl'
gff_resolver.pickle_genome(genome_path, pickled_genome_path)
gff_resolver.pickle_gff(gff_path, pickled_gff_path)


seq_dic = pickle.load(open(pickled_genome_path, 'rb'))
gff_lis = pickle.load(open(pickled_gff_path, 'rb'))
gff_resolver.generate_hmm_model_from_pickled_gff(pickled_gff_path, pickled_genome_path, hmm_model_path)
hmm_model = pickle.load(open(hmm_model_path, 'rb'))
trans, emiss = training(generate_anno_dic(seq_dic, hmm_model))
print(trans, emiss)


pred_seq = ''
with open('/mnt/d/tinyHMM/data/MPVR/GCF_014621545.1/GCF_014621545.1_ASM1462154v1_genomic.fna', 'r') as f:
    for line in f:
        if line[0] == '>':
            continue
        pred_seq += line.strip()
    pred_seq = pred_seq.replace('\n', '')

print(predicting(pred_seq, trans, emiss))