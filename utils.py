'''
    Copyright (C) 2020 Isacchini Giulio
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.'''

import olga.load_model as olga_load_model
import olga.generation_probability as pgen
import olga.sequence_generation as seq_gen
from minimal_sonia import MinimalSonia 
import os
import numpy as np

options_of=['human_T_beta','human_T_alpha','human_B_heavy','mouse_T_beta']
names=['Human T Beta','Human T Alpha','Human B Heavy','Mouse T Beta']
norms=[[0.24566713516135608, 0.9777790833333333],[0.2877415063096418, 0.49879584375], [0.15107851669614455, 0.9770069166666666],[0.27744730165886944 ,0.911731328125]]
qfiles=['features_'+i+'.tsv' for i in options_of]

tabs_styles = {}
tab_style = {}
tab_selected_style = {'backgroundColor': '#ffffff'}

def return_genes(index):
    main_folder=os.path.join(os.path.dirname(olga_load_model.__file__), 'default_models', options_of[index])
    params_file_name = os.path.join(main_folder,'model_params.txt')
    marginals_file_name = os.path.join(main_folder,'model_marginals.txt')
    V_anchor_pos_file = os.path.join(main_folder,'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = os.path.join(main_folder,'J_gene_CDR3_anchors.csv')

    if options_of[index]!='human_T_alpha':
        genomic_data = olga_load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
    else:
        genomic_data = olga_load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
    gene_v=np.unique([value[0].split('*')[0] for value in genomic_data.genV])
    gene_j=np.unique([value[0].split('*')[0] for value in genomic_data.genJ])
    return gene_v,gene_j

def compute_pgen(index,seq):
    index_=int(index)
    main_folder=os.path.join(os.path.dirname(olga_load_model.__file__), 'default_models', options_of[index_])
    params_file_name = os.path.join(main_folder,'model_params.txt')
    marginals_file_name = os.path.join(main_folder,'model_marginals.txt')
    V_anchor_pos_file = os.path.join(main_folder,'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = os.path.join(main_folder,'J_gene_CDR3_anchors.csv')

    if options_of[index_]!='human_T_alpha':
        genomic_data = olga_load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = olga_load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)

    else:
        genomic_data = olga_load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = olga_load_model.GenerativeModelVJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = pgen.GenerationProbabilityVJ(generative_model, genomic_data)

    return pgen_model.compute_aa_CDR3_pgen(seq[0],seq[1],seq[2])/norms[index_][0]

def compute_pgen_for_ppost(index,seqs):
    index_=int(index)
    main_folder=os.path.join(os.path.dirname(olga_load_model.__file__), 'default_models', options_of[index_])
    params_file_name = os.path.join(main_folder,'model_params.txt')
    marginals_file_name = os.path.join(main_folder,'model_marginals.txt')
    V_anchor_pos_file = os.path.join(main_folder,'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = os.path.join(main_folder,'J_gene_CDR3_anchors.csv')

    if options_of[index_]!='human_T_alpha':
        genomic_data = olga_load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = olga_load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)

    else:
        genomic_data = olga_load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = olga_load_model.GenerativeModelVJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = pgen.GenerationProbabilityVJ(generative_model, genomic_data)

    try:
        if isinstance(seqs[0], str):
            seqs = [seqs]
    except:
        return None

    pgens=[pgen_model.compute_aa_CDR3_pgen(seq[0],seq[1],seq[2])/norms[index_][0] for seq in seqs]

    return np.array(pgens)


def compute_pgen_nt(index,cdr3, v_gene, j_gene,type_cdr3):
    index_=int(index)
    main_folder=os.path.join(os.path.dirname(olga_load_model.__file__), 'default_models', options_of[index_])
    params_file_name = os.path.join(main_folder,'model_params.txt')
    marginals_file_name = os.path.join(main_folder,'model_marginals.txt')
    V_anchor_pos_file = os.path.join(main_folder,'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = os.path.join(main_folder,'J_gene_CDR3_anchors.csv')

    if options_of[index_]!='human_T_alpha':
        genomic_data = olga_load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = olga_load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)

    else:
        genomic_data = olga_load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = olga_load_model.GenerativeModelVJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = pgen.GenerationProbabilityVJ(generative_model, genomic_data)
    cdr3=str(cdr3).upper()
    if type_cdr3<1: return pgen_model.compute_nt_CDR3_pgen(cdr3,v_gene,j_gene)/norms[index_][0]
    else:    return 1e-80

def sample_olga(num_gen_seqs=1,chain_index=0, ppost=False,seed=None):
    if seed is not None: np.random.seed(seed)
    else: np.random.seed()

    num_gen_seqs=np.min([num_gen_seqs,1000])
    chain_type=options_of[chain_index]
    main_folder = os.path.join(os.path.dirname(olga_load_model.__file__), 'default_models',chain_type)
    params_file_name = os.path.join(main_folder,'model_params.txt')
    marginals_file_name = os.path.join(main_folder,'model_marginals.txt')
    V_anchor_pos_file = os.path.join(main_folder,'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = os.path.join(main_folder,'J_gene_CDR3_anchors.csv')

    if chain_type=='human_T_alpha':
        genomic_data = olga_load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = olga_load_model.GenerativeModelVJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        sg_model = seq_gen.SequenceGenerationVJ(generative_model, genomic_data)
    else:
        genomic_data = olga_load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = olga_load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        sg_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)

    if not bool(ppost): 
        return [[seq[0],seq[1], genomic_data.genV[seq[2]][0].split('*')[0], genomic_data.genJ[seq[3]][0].split('*')[0]] for seq in [sg_model.gen_rnd_prod_CDR3() for _ in range(int(num_gen_seqs))]] 
    else:
        qm=MinimalSonia(qfiles[chain_index],norms[chain_index][1])
        seqs_post=[['a','b','c']] # initialize
        while len(seqs_post)<num_gen_seqs:
            seqs=[[seq[1], genomic_data.genV[seq[2]][0].split('*')[0], genomic_data.genJ[seq[3]][0].split('*')[0]] for seq in [sg_model.gen_rnd_prod_CDR3() for _ in range(int(11*num_gen_seqs))]] 
            Qs = qm.compute_sel_factor(seqs)
            random_samples=np.random.uniform(size=len(Qs)) # sample from uniform distribution
            #do rejection
            rejection_selection=random_samples < np.clip(Qs,0,10)/10.
            print (np.sum(rejection_selection)/float(len(rejection_selection)))
            seqs_post=np.concatenate([seqs_post,np.array(seqs)[rejection_selection]])
        return seqs_post[1:num_gen_seqs+1]

def nt2aa(ntseq):
    nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    aa_dict ='KQE*TPASRRG*ILVLNHDYTPASSRGCILVFKQE*TPASRRGWMLVLNHDYTPASSRGCILVF'
    
    return ''.join([aa_dict[nt2num[ntseq[i]] + 4*nt2num[ntseq[i+1]] + 16*nt2num[ntseq[i+2]]] for i in range(0, len(ntseq), 3) if i+2 < len(ntseq)])
