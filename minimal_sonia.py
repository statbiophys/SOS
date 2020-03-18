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

import os
import numpy as np

class MinimalSonia(object):

    def __init__(self, feature_file,Z=1):
        self.feature_file=feature_file
        self.add_features()
        self.Z=Z
        
    def add_features(self):
        with open(self.feature_file, 'r') as features_file:
            all_lines = features_file.read().strip().split('\n')[1:] #skip header
            features = np.array([l.split('\t')[0].split(';') for l in all_lines])
            feature_energies = np.array([float(l.split('\t')[-1]) for l in all_lines]).reshape((len(features), 1))
        self.features = features
        self.feature_dict = {tuple(f): i for i, f in enumerate(self.features)}
        self.feature_energies=feature_energies

    def find_seq_features(self, seq):
        seq_feature_lsts = [['l' + str(len(seq[0]))]]
        seq_feature_lsts += [['a' + aa + str(i)] for i, aa in enumerate(seq[0])]
        seq_feature_lsts += [['a' + aa + str(-1-i)] for i, aa in enumerate(seq[0][::-1])]
        v_genes = [gene for gene in seq[1:] if 'v' in gene.lower()]
        j_genes = [gene for gene in seq[1:] if 'j' in gene.lower()]
        #Allow for just the gene family match
        v_genes += [gene.split('-')[0] for gene in seq[1:] if 'v' in gene.lower()]
        j_genes += [gene.split('-')[0] for gene in seq[1:] if 'j' in gene.lower()]
        try:
            seq_feature_lsts += [['v' + '-'.join([str(int(y)) for y in gene.lower().split('v')[-1].split('-')])] for gene in v_genes]
            seq_feature_lsts += [['j' + '-'.join([str(int(y)) for y in gene.lower().split('j')[-1].split('-')])] for gene in j_genes]
            seq_feature_lsts += [['v' + '-'.join([str(int(y)) for y in v_gene.lower().split('v')[-1].split('-')]), 'j' + '-'.join([str(int(y)) for y in j_gene.lower().split('j')[-1].split('-')])] for v_gene in v_genes for j_gene in j_genes]
        except ValueError:
            pass
        seq_features = list(set([self.feature_dict[tuple(f)] for f in seq_feature_lsts if tuple(f) in self.feature_dict]))
        return seq_features

    def compute_seq_energy_from_parameters(self,seqs = None):
        try:
            if isinstance(seqs[0], str):
                seqs = [seqs]
        except:
            return None
        seqs_features = [self.find_seq_features(seq) for seq in seqs]
        return np.array([np.sum(self.feature_energies[seq_features]) for seq_features in seqs_features])
    
    def compute_sel_factor(self,seqs = None):
        energies=self.compute_seq_energy_from_parameters(seqs)
        return np.exp(-energies)/self.Z