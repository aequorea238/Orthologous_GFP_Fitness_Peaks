import numpy as np
import collections
import matplotlib.pyplot as plt
from Bio import pairwise2
from tqdm import tqdm

def flatten(l):
    """Flatten a list of lists"""
    return [item for sublist in l for item in sublist]

class Clone:
    """This class describes individual clones, i.e. primary barcodes. Associates their corresponding nucleotide sequence, protein sequence, and list of mutations."""
    def __init__(self, bc, sq):
        self.bc = bc
        self.sq = sq
         
        
    def extract_mutations(self, ref_sq, catch_indels=True, gapper='.'):
        """ Compares the clone's sequence with the reference wild-type sequence. Detects and labels mutations in the format 'AiB', where A is the reference state, B is the mutated state, and i is the position. If 'catch_indels' is enabled (default), employs a BioPython global pairwise alignment of the clone and reference sequences, otherwise, performs a simple and quicker comparison of states at each position which detects substitutions but not indels."""
        self.ntmut = []
        self.ntmut_pos = []
        
        if catch_indels==False:
            for i in range(min(len(self.sq), len(ref_sq))):
                if self.sq[i] != ref_sq[i]:
                    self.ntmut.append(ref_sq[i] + str(i) + self.sq[i])
                    self.ntmut_pos.append(i)
                    
        else:
            sq = self.sq
            alignment = pairwise2.align.globalms(ref_sq, sq, 2,-1,-2,-0,penalize_end_gaps=False, 
                                         one_alignment_only=True, gap_char=gapper)[0]
            aligned_ref, aligned_self = alignment[0], alignment[1]
            self.sq_alignment = {'ref':aligned_ref, 'self':aligned_self}
            
            indel_adjuster = 0
            muts = []
            muts_pos = []

            for i in range(len(aligned_ref)):
                if aligned_ref[i] != aligned_self[i]:
                    if aligned_self[i] == gapper and aligned_self[i-1] != gapper: # deletion
                        deletion_pos = i+indel_adjuster
                        muts.append(aligned_ref[i] + str(i+indel_adjuster) + aligned_self[i])
                        muts_pos.append(i+indel_adjuster)
                    elif aligned_self[i] == gapper and aligned_self[i-1] == gapper: # multi-nt deletion
                        muts.append(aligned_ref[i] + str(deletion_pos) + gapper)
                        muts_pos.append(deletion_pos)
                    else:
                        muts.append(aligned_ref[i] + str(i+indel_adjuster) + aligned_self[i])
                        muts_pos.append(i+indel_adjuster)
                    if aligned_ref[i] == gapper: # insertion, need to adjust index to keep the later positions consistent
                        indel_adjuster -= 1        
                        
            self.ntmut = muts
            self.ntmut_pos = muts_pos
                    
                    
    def group_tandem_indels(self, gapper='.'):
        """ In cases of multi-nucleotide indels, groups all affected nucleotides into one 'mutation' instead of having multiple adjacent single-nucleotide insertions or deletions"""
        muts = list(self.ntmut)
#         muts_pos = self.ntmut_pos_new
        unique_pos = sorted(set(self.ntmut_pos))    
        if len(self.ntmut_pos) != len(unique_pos):
            for i,pos in enumerate(unique_pos):
                if self.ntmut_pos.count(pos)>1:
                    multis = [mut for mut in muts if mut[0:-1]==gapper+str(pos) or mut[1:]==str(pos)+gapper]
                    muts = [mut for mut in muts if mut not in multis]
#                     print (multis)
                    assert set([mut[0] for mut in multis]) == {gapper} or set([mut[-1] for mut in multis]) == {gapper}
                    if set([mut[0] for mut in multis]) == {gapper}: # insertion
                        new_mut = gapper + str(pos) + ''.join([mut[-1] for mut in multis])
                    elif set([mut[-1] for mut in multis]) == {gapper}: # deletion
                        new_mut = ''.join([mut[0] for mut in multis]) + str(pos) + gapper
                    muts.insert(i, new_mut)
        
        self.ntmut_tandemgrouped = muts
        self.ntmut_pos_tandemgrouped = [int(''.join([x for x in mut if x in '1234567890'])) for mut in muts]
        
        
    def translate(self, genetic_code):
        """ Gets protein sequence from DNA sequence"""
        assert self.sq[0:3] == 'ATG'
        protein = ''
        for i in range(0, len(self.sq), 3):
            if i+3 > len(self.sq):
                break
            protein += genetic_code[self.sq[i:i+3]]
        self.protein = protein
        
        
    def extract_aa_subs(self, ref_prot, gapper='.'):
        """ Detects and labels amino acid substitutions and indels by comparing the clone's protein sequence with the reference."""
        self.aasub = []
        self.aasub_pos = []

        n_indels = len([mut for mut in self.ntmut if gapper in mut])
        if n_indels>0 and n_indels%3==0:
            sq = self.protein
            alignment = pairwise2.align.globalms(ref_prot, sq, 2,-1,-2,-0,penalize_end_gaps=False, 
                                         one_alignment_only=True, gap_char=gapper)[0]
            aligned_ref, aligned_self = alignment[0], alignment[1]
            self.prot_alignment = {'ref':aligned_ref, 'self':aligned_self}

            indel_adjuster = 0
            muts = []
            muts_pos = []

            for i in range(len(aligned_ref)):
                if aligned_ref[i] != aligned_self[i]:
                    if aligned_self[i] == gapper and aligned_self[i-1] != gapper: # deletion
                        deletion_pos = i+indel_adjuster
                        muts.append(aligned_ref[i] + str(i+indel_adjuster) + aligned_self[i])
                        muts_pos.append(i+indel_adjuster)
                    elif aligned_self[i] == gapper and aligned_self[i-1] == gapper: # multi-nt deletion
                        muts.append(aligned_ref[i] + str(deletion_pos) + gapper)
                        muts_pos.append(deletion_pos)
                    else:
                        muts.append(aligned_ref[i] + str(i+indel_adjuster) + aligned_self[i])
                        muts_pos.append(i+indel_adjuster)
                    if aligned_ref[i] == gapper: # insertion, need to adjust index to keep the later positions consistent
                        indel_adjuster -= 1        

            self.aasub = muts
            self.aasub_pos = muts_pos

        else:
            for i in range(min(len(self.protein), len(ref_prot))):
                if self.protein[i] != ref_prot[i]:
                    self.aasub.append(ref_prot[i] + str(i) + self.protein[i])
                    self.aasub_pos.append(i)          
                
                
class Library:
    """ This class describes whole mutant libraries, populated by Clones. Contains library-wide statistics such as number of unique genotypes, average number of mutations per clone, etc."""
    def __init__(self, fasta):
        self.clones = []
        with open(fasta, 'r') as f:
            lines = f.readlines()
        for i in range(0, len(lines), 2):
            self.clones.append(Clone(bc = lines[i][1:-1], sq = lines[i+1][:-1]))
            
        print ('There are %s unique barcodes and %s unique mutants.' %(len(self.clones), len(set([clone.sq for clone in self.clones]))))
            
            
    def extract_mutations(self, ref_sq, catch_indels=True, progress_bar=False):
        if progress_bar:
            for clone in tqdm(self.clones):
                clone.extract_mutations(ref_sq, catch_indels)
        else:
            for clone in self.clones:
                clone.extract_mutations(ref_sq, catch_indels)
        
#         max_count = len(self.clones)
#         f = IntProgress(min=0, max=max_count) # instantiate the bar
#         display(f) # display the bar
#         for clone in self.clones:
#             clone.extract_mutations(ref_sq)
#             f.value += 1 # signal to increment the progress bar
            
        unique_ntmut = set(flatten([sq.ntmut for sq in self.clones]))
        print('This library contains %d unique nucleotide mutations.' % len(unique_ntmut))
        if catch_indels==True:
            print('Indels have been properly labeled.')
        
        
    def group_tandem_indels(self,  gapper='.'):
        for clone in  self.clones:
            clone.group_tandem_indels( gapper)
        print('Tandem indels have been grouped together.')
        
     
    def get_mut_numbers_per_clone(self):
        muts = []
        for clone in self.clones:
            muts.append(len(clone.ntmut))
        print('There are %s mutations per clone on average (median: %s).' % (np.mean(muts), np.median(muts)))

        
        
    def translate(self, genetic_code):
        for clone in self.clones:
            clone.translate(genetic_code)
        print('There are %s unique protein sequences.' % len(set([clone.protein for clone in self.clones])))
        
        
    def extract_aa_subs(self, ref_prot, gapper='.', progress_bar=False):
        if progress_bar:
            for clone in tqdm(self.clones):
                clone.extract_aa_subs(ref_prot)
        else:
            for clone in self.clones:
                clone.extract_aa_subs(ref_prot)
                
        unique_aasub = set(flatten([sq.aasub for sq in self.clones]))
        print('This library contains %d unique amino acid changes' % len(unique_aasub))
    
    
    def mut_counts(self):
        self.ntmut_counts = {}
        self.aasub_counts = {}
        
        for clone in self.clones:
            for mut in clone.ntmut:
                if mut in self.ntmut_counts:
                    self.ntmut_counts[mut] += 1
                else:
                    self.ntmut_counts[mut] = 1
                    
            for mut in clone.aasub:
                if mut in self.aasub_counts:
                    self.aasub_counts[mut] += 1
                else:
                    self.aasub_counts[mut] = 1
        print('The number of times each unique mutation appears in the library has been quantified.')
        
        
    def get_aasub_numbers_per_clone(self):
        aamuts = []
        for clone in self.clones:
            aamuts.append(len(clone.aasub))
        print('There are %s amino acid substitutions per clone on average (median: %s).' % (np.mean(aamuts), np.median(aamuts)))
        