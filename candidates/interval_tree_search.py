'''
Created on 1. okt. 2014

@author: hakon
'''
import intervaltree
from intervaltree.bio import GenomeIntervalTree



GENOMENR_TO_FASTAFILE = {}
GENOMENR_TO_FASTAFILE["NC_000001.10"] = "chr1.fa"
GENOMENR_TO_FASTAFILE["NC_000002.11"] = "chr2.fa"
GENOMENR_TO_FASTAFILE["NC_000003.11"] = "chr3.fa"
GENOMENR_TO_FASTAFILE["NC_000004.11"] = "chr4.fa"
GENOMENR_TO_FASTAFILE["NC_000005.9"] = "chr5.fa"
GENOMENR_TO_FASTAFILE["NC_000006.11"] = "chr6.fa"
GENOMENR_TO_FASTAFILE["NC_000007.13"] = "chr7.fa"
GENOMENR_TO_FASTAFILE["NC_000008.10"] = "chr8.fa"
GENOMENR_TO_FASTAFILE["NC_000009.11"] = "chr9.fa"
GENOMENR_TO_FASTAFILE["NC_000010.10"] = "chr10.fa"
GENOMENR_TO_FASTAFILE["NC_000011.9"] = "chr11.fa"
GENOMENR_TO_FASTAFILE["NC_000012.11"] = "chr12.fa"
GENOMENR_TO_FASTAFILE["NC_000013.10"] = "chr13.fa"
GENOMENR_TO_FASTAFILE["NC_000014.8"] = "chr14.fa"
GENOMENR_TO_FASTAFILE["NC_000015.9"] = "chr15.fa"
GENOMENR_TO_FASTAFILE["NC_000016.9"] = "chr16.fa"
GENOMENR_TO_FASTAFILE["NC_000017.10"] = "chr17.fa"
GENOMENR_TO_FASTAFILE["NC_000018.9"] = "chr18.fa"
GENOMENR_TO_FASTAFILE["NC_000019.9"] = "chr19.fa"
GENOMENR_TO_FASTAFILE["NC_000020.10"] = "chr20.fa"
GENOMENR_TO_FASTAFILE["NC_000021.8"] = "chr21.fa"
GENOMENR_TO_FASTAFILE["NC_000022.10"] = "chr22.fa"
GENOMENR_TO_FASTAFILE["NC_000023.10"] = "chrX.fa"
GENOMENR_TO_FASTAFILE["NC_000024.9"] = "chrY.fa"
GENOMENR_TO_FASTAFILE["NC_001807.4"] = "chrM.fa"

MAX_HAIRPIN_LEN = 80
MIN_HAIRPIN_LEN = 46


def find_candidates(sequence_hits):
    ''' finds microRNA candidates from bowtie data using interval trees
    
        sequence_hits -- an iterable of lists on bowtie output format:
  0          1    2                                 3           4                           5                           6
['1-15830', '-', 'gi|224589818|ref|NC_000006.11|', '72113295', 'AGCTTCCAGTCGAGGATGTTTACA', 'IIIIIIIIIIIIIIIIIIIIIIII', '0']
        returns a list of candidates, and the interval tree with all sequences
    '''
    
    sequence_tree = GenomeIntervalTree()
    candidate_tree = GenomeIntervalTree() # only candidates here
    # add all intervals to the tree
    for prop in sequence_hits:
        print prop
        seq_name = prop[0]
        strand_dir = prop[1] # forward: + backward: -
        genome_nr = prop[2].split("|")[3] # which genome (and version)
        genome_offset = int(prop[3]) # offset into the genome, 0-indexed
        dna_sequence = prop[4] # the dna_sequence matching this position.
        sequence_info = [strand_dir, seq_name, dna_sequence]
        
        sequence_tree.addi(genome_nr, genome_offset, genome_offset + len(dna_sequence), sequence_info)
    

     
    for tree in sequence_tree:
        # test all intervals to find candidates

        for five_interval in sorted(sequence_tree[tree]):
            
            three_range_begin = five_interval.end
            three_range_end = five_interval.begin + MAX_HAIRPIN_LEN
            outside = five_interval.begin + MAX_HAIRPIN_LEN + 1
            
            three_sequences = sequence_tree[tree][three_range_begin:three_range_end]
        
            if three_sequences:
                legal_three_ends = [x for x in three_sequences if x.end < outside and x.data[0] == five_interval.data[0]]
                
                if not legal_three_ends:
                    continue

                three_interval = max(legal_three_ends, key=lambda x:x[1]) # 

#                 if this is a sub-interval, it is not added
                if tree in candidate_tree:
                    if candidate_tree[tree][five_interval.begin:three_interval.end]:
                        continue

#                 this interval is 5' in a new candidate
#                 [strand, 5'name, 5'sequence, 3'name, 3'sequence]
                candidate_data = five_interval.data[:]
                candidate_data.extend(three_interval.data[1:])
                
                candidate_tree[tree][five_interval.begin:three_interval.end] = candidate_data
               
    return candidate_tree, sequence_tree 





