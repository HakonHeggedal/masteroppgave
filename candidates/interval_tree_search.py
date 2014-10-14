'''
Created on 1. okt. 2014

@author: hakon
'''

from intervaltree.bio import GenomeIntervalTree
import intervaltree
from candidates import structure


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
    candidate_list = []
    # add all intervals to the tree
    for prop in sequence_hits:
#         print prop
        seq_name = prop[0]
        strand_dir = prop[1] # forward: + backward: -
        genome_nr = prop[2].split("|")[3] # which genome (and version)
        genome_offset = int(prop[3]) # offset into the genome, 0-indexed
        dna_sequence = prop[4] # the dna_sequence matching this position.
        sequence_info = [strand_dir, seq_name, dna_sequence]
        
        sequence_tree.addi(genome_nr, genome_offset, genome_offset + len(dna_sequence), sequence_info)
    
    
    for tree in sequence_tree:
        # test all intervals to find candidates
        print tree

        for five_interval in sorted(sequence_tree[tree]):
            
            three_range_begin = five_interval.end
            three_range_end = five_interval.begin + MAX_HAIRPIN_LEN
            outside = five_interval.begin + MAX_HAIRPIN_LEN + 1
            
            #TODO: all, then remove wrong strand, then look for 3's
            
            three_sequences = sequence_tree[tree][three_range_begin:three_range_end] 
            three_sequences = three_sequences - sequence_tree[tree][three_range_begin -1] # not partly same sequence
            three_sequences = [s for s in three_sequences if s.data[0] == five_interval.data[0]] # same strand direction
            print 
            print five_interval
            print "sequences:", len(three_sequences)
            print [c for c in three_sequences]
            
            if len(three_sequences) > 0:
                
                candidate_sequences = sequence_tree[tree][five_interval.begin:three_range_end]

#                 find 5' start and ends
                
                starts = {}
                ends = {}
                get_start = {}
                get_end = {}
                
                best_start = 0
                best_start_pos = 0
                best_end = 0
                best_end_pos = 0

                
                # find best start and end position
                for interval in candidate_sequences:
                    
                    start =  interval.begin - five_interval.begin
                    end = interval.end - five_interval.begin
                    name = interval.data[1]
                    frequency = int(name.split("-")[1])
                    
                    starts[start] = frequency if start not in starts else starts[start] + frequency
                    ends[end] = frequency if end not in ends else ends[end] + frequency
                    get_start[end] = start
                    get_end[start] = end
                    
                    if starts[start] > best_start:
                        best_start = starts[start]
                        best_start_pos = start

                                                
                    if ends[end] > best_end:
                        best_end = ends[end]
                        best_end_pos = end
                
                
                
                print "start positions:\n", starts, "goes to:", get_end
                print "ends positions:\n", ends, "starts from:", get_start
                
                print "best start position:", best_start_pos, best_start, get_end[best_start_pos],
                print "best ends position:", best_end_pos, best_end, get_start[best_end_pos]
                
                second_starts = [(s,val) for (s,val) in starts.iteritems() if s < best_start_pos-5 or s > get_end[best_start_pos] ]
                second_ends = [(s,v) for (s,v) in ends.iteritems() if s > best_end_pos+5 or s < get_start[best_end_pos]]
                
                
                if len(second_starts) == 0 or len(second_ends) == 0:
                    continue 
                
                second_start = max(second_starts, key=lambda (k,v): v)
                second_end = max(second_ends, key=lambda(k,v):v)
                
                print "!!!"
                print "second_starts", second_starts, second_start    
                print "second ends", second_ends, second_end
                
                
#                 TODO: assemble the intervals, full length format 999999999
                begin_5 = 0
                end_5 = 0
                begin_3 = 0
                end_3 = 0
                strand_dir = interval.data[0]
                chromosome = tree
#                 if this is a sub-interval, it is not added 
#                 TODO: same strand direction
                if tree in candidate_tree:
                    if candidate_tree[tree][begin_5:end_3]:
                        continue


                
                mapped_sequences = [x for x in candidate_sequences]
                candidate = structure.Candidate(chromosome, strand_dir, begin_5, end_5, begin_3, end_3, mapped_sequences)

                
                candidate_tree[tree][begin_5:end_3] = candidate
                candidate_list.append(candidate)
               
    return candidate_tree, sequence_tree, candidate_list


def get_nondup_intervals(sequence_hits, seqeunce_tree):
    
    dup_free = []
    for prop in sequence_hits:
        seq_name = prop[0]
        strand_dir = prop[1] # forward: + backward: -
        genome_nr = prop[2].split("|")[3] # which genome (and version)
        genome_offset = int(prop[3]) # offset into the genome, 0-indexed
        dna_sequence = prop[4] # the dna_sequence matching this position.
        sequence_info = [strand_dir, seq_name, dna_sequence]
        interval = intervaltree.Interval(genome_offset, genome_offset + len(dna_sequence), sequence_info)
         
        if genome_nr in seqeunce_tree:
            if interval not in seqeunce_tree[genome_nr]:
                dup_free.append(interval)
            
            
        

def add_seq_to_existing_candidates(seqs, candidate_tree):
    
    for sequence in seqs:
        pass
        
# 
# GENOMENR_TO_FASTAFILE = {}
# GENOMENR_TO_FASTAFILE["NC_000001.10"] = "chr1.fa"
# GENOMENR_TO_FASTAFILE["NC_000002.11"] = "chr2.fa"
# GENOMENR_TO_FASTAFILE["NC_000003.11"] = "chr3.fa"
# GENOMENR_TO_FASTAFILE["NC_000004.11"] = "chr4.fa"
# GENOMENR_TO_FASTAFILE["NC_000005.9"] = "chr5.fa"
# GENOMENR_TO_FASTAFILE["NC_000006.11"] = "chr6.fa"
# GENOMENR_TO_FASTAFILE["NC_000007.13"] = "chr7.fa"
# GENOMENR_TO_FASTAFILE["NC_000008.10"] = "chr8.fa"
# GENOMENR_TO_FASTAFILE["NC_000009.11"] = "chr9.fa"
# GENOMENR_TO_FASTAFILE["NC_000010.10"] = "chr10.fa"
# GENOMENR_TO_FASTAFILE["NC_000011.9"] = "chr11.fa"
# GENOMENR_TO_FASTAFILE["NC_000012.11"] = "chr12.fa"
# GENOMENR_TO_FASTAFILE["NC_000013.10"] = "chr13.fa"
# GENOMENR_TO_FASTAFILE["NC_000014.8"] = "chr14.fa"
# GENOMENR_TO_FASTAFILE["NC_000015.9"] = "chr15.fa"
# GENOMENR_TO_FASTAFILE["NC_000016.9"] = "chr16.fa"
# GENOMENR_TO_FASTAFILE["NC_000017.10"] = "chr17.fa"
# GENOMENR_TO_FASTAFILE["NC_000018.9"] = "chr18.fa"
# GENOMENR_TO_FASTAFILE["NC_000019.9"] = "chr19.fa"
# GENOMENR_TO_FASTAFILE["NC_000020.10"] = "chr20.fa"
# GENOMENR_TO_FASTAFILE["NC_000021.8"] = "chr21.fa"
# GENOMENR_TO_FASTAFILE["NC_000022.10"] = "chr22.fa"
# GENOMENR_TO_FASTAFILE["NC_000023.10"] = "chrX.fa"
# GENOMENR_TO_FASTAFILE["NC_000024.9"] = "chrY.fa"
# GENOMENR_TO_FASTAFILE["NC_001807.4"] = "chrM.fa"
# 




# 
# def find_candidates(sequence_hits):
#     ''' finds microRNA candidates from bowtie data using interval trees
#     
#         sequence_hits -- an iterable of lists on bowtie output format:
#   0          1    2                                 3           4                           5                           6
# ['1-15830', '-', 'gi|224589818|ref|NC_000006.11|', '72113295', 'AGCTTCCAGTCGAGGATGTTTACA', 'IIIIIIIIIIIIIIIIIIIIIIII', '0']
#         returns a list of candidates, and the interval tree with all sequences
#     '''
#     
#     sequence_tree = GenomeIntervalTree()
#     candidate_tree = GenomeIntervalTree() # only candidates here
#     # add all intervals to the tree
#     for prop in sequence_hits:
# #         print prop
#         seq_name = prop[0]
#         strand_dir = prop[1] # forward: + backward: -
#         genome_nr = prop[2].split("|")[3] # which genome (and version)
#         genome_offset = int(prop[3]) # offset into the genome, 0-indexed
#         dna_sequence = prop[4] # the dna_sequence matching this position.
#         sequence_info = [strand_dir, seq_name, dna_sequence]
#         
#         sequence_tree.addi(genome_nr, genome_offset, genome_offset + len(dna_sequence), sequence_info)
#     
#     
#     for tree in sequence_tree:
#         # test all intervals to find candidates
# 
#         for five_interval in sorted(sequence_tree[tree]):
#             
#             three_range_begin = five_interval.ends
#             three_range_end = five_interval.begin + MAX_HAIRPIN_LEN
#             outside = five_interval.begin + MAX_HAIRPIN_LEN + 1
#             
#             three_sequences = sequence_tree[tree][three_range_begin:three_range_end]
#         
#             if three_sequences != None:
#                 
#                 three_interval = None
#                 
#                 # find the legal 3' interval longest distance from 5'
#                 for seq in three_sequences:
#                     if seq.ends < outside and seq.begin > five_interval.ends:
#                         if three_interval:
#                             if seq.ends > three_interval.ends:
#                                 three_interval = seq
#                         else:
#                             three_interval = seq
#                         
#                 
# #                 legal_three_ends = [x for x in three_sequences if x.ends < outside and x.data[0] == five_interval.data[0]]
# #                 three_interval = max(legal_three_ends, key=lambda x:x.ends) # 
#                 if not three_interval:
#                     continue
# 
# #                 if this is a sub-interval, it is not added
#                 if tree in candidate_tree:
#                     if candidate_tree[tree][five_interval.begin:three_interval.ends]:
#                         continue
# 
#                 print three_interval
# #                 this interval is 5' in a new candidate
# #                 [strand, 5'name, 5'sequence, 3'name, 3'sequence]
#                 candidate_data = five_interval.data[:]
#                 candidate_data.extend(three_interval.data[1:])
#                 
#                 candidate_tree[tree][five_interval.begin:three_interval.ends] = candidate_data
#                
#     return candidate_tree, sequence_tree 

