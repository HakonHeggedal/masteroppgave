'''
Created on 1. okt. 2014

@author: hakon
'''

from intervaltree.bio import GenomeIntervalTree
import intervaltree
from candidates import structure
from gi.overrides.keysyms import percent


MAX_HAIRPIN_LEN = 80
MIN_HAIRPIN_LEN = 46

MIN_HAIRPIN_LOOP = 10
# MIN_PRIME_SEQ = 10

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
    all_mapped_sequences = set()
    seq_to_candidates = {}
    
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
    
    # test all intervals to find candidates
    for tree in sequence_tree:
        print tree

        for five_interval in sorted(sequence_tree[tree]):
            
            if five_interval in all_mapped_sequences:
                continue
            
            three_range_begin = five_interval.end
            three_range_end = five_interval.begin + MAX_HAIRPIN_LEN
            

            
            #test if this interval is present in another candidate

            # all candidates in allowed range
            candidate_sequences = sequence_tree[tree][five_interval.begin:three_range_end]
            
            if len(candidate_sequences) <= 1:
                continue

            outside_after = sequence_tree[tree][five_interval.begin+MAX_HAIRPIN_LEN+1]
            outside_before = sequence_tree[tree][five_interval.begin-1]           
            
            
            candidate_sequences.difference_update(outside_before, outside_after)
            candidate_sequences = [s for s in candidate_sequences if s.data[0] == five_interval.data[0]]
            
#             three_sequences = sequence_tree[tree][three_range_begin:three_range_end]
#             print sequence_tree[tree][five_interval.begin:three_range_begin -1]

            three_sequences = set(candidate_sequences).difference(sequence_tree[tree][five_interval.begin:three_range_begin -1])

#             three_sequences = [s for s in three_sequences if s.data[0] == five_interval.data[0]] # same strand direction

            
            if len(three_sequences) > 0:

#                 find 5' start and ends
                
                starts = {} # start positions -> frequency
                ends = {} # end positions -> frequency
                get_start = {} # end position -> first start position
                get_end = {} # start positon -> last end position
                
                best_start = 0
                best_start_pos = 0
                best_end = 0
                best_end_pos = 0

                
                # find best start and end position
                for interval in candidate_sequences:
                    
                    start =  interval.begin - five_interval.begin # start position
                    if start < 0:
                        continue
                    end = interval.end - five_interval.begin # end position
                    name = interval.data[1]
                    frequency = float(name.split("-")[1])
                    
                    starts[start] = frequency if start not in starts else starts[start] + frequency
                    ends[end] = frequency if end not in ends else ends[end] + frequency
                    
                    # find the largest interval from given start position
                    get_start[end] = start if end not in get_start else min(start, get_start[end])
                    get_end[start] = end if start not in get_end else max(end, get_end[start])
                    
                    if starts[start] > best_start:
                        best_start = starts[start]
                        best_start_pos = start
                                                
                    if ends[end] > best_end:
                        best_end = ends[end]
                        best_end_pos = end
                
                
#                 print "start positions:\n", starts, "goes to:", get_end
#                 print "ends positions:\n", ends, "starts from:", get_start
#                 
#                 print "best start position:", best_start_pos, best_start, get_end[best_start_pos],
#                 print "best ends position:", best_end_pos, best_end, get_start[best_end_pos]
                
                #TODO not overlapping or very close
                second_starts = set()
                second_ends = set()
                
                for (start, count) in starts.iteritems():
                    if start < best_start_pos: # this is 5'
                        if get_end[start] + MIN_HAIRPIN_LOOP < best_start_pos:
                            second_starts.add( (start, count) )
                    elif get_end[best_start_pos] + MIN_HAIRPIN_LOOP < start:
                        second_starts.add( (start, count) )
                
                for (end, count) in ends.iteritems():
                    if end < best_end_pos: # this is 5'
                        if end + MIN_HAIRPIN_LOOP < get_start[best_end_pos]: 
                            second_ends.add( (end, count) )
                    elif best_end_pos + MIN_HAIRPIN_LOOP < get_start[end]: 
                        second_ends.add( (end, count) )
                    
                
#                 second_starts = [(s,val) for (s,val) in starts.iteritems() if s < best_start_pos-5 or s > get_end[best_start_pos] ]
#                 second_ends = [(s,v) for (s,v) in ends.iteritems() if s > best_end_pos+5 or s < get_start[best_end_pos]]
                
                
                if len(second_starts) == 0 or len(second_ends) == 0:
                    continue 
                
                second_start = max(second_starts, key=lambda (k,v): v)
                second_end = max(second_ends, key=lambda(k,v):v)
                
#                 print "!!!"
#                 print "second_starts", second_starts, second_start    
#                 print "second ends", second_ends, second_end

                begin_5 = five_interval.begin + min(best_start_pos, second_start[0])
                end_5 = five_interval.begin + min(best_end_pos, second_end[0])
                begin_3 = five_interval.begin + max(best_start_pos, second_start[0])
                end_3 = five_interval.begin + max(best_end_pos, second_end[0])
                strand_dir = interval.data[0]
                chromosome = tree
                

                if tree in candidate_tree:
                    if candidate_tree[tree][begin_5:end_3]:
                        continue

                candidate = structure.Candidate(chromosome,
                                                 strand_dir,
                                                 begin_5,
                                                 end_5,
                                                 begin_3,
                                                 end_3,
                                                 candidate_sequences)
                
                for candidate_interval in candidate_sequences:
                    name = candidate_interval.data[1]
                    if name not in seq_to_candidates:
                        number_id = int(name.split("-")[0])
                        duplicates = float(name.split("-")[1])
                        s = structure.Sequence(number_id, duplicates, candidate_interval.data[2])
                        s.add_candidate(candidate)
                        seq_to_candidates[name] = s
                    else:
                        seq_to_candidates[name].add_candidate(candidate)
                    
                candidate_tree[tree][begin_5:end_3] = candidate
                candidate_list.append(candidate)

                if len(all_mapped_sequences) == 0:
                    all_mapped_sequences = candidate.all_mapped_sequences
               
    return candidate_tree, sequence_tree, candidate_list, seq_to_candidates





def align_miRNAs(mirna_hits, hairpinID_to_mature, candidate_tree, sequence_tree):
    print "align candidates with miRNA"
    candidate_to_miRNAid = {}
    derp = {}
#     finne kandidater som matcher miRNA
#     finne sekvenser fra bowtie som aligner med miRNA
#     sequence_tree[tree][five_interval.begin:three_range_end]
    candidate_already = 0
    seqs_count = 0
    no_seqs = 0
    candidate_count = 0
    has_candidates = []
    has_seqs = []
    for mirna_loki in mirna_hits:
        
        
        miRNAid = ">" + mirna_loki[0]
#         print "candidate", miRNAid, hairpinID_to_mature[miRNAid]
        strand_dir = mirna_loki[1]
        chromosome = mirna_loki[2].split("|")[3]
        genome_offset = int(mirna_loki[3])

        mature_pos = hairpinID_to_mature[miRNAid]
        begin_5p = mature_pos[0] + genome_offset if mature_pos[0] >= 0 else mature_pos[0]
        end_5p =  mature_pos[1] + genome_offset if mature_pos[1] >= 0 else mature_pos[1]
        begin_3p =  mature_pos[2] + genome_offset if mature_pos[2] >= 0 else mature_pos[2]
        end_3p =  mature_pos[3] + genome_offset if mature_pos[3] >= 0 else mature_pos[3]
        
        tree = candidate_tree[chromosome]
        if not tree:
            print "no:", chromosome
            continue
        candidates = tree[begin_5p:end_3p]
#         print tree, 
        
        #TODO: not -1 positions 5p, 3p

        if candidates:
            has_candidates.append(miRNAid)
            
            print len(candidates), begin_5p, end_3p
            candidate_already += 1
            
#             print
#             print "candidate match:"
#             print begin_5p
#             print end_5p
#             print begin_3p
#             print end_3p
            
            for candidate in candidates:
                print candidate.data.chromosome, "-", str(candidate.data.pos_5p_begin)
                candidate_count += 1
                similar_pos = 0
                hashval = candidate.data.chromosome + str(candidate.data.pos_5p_begin)
                if candidate.data.pos_5p_begin == begin_5p:
                    similar_pos += 1
                if candidate.data.pos_5p_end == end_5p:
                    similar_pos += 1
                if candidate.data.pos_3p_begin == begin_3p:
                    similar_pos += 1
                if candidate.data.pos_3p_end == end_3p:
                    similar_pos += 1
                if similar_pos >= 2:
                    hashval = candidate.data.chromosome + str(candidate.data.pos_5p_begin)
                    candidate_to_miRNAid[hashval] = miRNAid
#                 candidates_at_miRNA.add(candidate)
                derp[hashval] = miRNAid
            continue
        
#         # not candidate already
        tree = sequence_tree[chromosome]
        sequences = tree[begin_5p:end_3p]
        if sequences:
            has_seqs.append(miRNAid)

    print
    print "candidate is miRNA repeats:", candidate_already, candidate_count
    print "set of candidates:", len(set(candidate_to_miRNAid.iterkeys())), len(list(candidate_to_miRNAid.iterkeys()))
    print "derp candidates:", len(set(derp.iterkeys())), list(derp.iterkeys())[0:5]
    print "miRNA aligning candidates:", len(set(has_candidates)), has_candidates[0]
    print "miRNA aligning sequences:", len(set(has_seqs)), has_seqs[0]
    print "both (pls no):", set(has_seqs) & set(has_candidates)

#     for k,v in derp.iteritems():
#         print k,v
#     assert False
    return candidate_to_miRNAid


# 
#     def __init__(self, chromosome, strand_dir, begin_5, end_5,
#                  begin_3, end_3, mapped_sequences):








