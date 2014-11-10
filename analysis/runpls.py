# import os
# from bowtie import bowtie_get

import time
start_time = time.clock()
import numpy
from sklearn import svm

from inputs import merge
from inputs import mirbase
from genes import gene
from candidates import interval_tree_search
from candidates import heterogenity
from candidates import vienna
from candidates import tailing
from candidates import entropy
from candidates import quality
from candidates import structure
from candidates import overhang

from ml import vectorize
from ml import learn




def _align_bowtie(bowtie_output_file, collapsed_seq_file):
    from subprocess import check_output
    import os # hack: setting path to bowtie index
    os.environ['BOWTIE_INDEXES'] = "/home/hakon/Skrivebord/h_sapiens_37_asm.ebwt"

    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", "-v 0", "-a", "-m 10",
                   human_index, collapsed_seq_file, bowtie_output_file]

    print check_output(bowtie_cmds)





def main():
#     start_time = time.clock()
    print "starting"
    
    print "merging collapsed files"
    fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed",
                    "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed"]

    fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed"]
#     fasta_file = "SRR797062.fa"

#     merge collapsed input files
    dict_collapsed = merge.collapse_collapsed(fasta_files)
    
#     split small and larger sequences
    reads, reads_count, small_reads, small_reads_count = merge.filter_seqeunces(dict_collapsed, 18)
    print "reads:", len(reads), "small:", len(small_reads), len(small_reads_count)
#     write reads to file
    all_reads_file = "all.collapsed"
    merge.write_collapsed(all_reads_file, reads, reads_count)
    
    
#     aligning to genome using bowtie
    bowtie_output = "bowtie_out.map"
    _align_bowtie(bowtie_output, all_reads_file)
    print "finished bowtie in ", time.clock() - start_time, " seconds" 
    
#     read genome alignment from bowtie
    fixed_lines = [line.strip().split("\t") for line in open(bowtie_output)] 
    print "read positions in ", time.clock() - start_time, " seconds"
    
    
    
    hairpin_file = "hairpin.fa"
    mature_seq_file = "mature.fa"
    miRNA_file_name = "mirnas.fa"
    
    print "loading miRNA hairpins:"
    
    hsa_to_hairpin, other_to_hairpin = mirbase.read_miRNA_fasta(hairpin_file)
    hsa_to_mature, other_to_mature = mirbase.read_miRNA_fasta(mature_seq_file)
    
    miRNA_species = mirbase.similar_hairpins(hsa_to_hairpin, other_to_hairpin)
    
    hairpinID_to_mature = mirbase.combine_hairpin_mature(hsa_to_hairpin, hsa_to_mature)
    
    for k,v in miRNA_species.iteritems():
        if v > 0:
            print k,v
    
    
#     assert False
#     get mirbase entries
#     miRNAid_to_hairpin = mirbase.read_miRNA("mature.fa", "hairpin.fa")
#     print "got all miRNAs from mirbase files", len(miRNAid_to_hairpin)


    
#     miRNA_species = mirbase.mirna_copies(miRNAid_to_hairpin.values(), "hairpin.fa")
    
#     for k,v in miRNA_species.iteritems():
#         print v
    
#     run write micro rnas to file

    mirbase.write_miRNA(hsa_to_hairpin, miRNA_file_name)
    print "wrote human miRNAs to file", time.clock() - start_time, " seconds"
    
#     run bowtie to find miRNA positions
    miRNA_bowtie_output = "miRNA.map"
    _align_bowtie(miRNA_bowtie_output, miRNA_file_name)
    
    print "aligned miRNAs in", time.clock() - start_time, " seconds"
    
    miRNA_bowtie_hits = [line.strip().split("\t") for line in open(miRNA_bowtie_output)] 
    
    unique_mirna_hits = set([x[0] for x in miRNA_bowtie_hits])
    
    print "miRNA bowtie hits:", len(miRNA_bowtie_hits)
    print "unique miRNA hits:", len(unique_mirna_hits)

    
    print "create candidate structure from miRNAs"
    miRNA_candidates = []
    for mirna_loki in miRNA_bowtie_hits:
        miRNAid = ">" + mirna_loki[0]
#         print "candidate", miRNAid, hairpinID_to_mature[miRNAid]
        strand_dir = mirna_loki[1]
        chromosome = mirna_loki[2].split("|")[3]
        genome_offset = int(mirna_loki[3])
        hairpin = hsa_to_hairpin[miRNAid]
        mature_pos = hairpinID_to_mature[miRNAid]
        begin_5 = mature_pos[0] + genome_offset if mature_pos[0] >= 0 else mature_pos[0]
        end_5 =  mature_pos[1] + genome_offset if mature_pos[1] >= 0 else mature_pos[1]
        begin_3 =  mature_pos[2] + genome_offset if mature_pos[2] >= 0 else mature_pos[2]
        end_3 =  mature_pos[3] + genome_offset if mature_pos[3] >= 0 else mature_pos[3]
        candidate_sequences = None
        
        
        mirna_candidate = structure.Candidate(chromosome,
                                     strand_dir,
                                     begin_5,
                                     end_5,
                                     begin_3,
                                     end_3,
                                     candidate_sequences)

        mirna_candidate.set_hairpin_padding(hairpin, "", -1)
        mirna_candidate.set_hairpin_pos(genome_offset, genome_offset + len(hairpin) )


        miRNA_candidates.append(mirna_candidate)
    
    print "candidates:", len(miRNA_candidates)
    
    
    

    
#     assert False
    
#     using sequence tree to find possible candidates
    candidate_tree, sequence_tree, candidates, seq_to_candidates = interval_tree_search.find_candidates(fixed_lines)
    
    print "found candidates in ", time.clock() - start_time, " seconds"
    print "bowtie hits", len(fixed_lines)
    print "candidate tree", len(candidate_tree)
    print "candidates", len(candidates)
    print "sequence tree", len(sequence_tree)
    print "mapped seqs", len(candidates[0].all_mapped_sequences)
    
#     candidate_list = gene.find_all(candidates)

# 0            1   2[0] [1]      [2] [3]
# ['1-15830', '-', 'gi|224589818|ref|NC_000006.11|',
#         NC_000006.11
        
    gene.include_padding(candidates)
    print "padded all candidates in ", time.clock() - start_time, " seconds"
        
    
    print "align miRNAs to other sequences"
    candidate_to_miRNA = interval_tree_search.align_miRNAs(miRNA_candidates, candidate_tree, sequence_tree)
    

    
#     sequence_freq = reads.readcollapsed(fasta_file)
#     print len(sequence_freq)    
    
    
    # run and set vienna RNAfold + energy on all candidates
    vienna.energy_fold(candidates) # slow?
    
    
#     not_mapped_reads = [structure.Sequence(i,n,read) for i,(read,n) in 
#                         enumerate(zip(reads, reads_count))
#                         if read not in seq_to_candidates]
#     
#     print len(not_mapped_reads)
#     print len(seq_to_candidates)
#     print len(not_mapped_reads) + len(seq_to_candidates)
#     print len(reads)
#     
# #     A/U ends for all remaining candidates
#     tailing_au(candidates, not_mapped_reads)
# #     
#     #TODO: 5' and 3' alignment overhang
#     overhang.find_overhang(candidates)
#     
#     degree of entropy in structure and nucleotides
    entropy.entropy(candidates)
#      
#     heterogenity (position counting)
    heterogenity.heterogenity(candidates)
#      
#     candidate quality: nr of sequence hits / all candidate hits for given sequences
    quality.candidate_quality(candidates, seq_to_candidates)
    
    print "finished features in ", time.clock() - start_time, " seconds"
    
    
    miRNAs = []
    miRNA_annotations = [] 
    only_candidates = []
    
    for candidate in candidates:
        print candidate, candidate in candidate_to_miRNA
        if candidate in candidate_to_miRNA:
            print "ok mirna!"
            miRNAs.append(candidate)
            
            
    
    learn_miRNAs = vectorize.candidates_to_array(miRNAs)
    class_miRNAs = numpy.array(miRNA_annotations)
    candidate_array = vectorize.candidates_to_array(only_candidates)
    
    learner = svm.SVC()
    learner.fit(learn_miRNAs, class_miRNAs)
    learner.predict(only_candidates)


# if __name__ == "__main__":
#     main()

main()