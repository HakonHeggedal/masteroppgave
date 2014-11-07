# import os
# from bowtie import bowtie_get

import time
start_time = time.clock()

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

#     fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed"]
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
    
    

    
#     get mirbase entries
    miRNAid_to_hairpin = mirbase.read_miRNA("mature.fa", "hairpin.fa")
    print "got all miRNAs from mirbase files", len(miRNAid_to_hairpin)
    
#     run write micro rnas to file
    miRNA_file_name = "mirnas.fa"
    mirbase.write_miRNA(miRNAid_to_hairpin, miRNA_file_name)
    print "wrote human miRNAs to file", time.clock() - start_time, " seconds"
    
#     run bowtie to find miRNA positions
    miRNA_bowtie_output = "miRNA.map"
    _align_bowtie(miRNA_bowtie_output, miRNA_file_name)
    
    print "aligned miRNAs in", time.clock() - start_time, " seconds"
    
    miRNA_bowtie_hits = [line.strip().split("\t") for line in open(miRNA_bowtie_output)] 
    
    unique_mirna_hits = set([x[0] for x in miRNA_bowtie_hits])
    
    print "miRNA bowtie hits:", len(miRNA_bowtie_hits)
    print "unique miRNA hits:", len(unique_mirna_hits)

    
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
    
    
    print "create candidate structure from miRNAs"
    miRNA_candidates = []
    for mirna_loki in miRNA_bowtie_hits:
        miRNAid = mirna_loki[0]
        strand_dir = mirna_loki[1]
        chromosome = mirna_loki[2].split("|")[3]
        hairpin = miRNAid_to_hairpin[">"+miRNAid].hairpin
        begin_5 = int(mirna_loki[3])
        end_5 = 0
        begin_3 = 0
        end_3 = begin_5 + len(hairpin)
        candidate_sequences = None
        
        
        mirna_candidate = structure.Candidate(chromosome,
                                     strand_dir,
                                     begin_5,
                                     end_5,
                                     begin_3,
                                     end_3,
                                     candidate_sequences)

        mirna_candidate.set_hairpin_padding(hairpin, "", -1)
        print hairpin
        miRNA_candidates.append(mirna_candidate)
        
    
    print "align miRNAs to other sequences"
    candidate_union = interval_tree_search.align_miRNAs(miRNA_candidates, candidate_tree, sequence_tree)
    

    
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
    

    candidate_array = vectorize.candidates_to_array(candidates)
    
    


# if __name__ == "__main__":
#     main()

main()