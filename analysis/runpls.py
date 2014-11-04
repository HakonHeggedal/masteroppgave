# import os
# from bowtie import bowtie_get

import time
from inputs import merge
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



def _input_collapsed():
    pass


def _align_bowtie(bowtie_output, collapsed_seqs):
    import os
    from subprocess import check_output
    os.environ['BOWTIE_INDEXES'] = "/home/hakon/Skrivebord/h_sapiens_37_asm.ebwt"
    
    
    
    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", "-v 0", "-a", "-m 10",  human_index, collapsed_seqs, bowtie_output]
#     bowtie_str = "bowtie -f -v 0 -a -m 10 h_sapiens_37_asm SRR797062.fa test.map"
#     bowtie_str = "bowtie -f -v 0 h_sapiens_37_asm SRR797062.fa test.map"
#     
#     cmd_simple = ["bowtie", "-f", "-v 0", "h_sapiens_37_asm", "-c", "ACGTACGTACGT"]
#     cmd_str = "bowtie -f -v 0 h_sapiens_37_asm -c ACGTACGTACGT"
#     cmd_exitsts = "bowtie-inspect -n h_sapiens_37_asm"
#     
    
    print check_output(bowtie_cmds)


#     print check_output(cmd_exitsts, shell=True)
#     print check_output(cmd_str, shell=True)
    
#     print check_output(["echo", "testest"])
#     print check_output(["bowtie", "h_sapiens_37_asm", "-c", "ACGTACGTACGT"])
#     print check_output("bowtie h_sapiens_37_asm -c ACGTACGTACGT",  shell=True )
#     print "running bowtie"


#     outs = check_output(bowtie_str, shell=True).strip()
#     print bowtie_cmds
#     outs = check_output(bowtie_cmds, shell=True).strip()
    
    
    
#     outs = check_output(["bowtie", human_index, "-c", "ACGTACGT", "ts.txt"], shell=True).strip()
#     outs = check_output(bowtie_cmds, shell=True)
#     print "run commands2"        
# #     print answers
# #     print "errors", errors
#     parts = outs.split("\t")
#     print len(outs), outs
#     print len(parts), parts
#      
    
#     bowtie_get.runcommand(bowtie_cmds)
#     bowtie_get.runcommand(["bowtie","-f", "h_sapiens_37_asm", "SRR797062.fa", "test.map"])



def main():
    start_time =time.clock()
    print "starting"
    
    print "merging collapsed files"
    fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed",
                   "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed"]
    
#     fasta_file = "SRR797062.fa"

#     merge collapsed input files
    dict_collapsed = merge.collapse_collapsed(fasta_files)
    
#     split small and larger sequences
    reads, reads_count, small_reads, small_reads_count = merge.filter_seqeunces(dict_collapsed, 18)
    
#     write reads to file
    all_reads_file = "all.collapsed"
    merge.write_collapsed(all_reads_file, reads, reads_count)
    
    
#     aligning to genome using bowtie
    bowtie_output = "bowtie_out.map"
    _align_bowtie(bowtie_output, all_reads_file)
    print "finished bowtie in ", time.clock() - start_time, " seconds" 
    
    
    fixed_lines = [line.strip().split("\t") for line in open(bowtie_output)] 
    print "read positions in ", time.clock() - start_time, " seconds"
    
#     using sequence tree to find possible candidates
    candidate_tree, sequence_tree, candidates, seq_to_candidates = interval_tree_search.find_candidates(fixed_lines)
    print "found candidates in ", time.clock() - start_time, " seconds"
    print "bowtie hits", len(fixed_lines)
    print "candidate tree", len(candidate_tree)
    print "candidates", len(candidates)
    print "sequence tree", len(sequence_tree)
    print "mapped seqs", len(candidates[0].all_mapped_sequences)
    
#     candidate_list = gene.find_all(candidates)
    gene.include_padding(candidates)
    print "padded all candidates in ", time.clock() - start_time, " seconds"
    
#     TODO: run annotated miRNA through bowtie + others  
    
    
    
    
#     sequence_freq = reads.readcollapsed(fasta_file)
#     print len(sequence_freq)    
    
    
    # run and set vienna RNAfold + energy on all candidates
    vienna.energy_fold(candidates)
    
    
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
    entropy(candidates)
#      
#     heterogenity (position counting)
    heterogenity.heterogenity(candidates)
#      
#     candidate quality: nr of sequence hits / all candidate hits for given sequences
    quality.candidate_quality(candidates, seq_to_candidates)
    
    print "finished features in ", time.clock() - start_time, " seconds"
    

    
    
    


# if __name__ == "__main__":
#     main()

main()