# import os
# from bowtie import bowtie_get
from candidates import interval_tree_search, overhang, heterogenity
import time
from genes import gene
from subprocess import check_output
import reads
from candidates.vienna import energy_fold
from candidates.tailing import three_prime_au
from candidates.entropy import entropy

def main():
    start_time =time.clock()
    print "starting"
    fasta_file = "SRR797062.fa"
    bowtie_output = "test.map"
    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", "-v 0", "-a", "-m 10",  human_index, fasta_file, bowtie_output]
    bowtie_str = "bowtie -f -v 0 -a -m 10 h_sapiens_37_asm SRR797062.fa test.map"
    bowtie_str = "bowtie -f -v 0 h_sapiens_37_asm SRR797062.fa test.map"
    
    
    
#     print "running bowtie"
#     outs = check_output(bowtie_str, shell=True).strip()
# #     outs = check_output(["bowtie", human_index, "-c", "ACGTACGT", "ts.txt"], shell=True).strip()
#     print "run commands2"        
# #     print answers
# #     print "errors", errors
#     parts = outs.split("\t")
#     print len(outs), outs
#     print len(parts), parts
#      
    
#     bowtie_get.runcommand(bowtie_cmds)
#     bowtie_get.runcommand(["bowtie","-f", "h_sapiens_37_asm", "SRR797062.fa", "test.map"])
    
    print "finished bowtie in ", time.clock() - start_time, " seconds" 
#     unfixed_lines = open(bowtie_output).readlines() # read
#     fixed_lines = [line.strip().split("\t") for line in unfixed_lines] 
    fixed_lines = [line.strip().split("\t") for line in open(bowtie_output)] 
    print "read positions in ", time.clock() - start_time, " seconds" 
    candidate_tree, sequence_tree, candidates, seq_to_candidates = interval_tree_search.find_candidates(fixed_lines)
    
    print "bowtie hits", len(fixed_lines)
    print "candidate tree", len(candidate_tree)
    print "candidates", len(candidates)
    print "sequence tree", len(sequence_tree)
    print "mapped seqs", len(candidates[0].all_mapped_sequences)
    
    print "found candidates in ", time.clock() - start_time, " seconds"
    
#     candidate_list = gene.find_all(candidates)
    gene.include_padding(candidates)
    print "candidates??", len(candidates)
    
    

    print "found all loki in ", time.clock() - start_time, " seconds"
    
    
    sequence_freq = reads.readcollapsed(fasta_file)
    print len(sequence_freq)
    
    # run and set vienna RNAfold + energy on all candidates
    energy_fold(candidates)
    
#     # A/U ends for all candidates
#     three_prime_au(candidates, sequence_freq)
#     
#     # 5' and 3' alignment overhang
#     overhang(candidates)
#     
#     # degree of entropy in structure and nucleotides
#     entropy(candidates)
    
    # heterogenity (position counting)
    heterogenity.frequency_counting(candidates, 5)
    
    
    print "finished all in ", time.clock() - start_time, " seconds"
    
    
    
    


# if __name__ == "__main__":
#     main()

main()