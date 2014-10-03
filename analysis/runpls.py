import os
from bowtie import bowtie_get
from candidates import interval_tree_search



def main():
    fasta_file = "SRR797062.fa"
    bowtie_output = "test.map"
    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", human_index, fasta_file, bowtie_output]
    
    bowtie_get.runcommand(bowtie_cmds)
#     bowtie_get.runcommand(["bowtie","-f", "h_sapiens_37_asm", "SRR797062.fa", "test.map"])
    
    unfixed_lines = open(bowtie_output).readlines() # read
    fixed_lines = [line.strip().split("\t") for line in unfixed_lines] # make 
    
    candidates, sequences = interval_tree_search.find_candidates(fixed_lines)
    
    print len(candidates)
    print len(sequences)
    
    





















# if __name__ == "__main__":
#     main()

main()