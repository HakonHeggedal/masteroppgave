'''
Created on 29. okt. 2014

@author: hakon
'''


def merge_collapsed_fasta(fasta_files):
    ''' merge files in files list, return as lists'''
    
    all_sequences = {}
    
    for filename in fasta_files:
        
        with open(filename) as fasta_file:
            
            for line in fasta_file:
                line = line.strip()
                if line[0] == ">":
                    count = int(line.split("-")[1])
                elif len(line) > 10:
                    
                    if line in all_sequences:
                        all_sequences[line] += count
                    else:
                        all_sequences[line] = count
    
    small_reads = []
    small_reads_count = []
    reads = []
    reads_count = []
    
    for seq, count in all_sequences.iteritems():
        
        if len(seq) < 18:
            small_reads.append(seq)
            small_reads_count.append(count)
        else:
            reads.append(seq)
            reads_count.append(count)
            
    
    return reads, reads_count, small_reads, small_reads_count
    
