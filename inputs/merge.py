'''
Created on 29. okt. 2014

@author: hakon
'''


def merge_collapsed_fasta(fasta_files):
    ''' merge files in files list, return as lists'''
    
    all_sequences = {}
    
    for filename in fasta_files:
        
        with open(filename) as fasta_file:
            
            count = 0
            for line in fasta_file:
                line = line.strip()
                if line[0] == ">":
                    count = int(line.split("-")[1])
                elif len(line) > 10 and count > 0 and _legal_DNA(line):
                    if line in all_sequences:
                        all_sequences[line] += count
                    else:
                        all_sequences[line] = count
                    count = 0
    
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
    
    print len(reads), len(small_reads), len(reads_count), len(small_reads_count)
    
    return reads, reads_count, small_reads, small_reads_count




def _legal_DNA(sequence, chars=set("ACGT")):
    
    for char in sequence:
        if char not in chars:
#             print char
            return False
    
    return True


fastas = ["SRR797060.collapsed", "SRR797061.collapsed", "SRR797062.collapsed", "SRR797063.collapsed"]

merge_collapsed_fasta(fastas)