'''
Created on 29. okt. 2014

@author: hakon
'''


def collapse_collapsed(collapsed_files, min_len=7):
    ''' merge files in files list, return as lists'''
    
    all_sequences = {}
    
    for filename in collapsed_files:
        
        current_seqs = {}
        total_seqs = 0
        
        with open(filename) as fasta_file:

            count = 0
            for line in fasta_file:
                line = line.strip()
                if line[0] == ">":
                    count = int(line.split("-")[1])
                elif len(line) >= min_len and count > 0 and _legal_DNA(line):
                    current_seqs[line] = count
                    total_seqs += count
                    count = 0
        
        for s, c  in current_seqs.iteritems():
            c_norm = c * 1000000.0 / total_seqs
            
            if s in all_sequences:
                all_sequences[s] += c_norm
            else:
                all_sequences[s] = c_norm
    
    return all_sequences

#     ''' merge files in files list, return as lists'''
#     
#     all_sequences = {}
#     
#     for filename in collapsed_files:
#         
#         with open(filename) as fasta_file:
#             
#             count = 0
#             for line in fasta_file:
#                 line = line.strip()
#                 if line[0] == ">":
#                     count = int(line.split("-")[1])
#                 elif len(line) >= min_len and count > 0 and _legal_DNA(line):
#                     if line in all_sequences:
#                         all_sequences[line] += count
#                     else:
#                         all_sequences[line] = count
#                     count = 0
#     
#     return all_sequences
#     

def filter_seqeunces(all_sequences, min_read_len=18):
    
    small_reads = []
    small_reads_count = []
    reads = []
    reads_count = []
    
    for seq, count in all_sequences.iteritems():
        
        if len(seq) < min_read_len:
            small_reads.append(seq)
            small_reads_count.append(count)
        else:
            reads.append(seq)
            reads_count.append(count)
    
    print len(reads), len(small_reads), len(reads_count), len(small_reads_count)
    
    return reads, reads_count, small_reads, small_reads_count


def write_collapsed(output_name, reads, reads_count):
    
    with open(output_name, "w") as file_out:
        
        for i, (read, count) in enumerate(zip(reads, reads_count)):
            name = ">" + str(i) + "-" + str(count)
            
            file_out.write(name + "\n")
            file_out.write(read + "\n")
            


def _legal_DNA(sequence, chars=set("ACGT")):
    
    for char in sequence:
        if char not in chars:
#             print char
            return False
    
    return True


