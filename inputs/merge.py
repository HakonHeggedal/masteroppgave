'''
Created on 29. okt. 2014

@author: hakon
'''


def collapse_collapsed(collapsed_files, min_len, min_count):
    ''' merge files in files list, return as lists'''
    
    all_sequences = {}
    
    allll = 0
    
    for filename in collapsed_files:
        print filename
        
        current_seqs = {}
        
        total_count = 0
        used_count = 0
        
        total_norm_seqs = 0
        used_seq = 0
        
        
        all_seqs = 0
#         legals = 0
#         large = 0
#         counted = 0

        
        with open(filename) as fasta_file:
            
            count = 0
            for line in fasta_file:
                line = line.strip()
                if not line:
                    count = 0
                    continue
                if line[0] == ">":
                    count = int(line.split("-")[1])
                    total_count += count
                    allll += count
                else:
                    
                    all_seqs += 1
#                     if _legal_DNA(line):
#                         legals += 1
#                         
#                     if len(line) >= min_len:
#                         large += 1
#                     
#                     if count >= min_count:
#                         counted += 1
                        
                    
                    if len(line) >= min_len and count >= min_count and _legal_DNA(line):
                        current_seqs[line] = count
                        count = 0
        
        for s, c  in current_seqs.iteritems():
            c_norm = c * 1000000.0 / total_count
            
            total_norm_seqs += c_norm
            used_seq += 1
            used_count += c
            
            if s in all_sequences:
                all_sequences[s] += c_norm
            else:
                all_sequences[s] = c_norm
                
                
#         print all_seqs, used_seq*1.0/ all_seqs, "\t", total_count, used_count*1.0/total_count, "\t", filename
#         print total_count, all_seqs, "\t", legals, large, counted, "\t", used_seq, filename
#         total_count, total_norm_seqs,
    print allll
    assert 0
    return all_sequences


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
    
    print "\tfilter by length: min", min_read_len
    print "\t", len(reads), len(small_reads), len(reads_count), len(small_reads_count)
    
    return reads, reads_count, small_reads, small_reads_count


def write_collapsed(output_name, reads, reads_count):
    
    with open(output_name, "w") as file_out:
        
        for i, (read, count) in enumerate(zip(reads, reads_count)):
            name = ">" + str(i) + "-" + str(count)
            
            file_out.write(name + "\n")
            file_out.write(read + "\n")
            

def one_large_fasta(fasta_files, write_file_name):

    all_sequences = {}
    print 123
    
    for filename in fasta_files:
        
        current_seqs = {}
        
        total_count = 0
        total_norm_seqs = 0


        print filename
        
        with open(filename) as fasta_file:
            
            count = 0
            for line in fasta_file:
                line = line.strip()
                if not line:
                    count = 0
                    continue
                if line[0] == ">":
                    count = int(line.split("-")[1])
                    total_count += count
                else:                        
                    
                    if len(line) >= 10 and count >= 0 and _legal_DNA(line):
                        current_seqs[line] = count
                        count = 0
        
        for s, c  in current_seqs.iteritems():
#             c_norm = c * 1000000.0 / total_count
            c_norm = c 
            
            total_norm_seqs += c_norm

            
            if s in all_sequences:
                all_sequences[s] += c_norm
            else:
                all_sequences[s] = c_norm
    
    print "writing to file"
    i = 0
    with open(write_file_name, "w") as write_file:
        
        for i, (seq, count) in enumerate(all_sequences.iteritems()):
            
            identifier = ">" + str(i) + "_x" + str(count) + "\n"
            
            write_file.write(identifier)
            write_file.write(seq + "\n")
            
    print "one large file. nr of sequences:", i
    
    return
 
             
    

def _legal_DNA(sequence, chars=set("ACGT")):
    
    for char in sequence:
        if char not in chars:
#             print char
            return False
    
    return True


