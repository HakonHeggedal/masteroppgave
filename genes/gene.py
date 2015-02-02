'''
Created on 6. okt. 2014

@author: hakon
'''
import os

def include_padding(candidate_list, padding=40 ):
    
    # such hack...
    p =  os.getcwd()
    os.chdir(os.pardir)
    os.chdir("genes")
        
    i = 0
    print "candidates:", len(candidate_list)

    while i < len(candidate_list):
        
        candidate = candidate_list[i]
        gene_name = candidate.chromosome

        file_name = gene_name + ".fa"
        
        print "reading from", file_name
        
        with open(file_name) as chr_file:
#             print "reading from", file_name
            
            _ = chr_file.readline() # first line is info, not used
            chr_lines = chr_file.readlines() # current genome
            line_len = len(chr_lines[0].strip())
            
#             print line_len
#             print len(chr_lines) * line_len, len(chr_lines)
            
            while True:
                # fix every candidate
                
                read_start = candidate.pos_5p_begin - padding
                start_line = read_start // line_len
                start_pos = read_start % line_len

                read_end = candidate.pos_3p_end + padding
                end_line = read_end // line_len
                end_pos = read_end % line_len
                
                #assembling sequence
                padded = chr_lines[start_line][start_pos:].strip()
                for x in xrange(start_line+1, end_line):
#                     print "-", start_line, end_line
#                     print "-"
                    padded += chr_lines[x].strip()
                padded += chr_lines[end_line][:end_pos].strip() if start_line != end_line else ""
                
                padded = padded.upper()
                padded = padded.strip()
                padded = padded.strip("\n")
                
                hairpin = padded[padding:-padding]
                
#                 print "hairpin",hairpin


#                 has_seqs = False
#                 for s in candidate.mapped_sequences:
#                     has_seqs = has_seqs or s.data[2] in padded
#                 assert has_seqs
                  
#                 print s.data[2]
#                 if s.data[2] not in padded:
#                     print "\t NO MATCH!!!!!!", s.data[2], padded[:10], read_start

#                 if interval.data[2] not in padded or interval.data[4] not in padded:
#                     fails += 1


                if not hairpin:
                    print "hairpin", i, hairpin
                    assert hairpin
                assert padded
                assert padding
                
#                 if len(padded) != read_end - read_start:
#                     print candidate.chromosome, file_name
#                     print read_end - read_start, len(padded)
#                     print padded
#                     print len(chr_lines), end_line
#                     print len(chr_lines[start_line].strip()), chr_lines[start_line].strip()
#                     print len(chr_lines[end_line].strip()), chr_lines[end_line].strip()
#                 assert len(padded) == read_end - read_start
                candidate.set_hairpin_padding(hairpin, padded)
                
#                 print "candidate.hairpin", candidate.hairpin
                i += 1
                if  i >= len(candidate_list):
                    break
                
                if candidate_list[i].chromosome != gene_name:
                    print "wrong chromosome error", candidate_list[i].chromosome, gene_name
                    assert False
                    break
                
                candidate = candidate_list[i]
                gene_name = candidate.chromosome
                
#                 print "123"
                
                
    os.chdir(p)
    
#     assert False
    return 
        


        
        