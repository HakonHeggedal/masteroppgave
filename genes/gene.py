'''
Created on 6. okt. 2014

@author: hakon
'''
import os

def include_padding(candidate_list_unsort, padding=40 ):
    
    # such hack...
    p =  os.getcwd()
    os.chdir(os.pardir)
    os.chdir("genes")
        
    i = 0
    print "candidates:", len(candidate_list_unsort)
    
#    sort the candidate list by genome... much faster
    candidate_list = candidate_list_unsort[:]
    candidate_list = sorted(candidate_list, key=lambda x: x.chromosome)

    while i < len(candidate_list):
        
        candidate = candidate_list[i]
        gene_name = candidate.chromosome

        file_name = gene_name + ".fa"
        
        print "reading from", file_name
        
        with open(file_name) as chr_file:
            
            _ = chr_file.readline() # first line is info, not used
            chr_lines = chr_file.readlines() # current genome
            line_len = len(chr_lines[0].strip())
            
#             print line_len
#             print len(chr_lines) * line_len, len(chr_lines)
            
            while True:
                
                read_start = candidate.hairpin_start - padding
                start_line = read_start // line_len
                start_pos = read_start % line_len

                read_end = candidate.hairpin_end + padding
                end_line = read_end // line_len
                end_pos = read_end % line_len
                
                #assembling sequence
                padded = chr_lines[start_line][start_pos:].strip()
                for x in xrange(start_line+1, end_line):
                    padded += chr_lines[x].strip()
                padded += chr_lines[end_line][:end_pos].strip() if start_line != end_line else ""
                
                padded = padded.upper()
                padded = padded.strip()
                padded = padded.strip("\n")
                
                hairpin = padded[padding:-padding]

                if not hairpin:
                    print "hairpin", i, hairpin
                    assert hairpin
                assert padded
                assert padding
                
                if candidate.hairpin:
                    if hairpin != candidate.hairpin:
                        print len(hairpin), hairpin
                        print len(candidate.hairpin), candidate.hairpin
                    assert hairpin == candidate.hairpin
                
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
#                     print "must change chromosome"
#                     assert False
                    break
                
                candidate = candidate_list[i]
                gene_name = candidate.chromosome
                
                
    os.chdir(p) # change folder again
    
    return 
        


        
        