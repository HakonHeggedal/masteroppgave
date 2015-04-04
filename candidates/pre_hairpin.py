'''
Created on 4. apr. 2015

@author: hakon
'''




def hairpin_cutoff_seqs(candidates):
    
    print "hairpin cutoff seqs", len(candidates)
    counter = 0
    for c in candidates:
        if len(c.sequences_before):
             
            cutoff_pos = c.pos_5p_begin - 1
            
             
            seqs = [(s.end, float(s.data[1].split("-")[1]) ) for s in c.sequences_before]
            seqs = [(start, val) for (start, val) in seqs if start <= cutoff_pos]
            
#             teste = [(s.begin, s.end, float(s.data[1].split("-")[1]) ) for s in c.sequences_before]
 
            sum_sequences = sum( copies for _pos, copies in seqs)
            sum_cutoff = sum( copies for pos, copies in seqs if pos == cutoff_pos)
             
            score = sum_cutoff / sum_sequences if sum_sequences else 0.0
             
            if sum_cutoff:
                print "\n5p is happening", cutoff_pos
                print c.pos_5p_begin, c.pos_3p_end
                print sum_cutoff, sum_sequences, score
                print sorted(seqs)
                counter += 1
            c.stops_before_5p = score
            
        if len(c.sequences_after):
            
            cutoff_pos = c.pos_3p_end + 1
            
            seqs = [(s.begin, float(s.data[1].split("-")[1]) ) for s in c.sequences_after]
            seqs = [(begin, val) for (begin, val) in seqs if begin >= cutoff_pos]
            
#             teste = [(s.begin, s.end, float(s.data[1].split("-")[1]) ) for s in c.sequences_after]
            
            sum_sequences = sum( copies for _pos, copies in seqs)
            sum_cutoff = sum( copies for pos, copies in seqs if pos == cutoff_pos)
            
            score = sum_cutoff / sum_sequences if sum_sequences else 0.0
        
        
        
            if sum_cutoff:
                print "\n3p is happening", cutoff_pos
                print c.pos_5p_begin, c.pos_3p_end
                print sum_cutoff, sum_sequences, score
                print sorted(seqs)
                counter += 1
                    
            c.starts_after_3p = score
            
            
    print "count cutoffs, total candidates:", counter, len(candidates)
    
    