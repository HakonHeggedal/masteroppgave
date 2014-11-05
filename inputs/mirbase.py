'''
Created on 29. okt. 2014

@author: hakon
'''

class Micro_RNA:
    def __init__(self,hairpin):
        self.hairpin = hairpin
        self.pos_5_begin = None
        self.pos_5_end = None
        self.pos_3_begin = None
        self.pos_3_end = None
    
    def set_five_prime(self,pos_5_begin, pos_5_end):
        self.pos_5_begin = pos_5_begin
        self.pos_5_end = pos_5_end
    
    def set_three_prime(self, pos_3_begin, pos_3_end):
        self.pos_3_begin = pos_3_begin
        self.pos_3_end = pos_3_end

def _add_position(miRNA, startpos, endpos, is_5p, is_3p):
    assert startpos >= 0
    if is_5p:
        miRNA.set_five_prime(startpos, endpos)
    elif is_3p:
        miRNA.set_three_prime(startpos, endpos)
    elif startpos+5 < len(miRNA.hairpin) / 2: # not given, but aligns with 5' end
        miRNA.set_five_prime(startpos, endpos)
    else:
        miRNA.set_three_prime(startpos, endpos)
        


def read_miRNA(mature_file, hairpin_file):
    
    micro_rnas = {}
    
    with open(hairpin_file) as hairpins:
        hairpin = ""
        is_human = False
        name = ""
        for hairpin_seq in hairpins:
            hairpin_seq = hairpin_seq.lower().strip()
            
            if hairpin_seq[0] == ">":

                if is_human:  # previous entry
                    micro_rnas[name] = Micro_RNA(hairpin.replace("u", "t"))
                
                if hairpin_seq[1:4] == "hsa":
#                     print mature_seq[1:4]
                    is_human = True
                    name = hairpin_seq.split()[0]
#                     print "hairpin \"%s\"" % name, len(name)
                else:
                    is_human = False
                hairpin = ""
                    
            elif is_human:
                hairpin += hairpin_seq



    no_hits = 0
#     print "finding mature seqs"
    with open(mature_file) as mature_seqs:
        
        is_human = False
        name = ""
        for mature_seq in mature_seqs:
            mature_seq = mature_seq.lower().strip()
            mature_seq = mature_seq.replace("u", "t")
            if mature_seq[0] == ">":
                if mature_seq[1:4] == "hsa":
                    is_human = True
                    name = mature_seq.split()[0]
#                     print name
                    
                else:
                    is_human = False
                    
            elif is_human:
                found_hairpin = False
                is_3p = False
                is_5p = False
                
                if name.endswith("3p"):
                    is_3p = True
                    name = name[:-3]
                elif name.endswith("5p"):
                    is_5p = True
                    name = name[:-3]
                
#                 print name, len(name)
#                 print "after: %s" %name
#                 print "is \"%s\" in dict?" % name, name in micro_rnas
                if name in micro_rnas:
                    pos = micro_rnas[name].hairpin.find(mature_seq) # mature seq. must be in hairpin.
                    if pos >= 0:
                        found_hairpin = True
                        endpos = pos + len(mature_seq)
                        _add_position(micro_rnas[name], pos, endpos, is_5p, is_3p)

#                 elif not found_hairpin:

                x = 1
                
                if name+"-"+str(x) not in micro_rnas:
                    x = 2
                
                while name+"-"+str(x) in micro_rnas:
                    pos = micro_rnas[name+"-"+str(x)].hairpin.find(mature_seq)
                    if pos >= 0:
                        found_hairpin = True
                        endpos = pos + len(mature_seq)
                        _add_position(micro_rnas[name+"-"+str(x)],
                                      pos, endpos, is_5p, is_3p)
                    x += 1
#                     else:
#                         if not found_hairpin:
                c = "a"
                while name+c in micro_rnas:
                    pos = micro_rnas[name+c].hairpin.find(mature_seq)
                    if pos >= 0:
                        found_hairpin = True
                        endpos = pos + len(mature_seq)
                        _add_position(micro_rnas[name+c],pos, endpos, is_5p, is_3p)
                    c += chr(ord(c) + 1)
                                     
                if not found_hairpin:
                    no_hits += 1
                    
                    print "\tno hit", name, x, c, is_3p, is_5p
    
    print "finished assembling miRNAs"
    return micro_rnas
    


def write_miRNA(mi_rnas, filename):
    

    with open(filename, "w") as outfile:
        for name, miRNA in mi_rnas.iteritems():
            outfile.write( name + "\n")
            outfile.write(miRNA.hairpin.upper() + "\n")     
            

    

# read_miRNA("mature.fa", "hairpin.fa")
























