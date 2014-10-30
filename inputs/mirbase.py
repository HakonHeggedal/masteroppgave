'''
Created on 29. okt. 2014

@author: hakon
'''

class Micro_RNA:
    
    def __init__(self,hairpin):
        self.hairpin = hairpin
    
    def set_five_prime(self,pos_5_begin, pos_5_end):
        self.pos_5_begin = pos_5_begin
        self.pos_5_end = pos_5_end
    
    def set_three_prime(self, pos_3_begin, pos_3_end):
        self.pos_3_begin = pos_3_begin
        self.pos_3_end = pos_3_end

def _add_position(miRNA, mature_seq, is_5p, is_3p, name=None):
    
    #TODO: MOVE
    pos = miRNA.hairpin.find(mature_seq) # mature seq. must be in hairpin.
    
    if pos >= 0:
        endpos = pos + len(mature_seq)
        if is_5p:
            miRNA.set_five_prime(pos, endpos)
        elif is_3p:
            miRNA.set_three_prime(pos, endpos)
        elif pos+5 < len(miRNA.hairpin) / 2: # not given, but aligns with 5' end
            miRNA.set_five_prime(pos, endpos)
        else:
            miRNA.set_three_prime(pos, endpos)
            
    else:
        print "retard hairpin, or missing suffix..."
        print "\thairpin", miRNA.hairpin
        print "\tmature:", mature_seq
        print "\tmature in hairpin?", mature_seq in miRNA.hairpin
        print "\tname?", name

        


def read_miRNA(mature_file, hairpin_file):
    
    micro_rnas = {}
    key_list = []
    
    with open(hairpin_file) as hairpins:
        print "ok"
        hairpin = ""
        is_human = False
        name = ""
        for mature_seq in hairpins:
            mature_seq = mature_seq.lower().strip()
            
            if mature_seq[0] == ">":

                if is_human:
                    micro_rnas[name] = Micro_RNA(hairpin.replace("U", "T"))
                
                if mature_seq[1:4] == "hsa":
#                     print mature_seq[1:4]
                    is_human = True
                    hairpin = ""
                    name = mature_seq.split()[0]
                    print "hairpin \"%s\"" % name, len(name)
                    
            elif is_human:
                hairpin += mature_seq
    
    print len(micro_rnas)

#     return
    for k,v in micro_rnas.iteritems():
        print k, "!", len(k)
        key_list.append(k)

    print len(key_list)
    


    print ">hsa-mir-4747" in micro_rnas
    print ">hsa-miR-17" in micro_rnas, ">hsa-miR-17"
#     return
#     return


    mature_keys = []
    i = 0
    no_hits = 0
    print "finding mature seqs"
    with open(mature_file) as mature_seqs:
        
        is_human = False
        name = ""
        for mature_seq in mature_seqs:
            mature_seq = mature_seq.lower().strip()
            if mature_seq[0] == ">":
                if mature_seq[1:4] == "hsa":
                    is_human = True
                    name = mature_seq.split()[0]
#                     print name
                    
                else:
                    is_human = False
                    
            elif is_human:
#                 print "before %s" %name
                is_3p = False
                is_5p = False
                if name.endswith("3p"):
                    is_3p = True
#                     print "3p"
                    name = name[:-3]
                elif name.endswith("5p"):
                    is_5p = True
#                     print "5p", name, name[:-3]
                    name = name[:-3]
                mature_keys.append(name)
                
#                 print name, len(name)
#                 print "after: %s" %name
                    
#                 print "is \"%s\" in dict?" % name, name in micro_rnas
                if name in micro_rnas:
                    _add_position(micro_rnas[name], mature_seq, is_5p, is_3p, name)

                else:
                    x = 1
                    while name+"-"+str(x) in micro_rnas:
                        _add_position(micro_rnas[name+"-"+str(x)], mature_seq, is_5p, is_3p, name+"-"+str(x))
                        x += 1
                    else:
                        if x == 1:
                            c = "a"
                            while name+c in micro_rnas:
                                _add_position(micro_rnas[name+c], mature_seq, is_5p, is_3p, name+c)
                                c += chr(ord(c) + 1) 
                                
                            else:                        
                                if c == "a":
                                    no_hits += 1
                                    
                                    print "\tno hit", name, x, c, is_3p, is_5p
                    
                    


    print "no hits:", no_hits

    print "testing only"
    
    key_list = set(key_list)
    mature_keys = set(mature_keys)
    
    print
    print "hairpins", len(key_list)
    print "mature seqs", len(mature_keys)
    
    hp_nrs = set( x.split("-")[2] for x in key_list)

    mat_nrs = set( x.split("-")[2] for x in key_list)

    print "hairpins minus mature seqs:", len(key_list - mature_keys), key_list - mature_keys
    print "mature seqs minus hairpins:", len(mature_keys - key_list), mature_keys - key_list
    print "both", len(mature_keys & key_list)
    
    print '>hsa-let-7a-2' in mature_keys
    print '>hsa-let-7a-2' in key_list
    
    
    
    
    print
    print "hairpins:", len(hp_nrs)
    print "only in hairpins", len(hp_nrs - mat_nrs)
    print
    print "mature seqs:", len(mat_nrs)
    print "only in mature seqs:", len(mat_nrs - hp_nrs)

read_miRNA("mature.fa", "hairpin.fa")




