'''
Created on 29. okt. 2014

@author: hakon
'''

# class Micro_RNA:
#     def __init__(self,hairpin):
#         self.hairpin = hairpin
#         self.pos_5_begin = None
#         self.pos_5_end = None
#         self.pos_3_begin = None
#         self.pos_3_end = None
#     
#     def set_five_prime(self,pos_5_begin, pos_5_end):
#         self.pos_5_begin = pos_5_begin
#         self.pos_5_end = pos_5_end
#     
#     def set_three_prime(self, pos_3_begin, pos_3_end):
#         self.pos_3_begin = pos_3_begin
#         self.pos_3_end = pos_3_end
#         

def read_miRNA_fasta(fasta_file):
    ''' reads miRNAs from mirbase file, and returns {id: dna_seq}'''
    
    human_to_seq = {}
    other_to_seq = {}
    
    with open(fasta_file) as miRNAs:
        
        is_human = False
        dna_seq = ""
        miRNAid = ""
        
        for miRNA_part in miRNAs:
            
            if miRNA_part[0] == ">":
                # save previous
                if is_human:
                    assert miRNAid not in human_to_seq
                    human_to_seq[miRNAid] = dna_seq
                else:
                    assert miRNAid not in other_to_seq
                    other_to_seq[miRNAid] = dna_seq
                
                # save info about current dna_seq, empty dna_seq
                is_human = miRNA_part[1:4] == "hsa"    
                miRNAid = miRNA_part.split()[0].lower()
                dna_seq = ""
            
            else:
                miRNA_part = miRNA_part.upper().strip()
                dna_seq += miRNA_part.replace("U", "T")

        if is_human:
            assert miRNAid not in human_to_seq
            human_to_seq[miRNAid] = dna_seq
        else:
            assert miRNAid not in other_to_seq
            other_to_seq[miRNAid] = dna_seq
            
    
    print
    print "last one human?", is_human, dna_seq
    print "human miRNAs", len(human_to_seq)
    print "other miRNAs", len(other_to_seq)        
    
    return human_to_seq, other_to_seq





def _add_mature_pos(miRNAid, miRNA_dict, hairpin_dict, startpos, endpos, is_5p, is_3p):
    
    if miRNAid not in miRNA_dict:
        miRNA_dict[miRNAid] = [-1,-1,-1,-1]
    
    assert startpos >= 0
    assert startpos < endpos
    

    old = miRNA_dict[miRNAid]
    
    if is_5p:
        miRNA_dict[miRNAid] = [startpos, endpos, old[2], old[3]]
        if startpos+5 >= len(hairpin_dict[miRNAid]) / 2:
            print "strange behavior 5p", miRNA_dict[miRNAid], miRNAid, startpos, len(hairpin_dict[miRNAid])
    elif is_3p:
        miRNA_dict[miRNAid] = [old[0], old[1], startpos, endpos]
        if startpos+5 < len(hairpin_dict[miRNAid]) / 2:
            print "strange behavior 3p", miRNA_dict[miRNAid], miRNAid, startpos, len(hairpin_dict[miRNAid])
    elif startpos+5 < len(hairpin_dict[miRNAid]) / 2: # not given, but aligns with 5' end
        miRNA_dict[miRNAid] = [startpos, endpos, old[2], old[3]]
    else:
        miRNA_dict[miRNAid] = [old[0], old[1], startpos, endpos]
    
    if miRNA_dict[miRNAid][1] != -1 and miRNA_dict[miRNAid][2] != -1:
        if miRNA_dict[miRNAid][1] >= miRNA_dict[miRNAid][2]:
            print "!error!", miRNA_dict[miRNAid], miRNAid, is_5p, is_3p, startpos, len(hairpin_dict[miRNAid])
            print "-!!!", startpos+5, len(hairpin_dict[miRNAid]) / 2
#         assert miRNA_dict[miRNAid][1] < miRNA_dict[miRNAid][2]


def combine_hairpin_mature(id_to_hairpin, id_to_mature):
    ''' combines hairpins with corresponding mature seqs into miRNA structs 
        hairpin identifiers and mature seq identifiers do not map directly
        returns {hairpinID: [5p,3p]}
    '''
    
    harpinID_to_mature = {}

    
    for mature_id, mature_seq in id_to_mature.iteritems():
#         print mature_id, mature_seq
                
        found_hairpin = False
        is_3p = False
        is_5p = False
        
        if mature_id.endswith("3p"):
            is_3p = True
            mature_id = mature_id[:-3]
        elif mature_id.endswith("5p"):
            is_5p = True
            mature_id = mature_id[:-3]
        
#                 print mature_id, len(mature_id)
#                 print "after: %s" %mature_id
#                 print "is \"%s\" in dict?" % mature_id, mature_id in id_to_hairpin

        if mature_id in id_to_hairpin:
            startpos = id_to_hairpin[mature_id].find(mature_seq) # mature seq. must be in hairpin.
            if startpos >= 0:
                found_hairpin = True
                endpos = startpos + len(mature_seq)
                _add_mature_pos(mature_id, harpinID_to_mature, id_to_hairpin,
                                startpos, endpos, is_5p, is_3p)
                
        x = 1

        if mature_id+"-"+str(x) not in id_to_hairpin:
            x = 2

        while mature_id+"-"+str(x) in id_to_hairpin:
            startpos = id_to_hairpin[mature_id+"-"+str(x)].find(mature_seq)
            if startpos >= 0:
                found_hairpin = True
                endpos = startpos + len(mature_seq)
                _add_mature_pos(mature_id+"-"+str(x), harpinID_to_mature, id_to_hairpin,
                              startpos, endpos, is_5p, is_3p)
            x += 1

        c = "a"
        while mature_id+c in id_to_hairpin:
            startpos = id_to_hairpin[mature_id+c].find(mature_seq)
            if startpos >= 0:
                found_hairpin = True
                endpos = startpos + len(mature_seq)
                _add_mature_pos(mature_id+c, harpinID_to_mature, id_to_hairpin,
                                startpos, endpos, is_5p, is_3p)
            c += chr(ord(c) + 1)
                             
        if not found_hairpin:
            print "\tno hit", mature_id, x, c, is_3p, is_5p
    
    print "FINISHED assembling miRNAs"
    return harpinID_to_mature








def write_miRNA(mi_rnas, filename):

    with open(filename, "w") as outfile:
        for name, hairpin in mi_rnas.iteritems():
            outfile.write(name + "\n")
            outfile.write(hairpin + "\n")     
            


def similar_hairpins(specie_to_hairpin, others_to_hairpin):
    
    id_to_speciecount = {}
    other_species = {}
    
    for miRNAid, hairpin in others_to_hairpin.iteritems():
        
        if hairpin in other_species:
            other_species[hairpin].add(miRNAid[1:4])
        else:
            other_species[hairpin] = set([miRNAid[1:4]])
    
    for miRNAid, hairpin in specie_to_hairpin.iteritems():
        
        if hairpin in other_species:
            id_to_speciecount[miRNAid] = len(other_species[hairpin])
        else:
            id_to_speciecount[miRNAid] = 0
            
    return id_to_speciecount
        
        
    


# 
# def read_miRNA(mature_file, hairpin_file):
#     
#     miRNAid_to_hairpin = {}
#     
#     with open(hairpin_file) as hairpins:
#         hairpin = ""
#         is_human = False
#         name = ""
#         for hairpin_seq in hairpins:
#             hairpin_seq = hairpin_seq.lower().strip()
#             
#             if hairpin_seq[0] == ">":
# 
#                 if is_human:  # previous entry
#                     miRNAid_to_hairpin[name] = Micro_RNA(hairpin.replace("u", "t"))
#                 
#                 if hairpin_seq[1:4] == "hsa":
# #                     print mature_seq[1:4]
#                     is_human = True
#                     name = hairpin_seq.split()[0]
# #                     print "hairpin \"%s\"" % name, len(name)
#                 else:
#                     is_human = False
#                 hairpin = ""
#                     
#             elif is_human:
#                 hairpin += hairpin_seq
# 
#     no_hits = 0
# #     print "finding mature seqs"
#     with open(mature_file) as mature_seqs:
#         
#         is_human = False
#         name = ""
#         for mature_seq in mature_seqs:
#             mature_seq = mature_seq.lower().strip()
#             mature_seq = mature_seq.replace("u", "t")
#             if mature_seq[0] == ">":
#                 if mature_seq[1:4] == "hsa":
#                     is_human = True
#                     name = mature_seq.split()[0]
# #                     print name
#                     
#                 else:
#                     is_human = False
#                     
#             elif is_human:
#                 found_hairpin = False
#                 is_3p = False
#                 is_5p = False
#                 
#                 if name.endswith("3p"):
#                     is_3p = True
#                     name = name[:-3]
#                 elif name.endswith("5p"):
#                     is_5p = True
#                     name = name[:-3]
#                 
# #                 print name, len(name)
# #                 print "after: %s" %name
# #                 print "is \"%s\" in dict?" % name, name in miRNAid_to_hairpin
#                 if name in miRNAid_to_hairpin:
#                     pos = miRNAid_to_hairpin[name].hairpin.find(mature_seq) # mature seq. must be in hairpin.
#                     if pos >= 0:
#                         found_hairpin = True
#                         endpos = pos + len(mature_seq)
#                         _add_position(miRNAid_to_hairpin[name], pos, endpos, is_5p, is_3p)
# 
#                 x = 1
#                 if name+"-"+str(x) not in miRNAid_to_hairpin:
#                     x = 2
#                 
#                 while name+"-"+str(x) in miRNAid_to_hairpin:
#                     pos = miRNAid_to_hairpin[name+"-"+str(x)].hairpin.find(mature_seq)
#                     if pos >= 0:
#                         found_hairpin = True
#                         endpos = pos + len(mature_seq)
#                         _add_position(miRNAid_to_hairpin[name+"-"+str(x)],
#                                       pos, endpos, is_5p, is_3p)
#                     x += 1
# 
#                 c = "a"
#                 while name+c in miRNAid_to_hairpin:
#                     pos = miRNAid_to_hairpin[name+c].hairpin.find(mature_seq)
#                     if pos >= 0:
#                         found_hairpin = True
#                         endpos = pos + len(mature_seq)
#                         _add_position(miRNAid_to_hairpin[name+c],pos, endpos, is_5p, is_3p)
#                     c += chr(ord(c) + 1)
#                                      
#                 if not found_hairpin:
#                     no_hits += 1
#                     
#                     print "\tno hit", name, x, c, is_3p, is_5p
#     
#     print "finished assembling miRNAs"
#     return miRNAid_to_hairpin
#     




    

# read_miRNA("mature.fa", "hairpin.fa")



# def mirna_copies(human_mirna, hairpin_file):
#     
#     human_hairpins = [x.hairpin.upper().strip() for x in human_mirna ]
#     
#     miRNA_to_dupl = {}
#     for h in human_hairpins:
#         miRNA_to_dupl[h] = 0
# 
#     
#     print len(human_hairpins), human_hairpins[0]
#     print len(set(human_hairpins))
#     for h in human_hairpins:
#         print "human hps:",h, len(h)
#         
#     hh = set(human_hairpins)
#     ah = set()
#     
#     with open(hairpin_file) as hairpins:
#         name = ""
#         for hairpin_part in hairpins:
#             hairpin_part = hairpin_part.strip()
#             if hairpin_part[0] == ">":
#                 if hairpin_part[1:4] == "hsa":
#                     print "human:", name
#                     ah.add(name)    
# #                 print "all hps", name, len(name)
#                 if name in miRNA_to_dupl:
#                     print "\tduplicate", hairpin_part
#                     miRNA_to_dupl[name] += 1
#                 name = ""
#             else:
#                 name += hairpin_part.replace("U", "T")
#     
#     print "-----"
# #     hh = [x for x in hh]
# #     ah = [x for x in ah]
# #     
# #     x = 0
# #     for h in sorted(hh):
# #         print h
# #         x+= 1
# #         if x == 10:
# #             break
# #         
# #     x = 0
# #     for h in sorted(ah):
# #         print h
# #         x+= 1
# #         if x == 10:
# #             break
# 
#     print "human:", len(hh)
#     print "all", len(ah)
#     print "both", len(ah & hh)
#     print "distinct", len(ah ^ hh)
#     assert False
#     
#     return miRNA_to_dupl
#                 
#     

