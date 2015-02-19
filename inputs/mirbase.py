'''
Created on 29. okt. 2014

@author: hakon
'''
from candidates.interval_tree_search import MIN_HAIRPIN_LOOP


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
            
    
#     print
#     print "last one human?", is_human, dna_seq
    print "\thuman miRNAs", len(human_to_seq)
    print "\tother miRNAs", len(other_to_seq)        
    
    return human_to_seq, other_to_seq





def _add_mature_pos(miRNAid, miRNA_dict, hairpin_dict, startpos, endpos, is_5p, is_3p):
    
    if miRNAid not in miRNA_dict:
        miRNA_dict[miRNAid] = [-1,-1,-1,-1]
    
    assert startpos >= 0
    assert startpos < endpos
    
    is_added = False
#     Todo: a<b<c<d

    mipos = miRNA_dict[miRNAid]
    
#     print mipos
    if is_5p:
        is_added = True
        miRNA_dict[miRNAid] = [startpos, endpos, mipos[2], mipos[3]]
#         if startpos+5 >= len(hairpin_dict[miRNAid]) / 2:
#             print "strange behavior 5p", miRNA_dict[miRNAid], miRNAid, startpos, len(hairpin_dict[miRNAid])
    elif is_3p:
        is_added = True
        miRNA_dict[miRNAid] = [mipos[0], mipos[1], startpos, endpos]
#         if startpos+5 < len(hairpin_dict[miRNAid]) / 2:
#             print "strange behavior 3p", miRNA_dict[miRNAid], miRNAid, startpos, len(hairpin_dict[miRNAid])
    
    elif startpos+5 < len(hairpin_dict[miRNAid]) / 2: # not given, but starts in 5' end
        if mipos[0] == -1 and (endpos + MIN_HAIRPIN_LOOP < mipos[2] or mipos[2] == -1):
            is_added = True
            miRNA_dict[miRNAid] = [startpos, endpos, mipos[2], mipos[3]]
    elif mipos[2] == -1 and (mipos[1] + MIN_HAIRPIN_LOOP < startpos or mipos[1] == -1):
        is_added = True
        miRNA_dict[miRNAid] = [mipos[0], mipos[1], startpos, endpos]
        
#     else:
#         print "no match---"
#         print miRNA_dict[miRNAid], miRNAid, is_5p, is_3p, startpos, len(hairpin_dict[miRNAid])
    
#     if not is_added:
#         print "\tno match:",miRNA_dict[miRNAid], miRNAid, is_5p, is_3p, startpos, len(hairpin_dict[miRNAid])

    
#     if miRNA_dict[miRNAid][1] != -1 and miRNA_dict[miRNAid][2] != -1:
#         if miRNA_dict[miRNAid][1] >= miRNA_dict[miRNAid][2]:
#             print "\tno hairpin, or overlapping!", miRNA_dict[miRNAid], miRNAid, is_5p, is_3p, startpos, len(hairpin_dict[miRNAid])
#             print "-!!!", startpos+5, len(hairpin_dict[miRNAid]) / 2
#             assert False
#     if ">hsa-mir-4456" in miRNA_dict:
#         print "!error!", miRNA_dict[miRNAid], miRNAid, is_5p, is_3p, startpos, len(hairpin_dict[miRNAid])
#         assert False
#     print mipos
    
    assert mipos[0] == -1 or miRNA_dict[miRNAid][0] < miRNA_dict[miRNAid][1]
    assert mipos[1] == -1 or mipos[2] == -1 or mipos[1] < mipos[2]
    assert mipos[2] == -1 or mipos[2] < mipos[3]
    
#     assert miRNA_dict[miRNAid][0] < miRNA_dict[miRNAid][1] < miRNA_dict[miRNAid][2] < miRNA_dict[miRNAid][3]
#         assert miRNA_dict[miRNAid][1] < miRNA_dict[miRNAid][2]


def combine_hairpin_mature(id_to_hairpin, id_to_mature):
    ''' combines hairpins with corresponding mature seqs into miRNA structs 
        hairpin identifiers and mature seq identifiers do not map directly
        returns {hairpinID: [5p,3p]}
    '''
    
    harpinID_to_mature = {}
    
    print "\n\t combining hairpin and mature seq miRNA..."
    
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
                             
#         if not found_hairpin:
#             print "\tno hit", mature_id, x, c, is_3p, is_5p
    
    
#     print ">hsa-mir-4456" in harpinID_to_mature
#     print harpinID_to_mature[">hsa-mir-4456"]
#     assert False
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
#             if len(other_species[hairpin]) > 5:
#                 id_to_speciecount[miRNAid] = 5
#             else:
# #                 id_to_speciecount[miRNAid] = 1
#                 id_to_speciecount[miRNAid] = len(other_species[hairpin])
        else:
#             print "NO SPECIES"
            id_to_speciecount[miRNAid] = 0
            
    return id_to_speciecount
        
        
















