'''
Created on 20. mai 2015

@author: hakon
'''
import itertools
from itertools import repeat, cycle




def fix_miRNA_training_test(annotated_data, annotations, low_confidence_data, hsa_to_hairpin, harpinID_to_matseqs, candidate_to_miRNA):
    
    print
    print "making mirna dataset split"
    
    test_data = annotated_data[0]
    training_data = list(itertools.chain.from_iterable(annotated_data[1:]))

#     def _is_miRNA(c):
#         hashval = c.chromosome+c.chromosome_direction+str(c.hairpin_start)
#         return hashval in candidate_to_miRNA



    # get the miRNA names:
    def _hairpin_ids(dataset):
        
        ids = []
        for c in dataset:
            hashval = c.chromosome+c.chromosome_direction+str(c.hairpin_start)
            
            if hashval in candidate_to_miRNA:
                ids.append(candidate_to_miRNA[hashval])
                
        return  ids
                
    
    test_ids = _hairpin_ids(test_data)
    training_ids = _hairpin_ids(training_data)
    
    print len(test_ids)
    print len(training_ids)
    
    
#     get the hairpins and mature seqs for the given mirnas

    
    
    training_hairpins = set([hsa_to_hairpin[x] for x in training_ids])
    test_hairpins = set([hsa_to_hairpin[x] for x in test_ids])
    
    
    print len(training_hairpins), list(training_hairpins)[0]
    print len(test_hairpins), list(test_hairpins)[0]
    
    training_matures = set()
    test_matures = set()
    
    
    for miRNAid in training_ids:
        for mature in harpinID_to_matseqs[miRNAid]:
            training_matures.add(mature)
            
    for miRNAid in test_ids:
        for mature in harpinID_to_matseqs[miRNAid]:
            test_matures.add(mature)
    
    print len(training_matures), training_matures
    print len(test_matures), test_matures
    print "^^^^^^^^^^^^^^^^^^^^^^"
    # filter and write to files
    
    
    with open("human_hairpins.fa") as hairpin_in:
 
        all_lines = hairpin_in.readlines()
        hp_ids = all_lines[0::2]
        hairpins = all_lines[1::2]
        
        print len(all_lines)
        print len(hp_ids), hp_ids[0]
        print len(hairpins), hairpins[0]
         
        with open("human_hairpins_training.fa", "w") as hairpin_training_out:
            with open("human_hairpins_testing.fa", "w") as hairpin_testing_out:
                 
                for hp_id, hairpin in zip(hp_ids, hairpins):
                    
                    hp_id = hp_id.strip()
                    hairpin = hairpin.strip()
                    
                    if hairpin in training_hairpins:
                        hairpin_training_out.write(hp_id + "\n")
                        hairpin_training_out.write(hairpin + "\n")
                         
                    elif hairpin in test_hairpins:
                        print "."
                        hairpin_testing_out.write(hp_id + "\n")
                        hairpin_testing_out.write(hairpin + "\n")
                    else:
                        pass
#                         print "-"
                        
            


            
    with open("human_matures.fa") as matures_in:

        all_lines = matures_in.readlines()
        mat_ids = all_lines[0::2]
        matures = all_lines[1::2]
        
        print len(all_lines)
        print len(mat_ids)
        print len(matures)
        
        with open("human_matures_training.fa", "w") as matures_training_out:
            with open("human_matures_testing.fa", "w") as matures_testing_out:
                
                for mat_id, mature in zip(mat_ids, matures):
                    
                    mat_id = mat_id.strip()
                    mature = mature.strip()

                    mature = mature.replace("U", "T")
                    
#                     print mature, mature in training_matures, mature in test_matures
                    
                    if mature in training_matures:
                        matures_training_out.write(mat_id + "\n")
                        matures_training_out.write(mature + "\n")
                        
                    if mature in test_matures:
                        matures_testing_out.write(mat_id + "\n")
                        matures_testing_out.write(mature + "\n")
                    else:
                        pass
#                         print "+"
        
        
        
        
        
        
        
        
#             for mat_id, mature_seq in zip(mat_ids, matures):
#                 if matseq in training_matures:
#                     pass        





        