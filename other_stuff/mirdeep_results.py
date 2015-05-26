'''
Created on 22. mai 2015

@author: hakon
'''
from sklearn import metrics
from matplotlib import pyplot as plot
import pickle


# novel:gi|224589821|ref|NC_000009.11|_21233

def _get_mirdeep_results():
    filename = "result_21_05_2015_t_22_02_09.bed"
    
    results = []
    with open(filename) as mirDeep_results:
        for line in mirDeep_results.readlines():
            
            parts = line.split()
    #         print parts
            
            identifier = parts[3]
            
            novel_type = identifier[:5] == "novel"
    
            
            if not novel_type:
                break
                print novel_type, identifier[:5]
                
            genome = identifier.split("|")[3]
            print genome,
            
            score = float(parts[4])
            strand_dir = parts[5]
            
            startpos = int(parts[6])
            endpos = int(parts[7])
            
            print score,
            print strand_dir,
            print startpos,
            print endpos
            results.append([score, genome, strand_dir, startpos, endpos])
            
    #         break
    
    print len(results)
    return results


def mirdeep_make_roc_data(candidate_tree, candidate_to_miRNA, miRNA_high_conf):
    
    mirdeep_results = _get_mirdeep_results()
    
    
    def _is_hc(c):
        hashval = c.data.chromosome + c.data.chromosome_direction + str(c.data.pos_5p_begin)
        if hashval in candidate_to_miRNA:
            if candidate_to_miRNA[hashval] in miRNA_high_conf:
                return True
        return False
    
    def _get_mirna_id(c):
        hashval = c.data.chromosome + c.data.chromosome_direction + str(c.data.pos_5p_begin)
        if hashval in candidate_to_miRNA:
            return candidate_to_miRNA[hashval]
        return None
    
    mirna_names = []
    hc_10 = []
    scores = []
    
    mirdeep_candidates = {}
    
    for [score, chromosome, strand_dir, startpos, endpos] in mirdeep_results:
        
        scores.append(score)
        

        candidates = candidate_tree[chromosome][startpos:endpos]
        
        candidates = [c for c in candidates if c.data.chromosome_direction == strand_dir]
        
        
        is_hc = False
        name = None

        for c in candidates:
            currentname = _get_mirna_id(c)
            name = currentname if currentname else name
            
            if _is_hc(c):
                is_hc = True
                break
        
        if name == None:
            hashval = c.data.chromosome + c.data.chromosome_direction + str(c.data.pos_5p_begin)
            if len(candidates):
                mirdeep_candidates[hashval] = candidates[0] 
            print len(candidates)
        
        mirna_names.append(name)
        hc_10.append(1 if is_hc else 0)
        
    
    print len(mirna_names), len(set(mirna_names))
    print len(hc_10)
    print len(scores)
    
    print mirna_names
    print hc_10
    print scores
    
    return mirdeep_candidates



names = ['>hsa-mir-21', '>hsa-let-7a-1', '>hsa-let-7a-3', '>hsa-let-7a-2', '>hsa-let-7f-2', None, '>hsa-mir-127', '>hsa-mir-150', None, None, '>hsa-mir-205', '>hsa-let-7b', '>hsa-mir-26b', '>hsa-mir-146b', None, None, '>hsa-mir-130a', '>hsa-let-7d', '>hsa-let-7c', '>hsa-mir-486-2', '>hsa-mir-486-1', '>hsa-mir-320a', None, '>hsa-mir-93', '>hsa-mir-23b', '>hsa-mir-107', '>hsa-mir-24-1', '>hsa-mir-410', '>hsa-mir-130b', '>hsa-mir-181c', '>hsa-mir-29c', None, '>hsa-mir-204', None, None, '>hsa-mir-106b', '>hsa-mir-501', '>hsa-mir-214', None, '>hsa-mir-372', '>hsa-mir-215', '>hsa-mir-889', '>hsa-mir-421', None, '>hsa-mir-373', None, '>hsa-mir-582', '>hsa-mir-452', '>hsa-mir-512-2', None, '>hsa-mir-502', '>hsa-mir-432', '>hsa-mir-1323', '>hsa-mir-429', None, '>hsa-mir-652', '>hsa-mir-185', '>hsa-mir-520a', None, '>hsa-mir-542', None, '>hsa-mir-324', '>hsa-mir-515-2', '>hsa-mir-515-1', '>hsa-mir-122', None, '>hsa-mir-518b', '>hsa-mir-664a', '>hsa-mir-206', '>hsa-mir-517a', '>hsa-mir-517b', None, '>hsa-mir-323b', '>hsa-mir-212', None, None, None, '>hsa-mir-518a-2', '>hsa-mir-184', '>hsa-mir-518a-1', '>hsa-mir-656', '>hsa-mir-487a', '>hsa-mir-2355', None, '>hsa-mir-376b', '>hsa-mir-1180', '>hsa-mir-380', '>hsa-mir-3909', None, '>hsa-mir-511', '>hsa-mir-944', '>hsa-mir-137', '>hsa-mir-550a-2', '>hsa-mir-550a-1', '>hsa-mir-2110', '>hsa-mir-4677', '>hsa-mir-216a', None, None, '>hsa-mir-549a', '>hsa-mir-518f', None, None, '>hsa-mir-598', None, None, '>hsa-mir-2682', '>hsa-mir-518c', '>hsa-mir-655', '>hsa-mir-550a-3', '>hsa-mir-2277', '>hsa-mir-20b', None, '>hsa-mir-217', None, None, '>hsa-mir-106a', None, '>hsa-mir-153-1', '>hsa-mir-3200', None, None, '>hsa-mir-520d', '>hsa-mir-4772', '>hsa-mir-1283-2', None, None, None, None, None, '>hsa-mir-519d', None, '>hsa-mir-3152', '>hsa-mir-517c', '>hsa-mir-371b', None, '>hsa-mir-153-2', None, '>hsa-mir-2114', '>hsa-mir-525', '>hsa-mir-4662a', '>hsa-mir-524', '>hsa-mir-449a', '>hsa-mir-3677', '>hsa-mir-1283-1', '>hsa-mir-627', None, '>hsa-mir-1-1', '>hsa-mir-526b', '>hsa-mir-1303', '>hsa-mir-541', None, None, '>hsa-mir-3912', None, '>hsa-mir-659', '>hsa-mir-1255a', None, None, '>hsa-mir-4433b', None, '>hsa-mir-1976', None, None, None, None, None, None, '>hsa-mir-520b', '>hsa-mir-585', None, None, None, None, None, None, '>hsa-mir-3940', None, None, None, None, None, '>hsa-mir-551b', None, None, None, None, '>hsa-mir-95', None, None, '>hsa-mir-2115', None, None, None, '>hsa-mir-653', None, None, '>hsa-mir-605', None, None, None, '>hsa-mir-216b', None, '>hsa-mir-580', '>hsa-mir-3155a', '>hsa-mir-4424', '>hsa-mir-760', '>hsa-mir-3174', '>hsa-mir-1245a', '>hsa-mir-302d', '>hsa-mir-551a', None, None, None, None, None, None, None, None, '>hsa-mir-4661', None, None, '>hsa-mir-3194', None, '>hsa-mir-6513', '>hsa-mir-3942', None, None, '>hsa-mir-3157', '>hsa-mir-4423', None, None, None, '>hsa-mir-3934', None, None, '>hsa-mir-3138', None, None, None, None, '>hsa-mir-1537', '>hsa-mir-643', '>hsa-mir-4745', '>hsa-mir-3177', None, None, None, None, None, '>hsa-mir-3919', None, '>hsa-mir-3187', None, None, None, None, None, None, None, '>hsa-mir-7854', '>hsa-mir-3127', '>hsa-mir-597', None, '>hsa-mir-3126', None, None, None, None, '>hsa-mir-3145', None, None, None, '>hsa-mir-1913', '>hsa-mir-3150b', None, None, None, None, '>hsa-mir-5579', None, '>hsa-mir-548at', None, None, None, '>hsa-mir-676', None, None, None, None, '>hsa-mir-4773-1', None, '>hsa-mir-4773-2', None, None, None, None, '>hsa-mir-4671', None, '>hsa-mir-4520-1', None, None, None, None, '>hsa-mir-3667', None, None, None, None, None, None, None, '>hsa-mir-548ag-2', None, None, '>hsa-mir-6501', None, None, None, None, None, None, '>hsa-mir-4482', None, '>hsa-mir-4753', None, '>hsa-mir-4766', None, None, None, None, None, None, '>hsa-mir-3622a', None, None, None, None, None, None, None, '>hsa-mir-4802', None, '>hsa-mir-3918', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-4731', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-4714', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-4467', None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-1245b', None, '>hsa-mir-4422', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-5583-1', '>hsa-mir-5583-2', None, None, None, '>hsa-mir-4433a', None, None, None, None, None, None, None, None, '>hsa-mir-3120', None, None, None, None, None, None, None, '>hsa-mir-520e', None, None, '>hsa-mir-610', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-4699', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-7850', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-4501', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-3612', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-641', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-934', None, None, None, None, None, None, None, None, None, None, '>hsa-mir-5189', None, None, None, None, None, None, None, None, '>hsa-mir-301b', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-665', '>hsa-mir-3620', None, None, None, '>hsa-mir-455', None, None, None, None, None, None, '>hsa-mir-552', None, None, None, None, None, None, None, None, '>hsa-mir-224', None, None, None, None, None, '>hsa-mir-371a', None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-142', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-527', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-145', None, None, None, None, None, '>hsa-mir-1305', '>hsa-mir-520g', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-302b', None, None, None, None, None, None, None, None, None, '>hsa-mir-1286', None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-520h', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-521-1', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-301a', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-665', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-1292', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-6837', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-486-1', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-let-7a-3', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-2277', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-6842', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, '>hsa-mir-3675', None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]
classes = [0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
scores = [9537144.4, 5835593.2, 5830523.2, 5827786.1, 3685170.5, 3614410.4, 2131644.3, 2042873.8, 1577541.3, 995005.9, 752860.7, 649438.5, 493753.6, 462758.0, 411086.4, 330941.5, 261502.0, 206552.7, 174341.2, 162542.6, 160565.0, 158173.0, 145121.1, 129440.7, 128566.8, 126193.8, 124811.1, 117085.6, 103743.5, 94966.2, 90735.2, 81919.4, 77833.7, 75258.4, 69068.3, 56949.8, 44338.4, 35175.0, 35120.6, 33284.9, 32787.9, 27702.2, 24719.4, 22650.4, 19679.3, 19376.5, 18566.2, 17469.1, 17446.3, 17446.3, 17380.8, 16832.9, 14743.4, 14102.1, 13305.0, 11640.6, 11406.1, 10656.5, 8751.0, 8137.0, 7813.7, 7480.5, 6625.9, 6625.9, 6608.8, 6563.1, 6422.5, 6207.5, 6190.8, 5978.0, 5977.4, 5406.5, 5242.3, 5227.5, 5038.6, 4633.5, 4556.0, 4437.6, 4315.2, 4286.7, 4225.5, 4116.1, 3974.5, 3793.5, 3774.4, 3767.8, 3675.4, 3658.6, 3318.8, 3318.8, 3075.8, 3056.0, 2937.9, 2937.1, 2893.8, 2815.3, 2770.7, 2756.5, 2717.6, 2699.3, 2593.3, 2544.3, 2340.3, 2199.9, 2176.4, 2092.9, 2082.1, 1992.9, 1960.5, 1823.1, 1815.8, 1783.0, 1723.0, 1677.6, 1655.2, 1653.9, 1589.3, 1541.2, 1530.6, 1512.3, 1489.7, 1401.4, 1368.9, 1356.7, 1336.1, 1291.1, 1274.7, 1274.7, 1274.6, 1274.6, 1238.7, 1213.3, 1204.9, 1201.0, 1200.1, 1163.9, 1097.5, 1090.0, 1076.1, 1075.1, 985.8, 975.9, 895.3, 863.4, 863.4, 823.1, 809.6, 776.7, 768.5, 767.8, 696.3, 662.9, 642.7, 637.8, 637.3, 636.6, 627.5, 624.8, 593.2, 590.1, 581.3, 564.6, 555.1, 553.2, 546.5, 526.0, 478.5, 477.1, 470.3, 466.9, 462.4, 462.4, 462.4, 462.4, 462.3, 461.0, 456.8, 455.9, 455.7, 455.7, 455.7, 455.7, 431.6, 429.0, 428.9, 420.3, 418.2, 416.1, 411.1, 411.1, 410.3, 406.1, 383.8, 366.1, 365.2, 363.6, 363.5, 362.3, 355.2, 355.2, 349.4, 346.9, 341.6, 310.6, 308.9, 303.6, 303.4, 299.6, 287.0, 283.3, 278.2, 276.8, 276.7, 274.8, 274.2, 272.4, 268.9, 263.9, 262.9, 261.3, 261.0, 260.9, 249.0, 244.2, 242.9, 241.5, 241.0, 236.7, 236.5, 231.4, 227.8, 224.9, 224.4, 216.5, 207.2, 206.4, 204.5, 200.5, 199.3, 196.1, 195.4, 189.7, 186.5, 186.4, 186.1, 185.0, 177.6, 175.4, 175.2, 175.1, 174.3, 172.6, 169.7, 166.6, 165.6, 160.2, 159.8, 157.1, 157.0, 152.7, 152.4, 144.5, 141.8, 136.8, 135.4, 135.4, 130.7, 129.9, 129.7, 128.6, 128.3, 125.0, 120.3, 119.9, 119.5, 119.0, 117.3, 117.2, 117.0, 115.9, 115.4, 113.8, 113.6, 113.3, 112.5, 109.1, 108.4, 105.9, 104.9, 104.1, 102.7, 102.2, 102.2, 100.6, 99.3, 99.3, 99.1, 98.8, 98.5, 98.4, 97.8, 97.1, 97.1, 96.8, 96.7, 96.6, 93.9, 93.6, 92.9, 92.7, 92.3, 91.1, 90.8, 89.1, 88.9, 88.6, 86.8, 86.5, 86.5, 86.3, 84.9, 84.8, 84.5, 82.6, 82.2, 80.0, 79.9, 79.5, 78.0, 77.4, 77.2, 77.2, 77.2, 76.7, 76.4, 75.6, 74.6, 73.9, 73.8, 73.7, 73.2, 71.6, 71.6, 71.5, 71.2, 70.8, 70.0, 69.4, 68.5, 68.4, 67.7, 67.4, 67.4, 66.5, 65.8, 65.5, 65.3, 65.0, 64.6, 64.4, 64.2, 64.2, 63.8, 63.8, 63.7, 63.0, 62.5, 62.2, 62.2, 62.1, 62.1, 61.7, 61.7, 60.6, 59.9, 59.6, 59.4, 58.3, 58.2, 58.2, 57.8, 57.6, 56.9, 56.6, 56.2, 56.2, 54.9, 54.7, 54.6, 53.8, 52.8, 52.3, 52.1, 52.0, 51.8, 51.0, 50.4, 50.2, 49.8, 49.8, 49.7, 49.5, 48.7, 48.6, 48.5, 48.4, 47.8, 47.7, 47.4, 47.2, 47.1, 46.1, 46.0, 45.8, 45.8, 45.4, 45.4, 45.1, 44.4, 44.1, 43.8, 43.8, 43.6, 43.2, 42.2, 41.9, 41.9, 41.6, 41.6, 41.6, 41.3, 41.2, 40.8, 40.7, 40.3, 39.8, 39.1, 39.1, 38.4, 38.4, 38.0, 37.6, 37.0, 36.9, 36.7, 36.7, 36.3, 36.1, 35.9, 35.8, 35.8, 35.7, 35.3, 35.3, 35.3, 35.2, 35.0, 35.0, 35.0, 34.8, 34.5, 34.2, 34.1, 34.1, 34.0, 33.9, 33.8, 33.3, 33.2, 32.6, 32.1, 32.0, 31.9, 31.8, 31.8, 31.7, 31.6, 31.6, 31.4, 31.3, 31.0, 30.6, 30.5, 29.7, 29.5, 29.3, 29.2, 29.1, 28.8, 28.8, 28.6, 28.6, 28.5, 28.3, 28.1, 28.1, 28.1, 28.0, 28.0, 28.0, 27.9, 27.8, 27.8, 27.8, 27.6, 27.5, 27.0, 26.9, 26.9, 26.8, 26.7, 26.5, 26.3, 26.2, 26.0, 26.0, 25.9, 25.8, 25.8, 25.8, 25.7, 25.2, 24.9, 24.9, 24.9, 24.9, 24.8, 24.6, 24.4, 24.3, 24.2, 24.0, 23.8, 23.8, 23.8, 23.6, 23.5, 23.1, 23.1, 22.8, 22.8, 22.7, 22.6, 22.6, 22.4, 22.3, 22.1, 22.0, 21.9, 21.7, 21.5, 21.5, 21.4, 21.4, 21.4, 21.3, 21.1, 21.1, 21.1, 21.0, 21.0, 20.9, 20.7, 20.6, 20.6, 20.5, 20.3, 20.1, 20.1, 19.9, 19.8, 19.7, 19.6, 19.4, 19.4, 19.3, 19.2, 19.2, 19.2, 19.0, 19.0, 18.9, 18.8, 18.8, 18.7, 18.6, 18.5, 18.5, 18.4, 18.4, 18.4, 18.3, 18.0, 18.0, 18.0, 18.0, 17.9, 17.8, 17.7, 17.6, 17.6, 17.6, 17.5, 17.5, 17.5, 17.5, 17.5, 17.4, 17.4, 17.3, 17.3, 17.3, 17.2, 17.1, 16.9, 16.9, 16.8, 16.8, 16.7, 16.6, 16.6, 16.6, 16.4, 16.4, 16.4, 16.3, 16.3, 16.3, 16.3, 16.2, 16.2, 16.1, 16.0, 15.9, 15.8, 15.8, 15.8, 15.6, 15.6, 15.5, 15.5, 15.4, 15.4, 15.3, 15.2, 15.2, 15.0, 14.8, 14.7, 14.7, 14.7, 14.7, 14.7, 14.7, 14.7, 14.6, 14.5, 14.4, 14.3, 14.3, 14.1, 14.1, 14.1, 13.8, 13.7, 13.7, 13.6, 13.6, 13.6, 13.5, 13.5, 13.5, 13.4, 13.4, 13.3, 13.3, 13.2, 13.1, 13.0, 12.9, 12.8, 12.8, 12.7, 12.6, 12.5, 12.5, 12.5, 12.4, 12.3, 12.3, 12.3, 12.2, 12.2, 12.2, 12.1, 12.1, 12.1, 12.0, 12.0, 12.0, 11.9, 11.8, 11.7, 11.7, 11.6, 11.6, 11.6, 11.4, 11.3, 11.3, 11.1, 11.1, 11.0, 11.0, 10.9, 10.9, 10.9, 10.9, 10.8, 10.5, 10.4, 10.3, 10.3, 10.2, 10.2, 10.2, 10.2, 10.0, 10.0, 10.0, 9.9, 9.9, 9.9, 9.8, 9.7, 9.4, 9.4, 9.4, 9.3, 9.3, 9.2, 9.1, 9.0, 8.9, 8.8, 8.7, 8.5, 8.4, 8.4, 8.2, 7.9, 7.8, 7.5, 7.4, 7.4, 7.1, 7.0, 7.0, 7.0, 6.6, 6.3, 6.3, 6.3, 6.2, 6.2, 6.1, 6.1, 6.1, 6.1, 6.1, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.1, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.4, 4.4, 4.4, 4.4, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.8, 3.8, 3.8, 3.8, 3.8, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6, 3.5, 3.5, 3.5, 3.4, 3.4, 3.4, 3.4, 3.4, 3.3, 3.3, 3.3, 3.3, 3.3, 3.2, 3.2, 3.2, 3.2, 3.2, 3.1, 3.1, 3.1, 3.1, 3.1, 3.0, 3.0, 3.0, 3.0, 3.0, 2.9, 2.9, 2.9, 2.9, 2.8, 2.8, 2.8, 2.8, 2.7, 2.7, 2.7, 2.6, 2.6, 2.6, 2.6, 2.6, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


small_reads = pickle.load( open("names_data.p", "rb"))




print "new:"
print len(filter(None, names))

# 
# print map(len, small_reads)
# print small_reads[0]
# # small_reads = set(small_reads[0])
# print len(small_reads), small_reads

ns = []
cs = []
ss = []
for n, c, s in zip(names, classes, scores):
    if (n and c) or not n:
        ns.append(n)
        cs.append(c)
        ss.append(s)

nameset = set(ns)

found_by_mirDeep = [1 if n in nameset else 0 for n in small_reads[0]]

pickle.dump(found_by_mirDeep, open("found_by_mirDeep.p", "wb"))

pickle.dump(cs, open("classes_mirdeep.p", "wb"))
pickle.dump(ss, open("scores_mirdeep.p", "wb"))



print sum(found_by_mirDeep), len(found_by_mirDeep)

# print 
# print "STATS:"
# print len(small_reads)
# print len(nameset)
# print len(small_reads - nameset)
# print len(nameset - small_reads)
# print len(nameset.union(small_reads))
# print small_reads.issuperset(nameset)

print 
print len(names), len(ns)

print len(filter(None, names))
print sum(classes)

# fpr, tpr, _thresholds = metrics.roc_curve(classes, scores)
fpr, tpr, _thresholds = metrics.roc_curve(cs, ss)
roc_auc = metrics.auc(fpr, tpr)

print 
print "ROC PLOT:"
# print fpr, tpr
print "area under curve:", roc_auc

title = "ROC plot mirdeep"

# 
# 
# plot.plot(fpr, tpr)
# plot.title(title)
# plot.xticks([0.1*x for x in range(0, 11)])
# plot.xlabel("False positive rate")
# plot.yticks([0.1*x for x in range(0, 11)])
# plot.ylabel("True positive rate")
# plot.show()
# 
# 
# # fpr, tpr, _thresholds = metrics.roc_curve(cs, ss)
# fpr, tpr, _thresholds = metrics.roc_curve(classes, scores)
# roc_auc = metrics.auc(fpr, tpr)
# 
# print 
# print "ROC PLOT:"
# # print fpr, tpr
# print "area under curve:", roc_auc
# 
# title = "ROC plot mirdeep"
# 
# 
# 
# plot.plot(fpr, tpr)
# plot.title(title)
# plot.xticks([0.1*x for x in range(0, 11)])
# plot.xlabel("False positive rate")
# plot.yticks([0.1*x for x in range(0, 11)])
# plot.ylabel("True positive rate")
# plot.show()



