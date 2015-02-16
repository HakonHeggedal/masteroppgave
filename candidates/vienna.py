
from subprocess import PIPE, Popen
import re
import numpy
from multiprocessing import Pool, current_process
import candidates



def energy_fold2(candidates):
    """ concurrent vienna energy+fold+entropy  """
    
    cores = 15
    pool = Pool(cores)
    
    def get_hp(candidate):
        return candidate.hairpin
    
    def get_hp_10(candidate):
        return candidate.hairpin_padded_40[30:-30]

    
    param_iter = map(get_hp, candidates)
    fold_energy = pool.map(_viennafold, param_iter)
    
    for c, fold, en in zip(candidates, fold_energy):
        c.set_fold_hairpin(fold, en)
        
    param_iter = map(get_hp_10, candidates)
    fold_energy_bitpair = pool.map(_wrap_fold_entropy, param_iter)
    
    for c, fold, en, bp in zip(candidates, fold_energy_bitpair):
        c.set_fold_10(fold, en)
        c.set_bitpair_entropy(bp)
    
    print "lol it may have worked..."
    assert False
    

def energy_fold(candidates):
    
    normal_filename = "pool"
    ps_filename = ">" + normal_filename # super retarded?
    derpy_filename = normal_filename + "_dp.ps"
    i = 0

    for candidate in candidates:
        i += 1
        if i % 100 == 0:
            print ".",
             
        fold, energy = _viennafold(candidate.hairpin, do_ps=False)
         
        fold_part = candidate.hairpin_padded_40[30:-30]
        
        fold10, energy10 = _viennafold(fold_part, filename=ps_filename, do_ps=True)
         
        bitpair_probs = _read_bitpair_probs(derpy_filename)
        bitpair_entropy_dict = _shannon_entropy(bitpair_probs)
         
#         fold40, energy40 = _viennafold(candidate.hairpin_padded_40)
#         candidate.set_viennafold(fold, energy, fold10, energy10, fold40, energy40)
         
        candidate.set_fold_hairpin(fold, energy)
        candidate.set_fold_10(fold10, energy10)
        candidate.set_bitpair_entropy(bitpair_entropy_dict)
        
        
def _wrap_fold_entropy(fold_part):
    
    file_name = current_process().name()
    
    fold, energy = _viennafold(fold_part, do_ps=True)

    normal_filename = "pool"
    
    ps_filename = ">" + normal_filename
    read_ps = normal_filename + "_dp.ps"
    
    bitpair_probs = _read_bitpair_probs(read_ps)
    bp_dict = _shannon_entropy(bitpair_probs)

    return fold, energy, bp_dict
    

def _viennafold(sequence, filename="", do_ps=False):
    """runs vienna folding, returns folding (string) and energy (double)"""
    
    cmds = ["RNAfold", "-p"] if do_ps else ["RNAfold"]
    seq_lines = filename+"\n"+sequence if do_ps else sequence
#     print cmds
#     print seq_lines
    
    mypipe = Popen(cmds, stdin=PIPE, stdout=PIPE, bufsize=-1)
    ans, errors = mypipe.communicate(seq_lines)

    assert errors is None


    ans = ans.strip()
    
    match_number = re.search("-?[0-9]+[.][0-9]+", ans)
#     print sequence
#     print match_number.group(0)
    energy = float(match_number.group(0))
    
    
    match_fold = re.search("[.\(\)]{5,}", ans)
#     print match_fold.group(0)
    fold = match_fold.group(0)

#     print fold, energy
    return fold, energy
    

def _read_bitpair_probs(filename):
    """ extract the bp probatilities from the RNAfold dotplot postscript file"""
    
    bitpair_score = {}
    
    with open(filename) as tmp_plot:
        for line in tmp_plot.readlines():
            line.strip()
            line = line.split()
            if(len(line) == 4 and line[3] == "ubox"):
                bit1 = int(line[0])
                bit2 = int(line[1])
                score = float(line[2])
                
                if bit1 not in bitpair_score: bitpair_score[bit1] = {}
                if bit2 not in bitpair_score: bitpair_score[bit2] = {}
                
                bitpair_score[bit1][bit2] = score
                bitpair_score[bit2][bit1] = score
    
    return bitpair_score



def _shannon_entropy(bitpair_probs):
    """ calculate the positional shannon entropy
    ( as described in Huynen, Gutell and Konings (1997) )"""
    entropy = {}
    for nt_id, bp_dict in bitpair_probs.iteritems():
        entropy[nt_id] = - sum([x * numpy.log(x) for x in bp_dict.values()])
    return entropy


# def _refold_rna(self):
# 
# #         self.mirna_graph = fgb.BulgeGraph()
# #         self.mirna_graph.from_dotbracket(self.dotbracket_seq)
#         self._mk_bp_propensities(self.name + "_dp.ps")
#         self._mk_shannon_entropy()
# #         os.remove(self.name + "_ss.ps")
# #         os.remove(self.name + "_dp.ps")




# # testing only:
# seq2 = "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACAUCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU"
#      
# print  
# print _viennafold(seq2, do_ps=True)



