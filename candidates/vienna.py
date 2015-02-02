
from subprocess import PIPE, Popen
import re


def energy_fold(candidates):
    
    i = 0
    for candidate in candidates:
        i += 1
        if i % 100 == 0:
            print ".",

        fold, energy = _viennafold(candidate.hairpin)
#         fold40, energy40 = _viennafold(candidate.hairpin_padded_40)
        fold10, energy10 = _viennafold(candidate.hairpin_padded_40[30:-30])

#         candidate.set_viennafold(fold, energy, fold10, energy10, fold40, energy40)
        candidate.set_viennafold(fold, energy, fold10, energy10, 0, 0)



def _viennafold(sequence):
    #TODO: ALL AT ONCE
    """runs vienna folding, returns folding (string) and energy (double)"""
    mypipe  = Popen("RNAfold", stdin=PIPE, stdout=PIPE, bufsize=-1)
    ans, errors = mypipe.communicate(sequence)
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
    
    


seq2 = "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACAUCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU"
  
print  
print _viennafold(seq2)