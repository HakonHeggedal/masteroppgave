
from subprocess import PIPE, Popen
import re


def energy_fold(candidates):
    
    for candidate in candidates:
        print
        fold, energy = _viennafold(candidate.hairpin)
        candidate.set_viennafold(fold, energy)

        fold, energy = _viennafold(candidate.hairpin_padded)
    assert False


def _viennafold(sequence):
    #TODO: ALL AT ONCE
    """runs vienna folding, returns folding (string) and energy (double)"""
    mypipe  = Popen("RNAfold", stdin=PIPE, stdout=PIPE, bufsize=-1)
    ans, errors = mypipe.communicate(sequence)
    assert errors is None
    
    ans = ans.strip()
    
    match_number = re.search("-?[0-9]+[.][0-9]+", ans)
    print sequence
#     print match_number.group(0)
    energy = float(match_number.group(0))
    
    
    match_fold = re.search("[.\(\)]{5,}", ans)
    print match_fold.group(0)
    fold = match_fold.group(0)

    
    return fold, energy
    
    


# seq2 = "uuggauguuggccuaguucuguguggaagacuagugauuuuguuguuuuuagauaacuaaaucgacaacaaaucacagucugccauauggcacaggccaugccucuacag"
#  
#  
# _viennafold(seq2)