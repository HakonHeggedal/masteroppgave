
from subprocess import PIPE, Popen



def energy_fold(candidates):
    
    for candidate in candidates:
        fold, energy = _viennafold(candidate.seq_hairpin)
        candidate.set_viennafold(fold, energy)




def _viennafold(sequence, ):
    """runs vienna folding, returns folding (string) and energy (double)"""
    mypipe  = Popen("RNAfold", stdin=PIPE, stdout=PIPE, bufsize=-1)
    ans, errors = mypipe.communicate(sequence)
    assert errors is None
    
    ans = ans.strip()
    if "\n" in ans:
        ans = "".
    print ans.strip().split(" ")
    ans = 

    energy = float(energy.strip("()"))
    
    return fold, energy
    
    


# seq2 = "uuggauguuggccuaguucuguguggaagacuagugauuuuguuguuuuuagauaacuaaaucgacaacaaaucacagucugccauauggcacaggccaugccucuacag"
#  
#  
# _viennafold(seq2)