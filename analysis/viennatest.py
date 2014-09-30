


# import subprocess
from subprocess import PIPE, Popen
# print subprocess.call("dir")

# subprocess.check_call("rnafold")


seq = "uuggauguuggccuaguucuguguggaagacuagugauuuuguuguuuuuagauaacuaaa" \
      "ucgacaacaaaucacagucugccauauggcacaggccaugccucuacag"

seq2 = "uuggauguuggccuaguucuguguggaagacuagugauuuuguuguuuuuagauaacuaaaucgacaacaaaucacagucugccauauggcacaggccaugccucuacag"



command = "rnafold"

mypipe  = Popen(command, stdin=PIPE, stdout=PIPE, bufsize=-1)

ans, errs = mypipe.communicate(seq2) 

fold, eng = ans.split(" ")
energy = float(eng.strip().strip("()"))

print fold
print energy

print errs

print "test"

def viennafold(sequence):
    """runs vienna folding, returns folding (string) and energy (double)"""
    mypipe  = Popen("rnafold", stdin=PIPE, stdout=PIPE, bufsize=-1)
    ans, errors = mypipe.communicate(seq2)
    assert errors is None
    fold, energy = ans.split(" ").strip()
    energy = float(energy.strip("()"))
    
    return fold, energy
    