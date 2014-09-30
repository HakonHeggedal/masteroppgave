
from subprocess import PIPE, Popen



def _build_index():
    pass



def bowtiemap(sequence):
    """runs vienna folding, returns folding (string) and energy (double)"""
    mypipe  = Popen("bowtie", stdin=PIPE, stdout=PIPE, bufsize=-1)
    answers, errors = mypipe.communicate(sequence)
#     fold, energy = ans.split(" ").strip()
#     energy = float(energy.strip("()"))
    
#     return fold, energy