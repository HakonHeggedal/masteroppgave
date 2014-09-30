
from subprocess import PIPE, Popen



def _build_index():
    pass

_seq = "bowtie -c /home/hakon/Skrivebord/h_sapiens_37_asm TTCACAGTGGCTAAGTTCCGC"

def bowtiemap(sequence):
    
    mypipe  = Popen(sequence, stdin=PIPE, stdout=PIPE, bufsize=-1)
    answers, errors = mypipe.communicate(sequence)
    
    print answers
    print "errors", errors
    
    


bowtiemap(_seq)