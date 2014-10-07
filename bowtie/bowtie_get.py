
from subprocess import PIPE, Popen, call, check_output






# _seq = "bowtie -c h_sapiens_37_asm TTCACAGTGGCTAAGTTCCGC"
# _seq = ["bowtie", "-c", "h_sapiens_37_asm", "TTCACAGTGGCTAAGTTCCGC"]
# _seq = ["bowtie","-f", "h_sapiens_37_asm", "mature.fa", "test.map"]
# _seq = ["bowtie","-f", "h_sapiens_37_asm", "SRR797062.fa", "test.map"]


def runcommand(params):
    print "run commands1"
#     mypipe  = Popen(params, stdin=PIPE, stdout=PIPE, bufsize=-1)
#     answers, errors = mypipe.communicate(params)
#     print check_output(["echo", "hello world"])
    outs = check_output(params).strip()
    print "run commands2"        
#     print answers
#     print "errors", errors
    parts = outs.split("\t")
    print len(outs), outs
    print len(parts), parts
    


# runcommand(_seq)