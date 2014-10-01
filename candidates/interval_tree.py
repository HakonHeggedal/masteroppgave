'''
Created on 1. okt. 2014

@author: hakon
'''
import intervaltree

def findall(unfixed_lines):
    fixed_lines = [line.strip().split("\t") for line in unfixed_lines]
    
    print len(fixed_lines)
    print fixed_lines[0]
    