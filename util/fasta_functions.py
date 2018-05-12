#!/usr/bin/env python
'''
Some shared Python functions for Opal helper scripts.
'''
import re
import string

def fasta_reader(f):
    '''Generator expression that returns a fasta sequence
        
        Ignores quality score string of FASTQ file
    '''
    seq = ''
    name = ''
    ignore_line = False
    first_line = True
    while True:
        line = f.readline()
        if line=='':
            if seq=='':
                break
            else:
                yield (name, seq)
                name = ''
                seq = ''
        elif line[0]=='>':
            if first_line:
                first_line = False
                pass
            else:
                yield (name, seq)
            name = line[1:].rstrip('\n')
            seq = ''
            ignore_line = False
        elif line[0]=='+':
            # Ignore quality score strings
            ignore_line = True
        else:
            if ignore_line:
                pass
            else:
                seq = seq + line.rstrip('\n')

trans = string.maketrans('ATGCatgc', 'TACGTACG')
def reverse_complement(dna):
    return dna[::-1].translate(trans)

def get_all_substrings(input_string, k):
    return [input_string[i:i+k] for i in xrange(len(input_string) - k + 1)]

pat = re.compile('^[ACGTacgt]*$')
def check_acgt(s):
    return pat.match(s)


