#!/usr/bin/env python3

import argparse,sys,regex,lzma,gzip
from collections import defaultdict
from pprint import pprint
from Bio.Seq import Seq

if __name__=="__main__":

    parse = argparse.ArgumentParser(description=__doc__)
    parse.add_argument("-t", "--trim_to",help="regular expression to trim all seqs to (for removing primers and extra seq",required=True)
    parse.add_argument("--min_seq_length",help="discard seqs whose length (AFTER trimming) is less than this",type=int,default=0)
    parse.add_argument("--max_seq_length",help="discard seqs whose length (AFTER trimming) is more than this",type=int,default=1000000)
    parse.add_argument("-i", "--input_filename",type=str,action="store",help="read from input file instead of stdin",default=sys.stdin)
    parse.add_argument("-m", "--mismatches",dest = "max_n_mismatches",action = "store",default = 2,
                        type=int,help = "Maximum number of mismatches (for regex primer search) (default = 2)")
    opt = parse.parse_args()

d = defaultdict(int)
trim_to = regex.compile('(' + opt.trim_to + ')' + '{s<=' + str(opt.max_n_mismatches) + '}')

if opt.input_filename:
    if opt.input_filename.rfind('.xz'):
        fh = lzma.open(opt.input_filename, mode='rt',errors='ignore')
    elif opt.input_filename.rfind('.gz'):
        fh = gzip.open(opt.input_filename, mode='rt',errors='ignore')

for line in map(str.rstrip, fh):
    d['# lines read'] += 1
    line_split = line.split('\t',maxsplit=2)
    Sequence = line_split[0]
    seq_trim_to = trim_to.findall(Sequence)
    if seq_trim_to:
        if len(seq_trim_to[0]) < opt.min_seq_length:
            d['# lines too short'] += 1
        elif len(seq_trim_to[0]) > opt.max_seq_length:
            d['# lines too long'] += 1
        else:
            d['# lines matched'] += 1
            array_to_print = list()
            array_to_print.append(seq_trim_to[0])
            print('\t'.join(array_to_print))

try:  # d['Total # read counts'] is 0 if not using starcode counts file as input
    d['% read counts matched'] = round( d['Total # read counts Matched'] / d['Total # read counts'] * 100 , 2 )
except:
    pass
d['% lines matched'] = round( d['# lines matched'] / d['# lines read'] * 100 , 2 )
sys.stdout.close() # maybe solves a problem on the cluster when piping into starcode
pprint(d,stream=sys.stderr)
