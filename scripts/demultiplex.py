#!/usr/bin/env python3

import argparse,sys,regex,lzma,gzip
from collections import defaultdict
from pprint import pprint
from Bio.Seq import Seq


if __name__=="__main__":
    
    parse = argparse.ArgumentParser(description=__doc__)
    
    parse.add_argument("-i", "--input_filename",type=str,action="store",help="read from input file instead of stdin",default=sys.stdin)
    parse.add_argument("-f", "--fwd_primer",type=str,action="store",help="forward primer sequence.",required=True)
    parse.add_argument("-r", "--rev_primer",type=str,action="store",help="end with _rc to use the reverse complement",required=True)
    parse.add_argument("-t", "--trim_to",help="regular expression to trim all seqs to (for removing primers and extra seq",required=False)
    parse.add_argument("--min_seq_length",help="discard seqs whose length (AFTER trimming) is less than this",type=int,default=80)
    parse.add_argument("--max_seq_length",help="discard seqs whose length (AFTER trimming) is more than this",type=int,default=300)
    parse.add_argument("-m", "--mismatches",dest = "max_n_mismatches",action = "store",default = 3,
                        type=int,help = "Maximum number of mismatches (for regex primer search) (default = 3)"
                        )
    parse.add_argument('-s',"--sampleID",action="append",
                        help = 'associate this primer pair with this sampleID. can pass multiple, e.g., -s glucose -s t0',
                        )
    opt = parse.parse_args()

if opt.rev_primer.endswith('_rc'):
    opt.rev_primer = str(Seq(opt.rev_primer[:-3]).reverse_complement())

# compile regex for search-with-mismatches
fwd_primer = regex.compile('(^' + opt.fwd_primer + ')' + '{s<=' + str(opt.max_n_mismatches) + '}')
rev_primer = regex.compile('(' + opt.rev_primer + '$)' + '{s<=' + str(opt.max_n_mismatches) + '}')
trim_to = regex.compile('(' + opt.trim_to + ')' + '{s<=' + str(opt.max_n_mismatches) + '}')

d = defaultdict(int)
d['fwd_primer'] = fwd_primer
d['rev_primer'] = rev_primer
d['sampleID'] = opt.sampleID
d['input_filename'] = opt.input_filename


# open compressed files; if not, stdin
if opt.input_filename:
    if opt.input_filename.rfind('.xz'):
        fh = lzma.open(opt.input_filename, mode='rt',errors='ignore')
    elif opt.input_filename.rfind('.gz'):
        fh = gzip.open(opt.input_filename, mode='rt',errors='ignore')

for line in map(str.rstrip, fh):
    d['# lines read'] += 1
    # first col is always sequence,  2nd column can be # reads; if not, skip it
    line_split = line.split('\t',maxsplit=2)
    Sequence = line_split[0] 
    
    # work for both 1 & 2 column inputs
    NReads = False
    if len(line_split)>1 and line_split[1].isnumeric():
        NReads = int(line_split[1])
        d['Total # read counts'] += NReads
    
    # if both our primers match, print the line
    if fwd_primer.match(Sequence) and rev_primer.findall(Sequence):
        d['# lines primer matched'] += 1
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
                if opt.sampleID: array_to_print.extend(opt.sampleID)
                if NReads:
                    array_to_print.insert(1,str(NReads))
                    d['Total # read counts Matched'] += NReads
                print('\t'.join(array_to_print))

                              
try:  # d['Total # read counts'] is 0 if not using starcode counts file as input
    d['% read counts matched'] = round( d['Total # read counts Matched'] / d['Total # read counts'] * 100 , 2 )
except:
    pass
d['% lines matched'] = round( d['# lines matched'] / d['# lines read'] * 100 , 2 )
sys.stdout.close() # maybe solves a problem on the cluster when piping into starcode
pprint(d,stream=sys.stderr)




