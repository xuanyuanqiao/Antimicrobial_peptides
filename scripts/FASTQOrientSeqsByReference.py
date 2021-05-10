#!/usr/bin/env python3
"""
FASTQOrientSeqsByReference.py :  designed to take FLASH-joined reads as stdin, and produce a single column output with reads correctly oriented. 
e.g.,                                           :
flash -t 5 -M 150 -O -c  R1.fq R2.fq \
| FASTQOrientSeqsByReference.py -m 5 -r 'TAACATATG.*TAATAAGTCGACCTGCA' \
| xz -e -t 5 \
> seqs.flipped.lst.xz

@author: LBC
"""
import re, regex
import argparse, sys
import gzip, lzma
from collections import defaultdict
from pprint import pprint
from Bio import SeqIO
from Bio.Seq import Seq
from Levenshtein import distance

if __name__=="__main__":
    
    parse = argparse.ArgumentParser(description=__doc__)
    
    parse.add_argument("-i", "--input",
                       dest = "infile",
                       action = "store",
                       help = "Input [zipped] FASTQ file"
                       )
    parse.add_argument("-m", "--mismatches",
                        dest = "max_n_mismatches",action = "store",default = 0,
                        type=int,help = "Maximum number of mismatches (for regex search) (default = 0)"
                        )
    parse.add_argument("--write_fastq",action='store_true',
                        help = 'flag:  write FASTQ or just text (only the sequences)',
                        )
    parse.add_argument("--tab",action='store_true',
                        help = 'flag:  input format is two columns (e.g., from flash -To | cut -f 1,2 | FASTQOrientSeqsByReference.py)',
                        )
    #parse.add_argument("-o","--output_filename",
    #                    action = "store",type=str,
    #                    help = 'filename to save. default = stdout'
    #                    )
    group = parse.add_mutually_exclusive_group(required=True)
    group.add_argument('-r','--regex',action='store',type=str,help="expected sequence as a regular expression")
    group.add_argument('-s','--reference_sequence',action='store',type=str,help="expected seq as a string - matches entire len")
    opt = parse.parse_args()

# compile regex for search-with-mismatches
    reg = regex.compile('[ACTGNactgn].*' + '(' + opt.regex + ')' + '{s<=' + str(opt.max_n_mismatches) + '}')

# read a fasta or fastq file with optional .gz	or .xz
#   if no arg passed, defaults to sys.stdin
    file_format = 'fastq'
    if opt.tab: file_format = 'tab'
    if opt.infile:
        if re.search('FASTA' , opt.infile , re.IGNORECASE) :
            file_format = 'fasta'
        if re.search('gz' , opt.infile , re.IGNORECASE) :
            inhandle = gzip.open(opt.infile, "r")
        elif re.search('xz' , opt.infile , re.IGNORECASE) :
            inhandle = lzma.open(opt.infile, "rt")
        else:
            inhandle = open(opt.infile, "r")
    else:
        inhandle = sys.stdin

    def write_output(rec):
        if opt.write_fastq:
            SeqIO.write( rec , sys.stdout ,  file_format )
        else:
            print( str(rec.seq) )

    total_read_seqs = 0
    matched_Fwd_seqs_counter = 0
    matched_RevComp_seqs_counter = 0
    d = defaultdict(int)
    for record in SeqIO.parse(inhandle, file_format ):
        d['total_read_seqs'] += 1
        seq_fwd = str(record.seq)
        seq_rc = str(record.seq.reverse_complement())
        if opt.reference_sequence:
            dist_fwd = distance( opt.reference_sequence , seq_fwd)
            dist_rev = distance( opt.reference_sequence , seq_rc)
            if dist_fwd  > dist_rev:
                record = record.reverse_complement(id=record.id)
                flipped_seqs_counter = flipped_seqs_counter + 1
        elif opt.regex:
            m_fwd = reg.match(seq_fwd)
            if m_fwd:
                d['matched_Fwd_seqs_counter'] += 1
                write_output(record)
            else:
                m_rc = reg.match(seq_rc)
                if m_rc:
                    d['matched_RevComp_seqs_counter'] += 1
                    write_output(record.reverse_complement(id=record.id + '_rc'))

    d['NumberMatched'] = d['matched_Fwd_seqs_counter']+d['matched_RevComp_seqs_counter']
    d['PctMatched'] = d['NumberMatched'] / d['total_read_seqs'] * 100
    d['PctFlipped'] = d['matched_RevComp_seqs_counter'] / d['NumberMatched']  * 100
    pprint(d,compact=True,stream=sys.stderr)