#!/usr/bin/env python3

import pandas as pd
import os

annotation_file='seqinfo.csv'
script='script/demultiplex.py'

cluster_run = 'pkurun-fat4way 1 6'
starcode_distance = 2
starcode_cmd = f'starcode -d {starcode_distance} -t 2 '

anno = pd.read_csv(annotation_file, sep=',',index_col=None)

for i,row in anno.iterrows():

    flash_output_filename = 'out.joined.flipped.xz' 
    
    # does the flash output file for this library exist? 
    if not os.path.isfile(flash_output_filename):
        raise Exception('Cannot find file ' , flash_output_filename)

    # print the cmd to run demultiplex.py
    script_params = f'-f {row.SeqF1} -r {row.SeqR1}_rc -s {row.sampleName}'
    library_meta_data = [str(x) for x in [row.sampleName,row.SeqF1,row.SeqR1]] 
    stc_output_filename = '_'.join(library_meta_data)+f'.d{starcode_distance}' + '.stc'
    # stc_output_filename looks like : Ara_rep2_t2_GGGCTAGCGAATTCGAGCT_CAAGCTTGCATGCCTGCA.d2.stc
    full_cmd = f'{script} {script_params} -i {flash_output_filename} | {starcode_cmd} -o {stc_output_filename}'
    print(f"{cluster_run} \'{full_cmd}\'")
    #print(f"'{full_cmd}\'")

