##process of BSubtilis sequencing data
# (1) FLASH join paired end reads to produce high confidence consensus sequences
# (2) flip sequences so that they have the correct orientation
# (3)use the sample meta data table (primers, replicates, conditions) to extract each experimental sample from the flash-joined & flipped list of sequences
#       this runs demultiplex.py
#       extracts seqs that begin and end with the sample primer pair
# (4)STARCODE to do PCR and sequencing error correction per-sample (after trimming)
# result:   one .stc file per sample;  all files should, in theory, have the same sequences
#

##Directories
DATADIR = /home/lucas_pkuhpc/lustre2/LABDATA/2019__AntiMicrobialPeptide/HHW_data/N2100438_HHW_80-588398309_ampliSEQ/combined
SCRIPTS = script

##Commands
XZ=xz --threads=2  --force --compress   --stdout
FLASHARGS = -t 4 -M 150 -O -m 50 -c -x 0.05
FLIP=$(SCRIPTS)/FASTQOrientSeqsByReference.py -m 5 -r 'TAAAGTGGTTATAACATAG.*CTCTTATACACATCTAATTAATTAATCGACCTGCAGG'

##run flash on all samples

flash:
	pkurun-fat4way 1 6 make _flash
_flash:
	flash --threads=2 $(DATADIR)/*_R1.fastq.gz $(DATADIR)/*_R2.fastq.gz $(FLASHARGS)\
		| $(FLIP) | $(XZ) > all.out.joined.flipped.xz

##starcode directly on flash joined reads
stc:
	pkurun-fat4way 4 24 make _stc T=6

_stc:
	xzcat all.out.joined.flipped.xz \
        | starcode -d 1 -t 8 -o out.joined.flipped.lst.stc.xz

#TRIM = $(SCRIPTS)/trim.py --min_seq_length=80 --max_seq_length=300 --trim_to='TAACATAGATGTGTAT.*TCTTATACACATCTAATTAATT'
STC = /home/lucas_pkuhpc/Develop/python/Sequence/starcode_split_into_samples_tidydata.py --trim_to='TAACATAGATGTGTAT.*TCTTATACACATCTAATTAATT' --min_seq_length=80 --max_seq_length=300

START:
	pkurun-cnlong 1 6 make _START
_START:
	$(STC) $(COMB8) -s Glu -s rep1 -s t0 -i all.out.joined.flipped.xz > t0_Glu_check_unique_reads.txt

## test the best distance for starcode
cut:
	pkurun-cnlong 1 6 make _cut

_cut:
	cut -f 1 t0_Glu_check_unique_reads.txt > t0_Glu_check_unique_reads_cut.txt

check_dist:
	pkurun-fat4way 1 8 make _check_dist

_check_dist:
	for i in {1..8} ; \
		do starcode -t 8 -d $$i -i t0_Glu_check_unique_reads_cut.txt -o check_distance2/out.count_$$i.stc; \
		done

##count the uniqe reads for each starcode distance
##the result shows that set the distance to 4 is OK.
wc:
	pkurun-cnnl 1 6 make _wc

_wc:
	for i in {1..8}; \
		do wc -l check_distance/out.count_$$i.stc >> check_distance/uniqreads.txt;\
		done

##combination of primers (parameters of SCRIPT for demultiplexing and pipeinto starcode)
COMB1 = -f GGGCTAGCGAATTCGAGCT -r CAAGCTTGCATGCCTGCA_rc
COMB2 = -f GGGCTAGCGAATTCGAGCT -r AACAGCCAAGCTTGCATGC_rc
COMB3 = -f GGGCTAGCGAATTCGAGCT -r GCCTGCAGGTCGATTAATTAATTAG_rc
COMB4 = -f GCGAATTCGAGCTCGGTACTA -r CAAGCTTGCATGCCTGCA_rc
COMB5 = -f GCGAATTCGAGCTCGGTACTA -r AACAGCCAAGCTTGCATGC_rc
COMB6 = -f GCGAATTCGAGCTCGGTACTA -r CATGCCTGCAGGTCGATTAATTAA_rc
COMB7 = -f GCGAATTCGAGCTCGGTACTA -r GCCTGCAGGTCGATTAATTAATTAG_rc
COMB8 = -f CGAGCTCGGTACTAAAGTGGTT -r CAAGCTTGCATGCCTGCA_rc
COMB9 = -f CGAGCTCGGTACTAAAGTGGTT -r AACAGCCAAGCTTGCATGC_rc
COMB10 = -f CGAGCTCGGTACTAAAGTGGTT -r CATGCCTGCAGGTCGATTAATTAA_rc

d=4
t=20
FILE = -i all.out.joined.flipped.xz


#parse_flash:
#	pkurun-fat4way 1 $t "$(STC) $(COMB1) $(FILE)| starcode -d $d -t $t -o STC_d3/Ara_rep2_t2.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB2) $(FILE)| starcode -d $d -t $t -o STC_d3/Ara_rep2_t3.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB3) $(FILE)| starcode -d $d -t $t -o STC_d3/Glu_t3.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB4) $(FILE)| starcode -d $d -t $t -o STC_d3/Ara_rep1_t1.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB5) $(FILE)| starcode -d $d -t $t -o STC_d3/Ara_rep1_t2.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB6) $(FILE)| starcode -d $d -t $t -o STC_d3/Ara_rep1_t3.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB7) $(FILE)| starcode -d $d -t $t -o STC_d3/Ara_rep2_t1.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB8) $(FILE)| starcode -d $d -t $t -o STC_d3/t0.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB9) $(FILE)| starcode -d $d -t $t -o STC_d3/Glu_t1.d$d.stc"
#	pkurun-fat4way 1 $t "$(STC) $(COMB10) $(FILE)| starcode -d $d -t $t -o STC_d3/Glu_t2.d$d.stc"
parse_flash:
	pkurun-fat4way 4 $t "$(STC) $(COMB1) $(FILE)| starcode -d $d -t $t -o STC_d4/Ara_rep2_t2.d$d.stc"
	pkurun-fat4way 4 $t "$(STC) $(COMB5) $(FILE)| starcode -d $d -t $t -o STC_d4/Ara_rep1_t2.d$d.stc"
	pkurun-fat4way 4 $t "$(STC) $(COMB6) $(FILE)| starcode -d $d -t $t -o STC_d4/Ara_rep1_t3.d$d.stc"

# xz -T 5 -f -z *.stc

##generate Putida_trimmed_parsed_samples.txt and analyze on jupyter

FILENAME = Putida_trimmed_parsed_samples.txt
parse_stc:
	pkurun-cnnl 1 $t make _parse_stc

_parse_stc:
	$(STC) $(COMB1) -s Ara -s rep2 -s t2 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB2) -s Ara -s rep2 -s t3 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB3) -s Glu -s rep1 -s t3 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB4) -s Ara -s rep1 -s t1 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB5) -s Ara -s rep1 -s t2 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB6) -s Ara -s rep1 -s t3 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB7) -s Ara -s rep2 -s t1 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB8) -s Glu -s rep1 -s t0 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB9) -s Glu -s rep1 -s t1 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	$(STC) $(COMB10) -s Glu -s rep1 -s t2 -i out.joined.flipped.lst.stc.xz >> $(FILENAME)
	wc -l $(FILENAME)

clean:
	rm -f *_*.out *_*.err job.srp*

xz:
	pkurun-cnnl 1 6 make _xz

_xz:
	xz -d -k all.out.joined.flipped.xz
