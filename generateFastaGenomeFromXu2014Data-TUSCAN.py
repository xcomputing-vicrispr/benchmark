'''
This script was written to generate a FASTA formatted genome file to test each 
CRISPR-Cas9 guide design tool against an experimentally validated dataset.

The dataset was obtained from:
    Xu, H., Xiao, T., Chen, C. H., Li, W., Meyer, C., Wu, Q., ... & Brown, M. (2015). Sequence determinants of improved CRISPR sgRNA design. Genome research, gr-191452.
    
    See Readme.txt for more information.

The format of the data provided by the citation above is the following:
    Note: an extra column was added to the end. It indicates whether the guide
    was determined to be efficient (1) or inefficient (0). See the paper for 
    more information.
    
    RPL10,chrX,153627684,-,GACTTTGGGTACGGCTTGTTCTTACAATACCGGTAACTGC,-1.001504666,-2.682526715,1

Expected input:
    A CSV file containing rows as described above.
        - The first line (header) is skipped in this script.
        - Expecting a blank line as the footer
        
Output:
    A FASTA formatted file, named: INPUT_FILE_NAME.fa
    
    <guideseq>NNNN...<guideseq>NNNN...<guideseq> etc
    
    There should be enough N's between each <guideseq> to prevent detection of
    overlapping sequences.

'''
import os

INPUT_FILE_NAME = "Xu-2015_Is-Efficient.csv"

NUM_N_BETWEEN_GUIDESEQ = 50

# Used for extracting the guide from the CSV format
GUIDE_INDEX = 4
GUIDE_START_POS = 10 - 4
GUIDE_LEN = 23 + 4 + 3 # TUSCAN considers nucleotides 4 upstream, 3 upstream from guide (?=([ACTG]{25}GG[ACTG]{3}))

with open(INPUT_FILE_NAME, 'r') as fRead, open('TUSCAN-%s.fa' % INPUT_FILE_NAME, 'w+') as fWrite:

    # read each line - file contains header and blank footer - then break apart by comma
    for line in [x.split(',') for x in fRead.read().split('\n')[1:-1]]:
        guideSeq = line[GUIDE_INDEX][GUIDE_START_POS:(GUIDE_START_POS + GUIDE_LEN)]
        guideGaps = 'N' * NUM_N_BETWEEN_GUIDESEQ
        fWrite.write('%s%s' % (guideSeq, guideGaps))