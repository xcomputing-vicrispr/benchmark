'''
A Benchmark of Computational CRISPR-Cas9 Guide Design Methods

Jacob Bradford, Dimitri Perrin. 2019.

The intentions of this script is to normalise the data produced in the tests using
the experimentally validated dataset provided in:

Xu, H., Xiao, T., Chen, C. H., Li, W., Meyer, C., Wu, Q., ... & Brown, M. (2015). Sequence determinants of improved CRISPR sgRNA design. Genome research, gr-191452.

Notes:
    - The difference between normalise-Xu2015.py and normalise.py is:
        o This script appends the score reported by the guide design tool to the
            normalised format.
Run:
    python normalise-Xu2015.py
    
Input:
    The arg parser requires the following:
        parser.add_argument('-i', help='Input file', default=None, required=True)
        parser.add_argument('-t', help='Tool name', default=None)
        parser.add_argument('-o', help='Output file', default=None)
    
Output:
    - Normalised data placed into directory OUTPUT_DIR
'''

import argparse, os, re, math

OUTPUT_DIR = r'normalisedOutputExp'
OUT_FILE_EXT = 'normalised' # do not modify

def main(inputFileUrl, toolName, outputFile):
    if toolName == None:
        toolName = raw_input('Enter the tool name: ')
        
    toolName = normaliseToolName(toolName)

    if outputFile == None:
        outputFile = os.path.join(OUTPUT_DIR, '%s.%s' %(inputFileUrl.split('\\')[-1], OUT_FILE_EXT))
    
    if os.path.exists(inputFileUrl) == False:
        print 'Input file does not exist: %s' % inputFileUrl
        exit()
    
    # GT-Scan is a special case (uses SQLite DB)
    if toolName == normaliseToolName('GT-Scan'):
        with open(outputFile, 'w') as fWrite:
            import sqlite3
            conn = sqlite3.connect(inputFileUrl)

            c = conn.cursor()

            c.execute("SELECT * FROM targets ORDER BY offset ASC")
            for target in c.fetchall():
                seq = target[2] # <span>5'-</span><span class='b'>TGGACCGAGAAT</span><span class='y'>CTCTG</span><span class='g'>T</span><span class='y'>GG</span>-3'
                seq = seq[32:44] + seq[67:75] + seq[98:99] + seq[122:124]
                startPos = int(target[0])
                endPos = startPos + 22
                strand = target[1]
                if 'N' not in seq:
                    fWrite.write(generateNormalisedRecord(toolName, seq, startPos, endPos, strand, None))
                
    else: # every other tool
        with open(inputFileUrl, 'r') as inputFile:
            with open(outputFile, 'w') as fWrite:
                
                if toolName == normaliseToolName('Cas-Designer'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first two lines as they are headers
                    fReadLines = fReadLines[2:]
                    
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 1:
                            start = int(lineSplit[1]) + 1
                            end = start + 22
                            
                            if 'N' not in lineSplit[0]:
                                fWrite.write(generateNormalisedRecord(toolName, lineSplit[0], start, end, lineSplit[3], lineSplit[5]))

                if toolName == normaliseToolName('CHOPCHOP'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    for line in fReadLines[1:]:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 1:
                            #   1	CACACACAGAGTCTACATTGAGG	ucsc-mm10-chr19-full-extract[10000000-10500000]:178	1	+	48	0	0	0	0	0	1.00
                            positionSplit = lineSplit[2].split(':')
                            start = int(positionSplit[1])
                            end = start + len(lineSplit[1]) - 1
                            strand = lineSplit[4]
                            
                            if 'N' not in lineSplit[1]:
                                fWrite.write(generateNormalisedRecord(toolName, lineSplit[1], start, end, strand, lineSplit[11]))
                
                if toolName == normaliseToolName('CRISPOR'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 2:
                            
                            m = re.search('(\d+)(rev|fw)', lineSplit[0])     
                            strand = {'fw': '+', 'rev': '-'}[m.group(2)]
                            if m.group(2) == 'fw':
                                start = int(m.group(1)) - 20
                                end = int(m.group(1)) + 2
                            else:
                                start = m.group(1)
                                end = int(m.group(1)) + 22
                        
                            if 'N' not in lineSplit[1]:
                                fWrite.write(generateNormalisedRecord(toolName, lineSplit[1], start, end, strand, lineSplit[2]))
                
                if toolName == normaliseToolName('GuideScan'): # good
                    for line in inputFile.readlines():
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 1:
                        
                            temp = lineSplit[1].split(':')
                            
                            start = int(temp[1]) + 1
                            end = start + 22
                            strand = temp[2].strip()
                            seq = lineSplit[0]
                            if seq[-3:] == 'NAG':
                                continue
              
                            if strand == '-':
                                # accurate:     chopchop	TAGGAGCCCACAGAGATTCTCGG	16	38	-
                                # inaccurate:   guidescan	TAGGAGCCCACAGAGATTCTNGG	38	60	-
                                # to fix, do this:
                                end = start
                                start = end - 22
                            
                            if 'N' not in seq[:-3]:
                                fWrite.write('%s' % generateNormalisedRecord(toolName, seq, start, end, strand, None))

                if toolName == normaliseToolName('Cas-Finder'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 5:
                        
                            temp = lineSplit[10].split(':')[1].split('-')
                            start = temp[0]
                            end = temp[1]
                            strand = lineSplit[11]
                            
                            if 'N' not in lineSplit[13]:
                                fWrite.write(generateNormalisedRecord(toolName, lineSplit[13], start, end, strand, None))
          
                if toolName == normaliseToolName('CRISPR-ERA'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 3:
                            # https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A10000019-10000038&hgsid=685137703_scDj1TrdeXNtK2LK0vucTilfREsW
                            # TAGGAGCCCACAGAGATTCT		18	+	10   (UCSC: - chr19:10,000,019-10,000,038)
                            # AGCTGGACCGAGAATCTCTG		8	-	10   (UCSC: + chr19:10,000,009-10,000,028)
                            #start = int(lineSplit[2]) + 1
                            #end = start + 19
                            
                            strand = {'-' : '+', '+' : '-'}[lineSplit[3]] # reverse the strand
                            seq = '%s...' % lineSplit[0]
                            if strand == '+':
                                start = int(lineSplit[2]) + 1
                                end = (int(lineSplit[2]) + 1) + 19 + 3
                            else:
                                start = (int(lineSplit[2]) + 1) - 3
                                end = (int(lineSplit[2]) + 1) + 19
                            
                            if 'N' not in seq:
                                fWrite.write(generateNormalisedRecord(toolName, seq, start, end, strand, None))
                
                if toolName == normaliseToolName('sgRNAScorer2.0'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 1:
                        
                            if False:
                                temp = lineSplit[0].split('_')
                                strand = {'Plus': '+', 'Minus': '-'}[temp[1]]
                                start = temp[2]
                                end = int(start) + len(lineSplit[1])
                            else:
                                temp = lineSplit[0].split('_')
                                strand = {'Plus': '+', 'Minus': '-'}[temp[1]] # they're also reversed
                                
                                start = int(temp[2]) + 1
                                end = start + 22
                             
                            if 'N' not in lineSplit[1]:
                                fWrite.write(generateNormalisedRecord(toolName, lineSplit[1], start, end, strand, lineSplit[2]))
                
                if toolName == normaliseToolName('CT-Finder'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = [x.split('\t') for x in fReadLines[1:-1]]
                    
                    # we need to determine what the input size was (tidy this?)
                    # 1	TAGGAGCCCACAGAGATTCTCGG	-	499963	499985	0.50	0.33	1
                    maxPositionValue = max([max(int(z[3]),int(z[4])) for z in fReadLines])
                    sizes = [500000, 1000000, 5000000, 61431566]
                    inputSize = sizes[0]
                    for z in xrange(len(sizes)):
                        if maxPositionValue > sizes[z]:
                            inputSize = sizes[min((len(sizes)-1),z+1)]

                    for lineSplit in fReadLines:
                        if len(lineSplit) > 1:
                            strand = lineSplit[2]
                            seq = lineSplit[1]
                            ogStart = int(lineSplit[3])
                            ogEnd = int(lineSplit[4])
                            
                            if strand == '-':
                                    # accurate:     chopchop	TAGGAGCCCACAGAGATTCTCGG	16	38	-
                                    # inaccurate:   ctfinder	TAGGAGCCCACAGAGATTCTCGG	499963	499985	-
                                    # to fix, do this:
                                    start = inputSize - ogEnd + 1
                                    end = inputSize - ogStart + 1
                            else:
                                start = ogStart
                                end = ogEnd
                        
                            if '[ACGT]' not in seq and 'N' not in seq:
                                fWrite.write(generateNormalisedRecord(toolName, seq, start, end, strand, None))

                if toolName == normaliseToolName('mm10-CRISPR-Database'): # good
                    # AATAGATGGATGGTACCCAC	GAUAGAUGGAUGGUACCCACGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU	(((.(((((((.....((..(((((((.((((....))))...)))))))..))...)))))))))).((((....))))(((((((...)))))))...	-28.90	100.000000	ucsc-mm10-chr19-full-extract[10000000-11000000]	675062	675081	+	Pga5	CDS	CACTATAGGATAGATGGATGGTACCCACgttttagagctaGAAAtagc	gggccTAATACGACTCACTATAGGATAGATGGATGGTACCCACg
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # we might run into duplicates because mm10db exports the same guide for every gene it hits
                    guidesSeen = {'+' : {}, '-' : {}} 
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 1:
                            strand = lineSplit[8]
                            seq = '%s...' % lineSplit[0]
                            if strand == '+':
                                start = int(lineSplit[6]) - 1
                                end = int(lineSplit[7]) + 2
                                hashKey = end
                            else:
                                start = int(lineSplit[6]) - 3
                                end = int(lineSplit[7])
                                hashKey = start
                                
                            if hashKey not in guidesSeen[strand]:
                                guidesSeen[strand][hashKey] = lineSplit
                                if 'N' not in seq:
                                    fWrite.write(generateNormalisedRecord(toolName, seq, start, end, strand, None))
                            else:
                                print 'mm10 repeat'

                if toolName == normaliseToolName('WU-CRISPR'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 1:
                            start = int(lineSplit[6])
                            end = int(start) + 22
                            strand = {'sense' : '+', 'antisense' : '-'}[lineSplit[4]]
                            if strand == '-':
                                start += 1
                                end += 1
                            score = math.floor((1 - float(lineSplit[1])) * 100 + 0.5)
                            if score >= 50:
                                fWrite.write(generateNormalisedRecord(toolName, lineSplit[8][:-1].upper(), start, end, strand, score))
                
                if toolName == normaliseToolName('PhytoCRISP-Ex'): # good
                    # ucsc-mm10-chr19-full-extract[10000000-10500000]_-_100575_100597_CCTGAATTTTGCAAAATCCAGCC,PASS,PASS,,HpyCH4V->9;Tsp509I->4
                    # ucsc-mm10-chr19-full-extract[10000000-10500000]_+_100638_100660_GGCTGGCACTTGGAGTGAGTCGG,PASS,PASS,MlyI->16;PleI->16,
                    
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        lineSplit = line.split(',')[0].split('_')
                        if len(lineSplit) >= 5:

                            start = int(lineSplit[2])
                            end = int(lineSplit[3])
                            strand = lineSplit[1]
                            
                            fWrite.write(generateNormalisedRecord(toolName, lineSplit[4], start, end, strand, None))       

                if toolName == normaliseToolName('FlashFry'): # good
                    # ucsc-mm10-chr19-full-extract[10000000-10500000]	8	31	AGCTGGACCGAGAATCTCTGTGG	ATGTACAGCTGGACCGAGAATCTCTGTGGGCTCCT	OK	FWD	1	AGCTGGACCGAGAATCTCTGTGG_1_0
                    
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 5:

                            start = int(lineSplit[1]) + 1
                            end = lineSplit[2]
                            
                            if lineSplit[6] == 'FWD':
                                strand = '+'
                            if lineSplit[6] == 'RVS':
                                strand = '-'
                            
                            fWrite.write(generateNormalisedRecord(toolName, lineSplit[3], start, end, strand, lineSplit[7]))
                            
                            
                if toolName == normaliseToolName('sgRNAcas9'): # good
                    # ucsc-mm10-chr19-full-extract[10000000-10500000]_S_4	178	200	CACACACAGAGTCTACATTGAGG	23	45.0 %
                    # ucsc-mm10-chr19-full-extract[10000000-10500000]_A_8457	400787	400809	AGGAAAGAAAGAGAAACCACTGG	23	40.0 %
                    
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 5:

                            start = lineSplit[1]
                            end = lineSplit[2]
                            
                            if '_S_' in lineSplit[0]:
                                strand = '+'
                            if '_A_' in lineSplit[0]:
                                strand = '-'
                            
                            fWrite.write(generateNormalisedRecord(toolName, lineSplit[3], start, end, strand, None))
                            
                if toolName == normaliseToolName('SSC'): # good
                    # AGCTGGACCGAGAATCTCTGTGGGCTCCTA	8	37	+	ucsc-mm10-chr19-full-extract[10000000-10500000]
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) >= 5:
                            
                            if 'N' in lineSplit[0][:23]:
                                continue
                            start = int(lineSplit[1]) + 1
                            end = start + 22
                            strand = lineSplit[3]
                            
                            fWrite.write(generateNormalisedRecord(toolName, lineSplit[0][:23], start, end, strand, lineSplit[5]))
                
                
                if toolName == normaliseToolName('CCTop'): # good
                    '''
                    We need to read two files for CCTop because the aggregated 
                    file is in a horrible format.
                    
                    1) read the *.bed file 
                        ucsc-mm10-chr19-full-extract[10000000-11000000]	8	31	T1	1000	+
                    2) read the *.fasta file
                        >T1
                        AGCTGGACCGAGAATCTCTGTGG
                    '''
                    
                    # we have already opened the .bed file, so now open the .fasta file too
                    seqs = []
                    with open('%s%s' % (inputFileUrl[:-3], 'fasta'), 'r') as fFastaRead:
                        lineNum = 0
                        for line in fFastaRead.readlines():
                            if lineNum % 2 == 1:
                                seqs.append(line.strip())
                            lineNum = lineNum + 1
                            
                    # now for each line in the .bed file
                    lineNum = 0
                    for line in inputFile.readlines():
                        lineSplit = line.strip().split('\t')
                        if len(lineSplit) > 4:
                            fWrite.write(generateNormalisedRecord(toolName, seqs[lineNum], int(lineSplit[1]) + 1, lineSplit[2], lineSplit[5], None))
                        
                        lineNum = lineNum + 1      

                if toolName == normaliseToolName('CRISPR-DO'): # good
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        lineSplit = line.split('\t')
                        if len(lineSplit) > 1:
                            seq = lineSplit[3][:-7] # remove the 3' flanking sequence
                            start = int(lineSplit[1]) + 1
                            end = start + len(seq) - 1
                            if lineSplit[4] == '-':
                                start += 7
                                end += 7
                            fWrite.write(generateNormalisedRecord(toolName, seq, start, end, lineSplit[4], lineSplit[5]))
                            
                if toolName == normaliseToolName('GT-Scan1.3'):
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    for line in fReadLines:
                        lineSplit = line.split(', ')
                        if len(lineSplit) > 1:
                            seq = lineSplit[2]
                            start = int(lineSplit[0])
                            end = start + 22
                            fWrite.write(generateNormalisedRecord(toolName, seq, start, end, lineSplit[1], None))    
                    
                 
                if toolName == normaliseToolName('TUSCAN'):
                    # ucsc-mm10-chr19-full-extract[10000000-10500000] 9          32         +        AGCTGGACCGAGAATCTCTGTGG         3.94072587300098
                    fRead = inputFile.read()
                    fReadLines = fRead.split('\n')
                    # skip first line (header)
                    fReadLines = fReadLines[1:]
                    for line in fReadLines:
                        r = re.search('[^\n\r]([^\s]*)\s*([^\s]*)\s*([^\s]*)\s*([^\s]*)\s*([^\s]*)\s*([^\s]*)', line)
                        if r is not None and len(r.groups()) > 1:
                            seq = r.group(5)
                            start = int(r.group(2))
                            end = int(r.group(3)) - 1
                            fWrite.write(generateNormalisedRecord(toolName, seq, start, end, r.group(4), r.group(6)))    
                            
                            
def generateNormalisedRecord(normalisedToolName, sgRNA, posStart, posEnd, strand, score):
    return '%s\n' % (','.join(map(str, [normalisedToolName, sgRNA, posStart, posEnd, strand, score])))

def normaliseToolName(toolName):
    return toolName.lower().replace('-', '').replace('_', '').replace('.', '')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Normalise CRISPR Guide Design Tool output')
    parser.add_argument('-i', help='Input file', default=None, required=True)
    parser.add_argument('-t', help='Tool name', default=None)
    parser.add_argument('-o', help='Output file', default=None)
    args = parser.parse_args()
    
    main(args.i, args.t, args.o)
   