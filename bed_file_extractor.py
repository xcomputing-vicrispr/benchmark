'''
A Benchmark of Computational CRISPR-Cas9 Guide Design Methods

Jacob Bradford, Dimitri Perrin. 2019.

The .bed files below were generated at:
	http://asia.ensembl.org/biomart
	
GRCm38-p6-mm10-chr19-exons.bed
	Dataset: Ensembl Genes 92 -> Mouse genes (GRCm38.p6)
	Filters: Chromosome/scaffold 19
	Attributes: Exon region start (bp), Exon region end (bp), Gene stable ID, Gene name
	
GRCm38-p6-mm10-chr19-genes.bed
	Dataset: Ensembl Genes 92 -> Mouse genes (GRCm38.p6)
	Filters: Chromosome/scaffold 19
	Attributes: Gene start (bp), Gene end (bp), Gene stable ID, Gene name
'''

files = [
	"C:\genomes\Biomart\GRCm38-p6-mm10-chr19-genes.bed",
	"C:\genomes\Biomart\GRCm38-p6-mm10-chr19-exons.bed"
]

startPos = 10000000
startPosWord = '10m' # uman readable version of startPos

lengths = [100000, 500000, 1000000, 5000000, 10000000, 20000000, None]
lengthsWords = ['100k', '500k', '1m', '5m', '10m', '20m', 'end']

chromosomeNames = [
    "ucsc-mm10-chr19-full-extract[10000000-10100000]",
    "ucsc-mm10-chr19-full-extract[10000000-10500000]",
    "ucsc-mm10-chr19-full-extract[10000000-11000000]",
    "ucsc-mm10-chr19-full-extract[10000000-15000000]",
    "ucsc-mm10-chr19-full-extract[10000000-20000000]",
    "ucsc-mm10-chr19-full-extract[10000000-30000000]",
    "ucsc-mm10-chr19-full",
    
]

i = 0
# we need a new file for 100k, 500k, 1m, 5m, etc
while i < len(lengths):
    length = lengths[i]
    newChromosomeName = chromosomeNames[i]
    
    # for both the exon and gene file
    for file in files:
    
        # read it
        with open(file, 'r') as fRead:
        
            # create a new file for the modified data
            with open('%s.%s-%s-adjusted.bed' % (file, startPosWord, lengthsWords[i]), 'w+') as fWrite:
                lines = fRead.read().split('\n')[:-1]
                
                # go line by line
                for line in lines:
                    lineSplit = line.split('\t')

                    # if length is null then we are adjusting from startPos to the end (eg: 10m to the end)
                    # if length is NOT null then we are adjusting from startPos to endPos (eg: 10m for length 10m (ending at 20m))
                    if (length == None and int(lineSplit[1]) >= startPos) or (length != None and int(lineSplit[1]) >= startPos and int(lineSplit[2]) < (startPos + length)):
                        
                        # make the adjustment and write to file
                        newStartPos = int(lineSplit[1]) - startPos
                        newEndPos = int(lineSplit[2]) - startPos
                        fWrite.write('%s\t%s\t%s\t%s\t%s\n' % (newChromosomeName, newStartPos, newEndPos, lineSplit[3], lineSplit[4]))           
    
    i = i + 1 # next length