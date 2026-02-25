'''
A Benchmark of Computational CRISPR-Cas9 Guide Design Methods

Jacob Bradford, Dimitri Perrin. 2019.

This scripts calculates the portion of guides reported by one tool, which would
have been rejected by mm10db. This applies to the exon-targeting guides reported
by tools when tested against the 500k dataset.

Run:
    1:  (normalised data for 500k dataset required; see: normalise.py)
    2:  (extract exon targeting guides using: normalised-extract-exon-guides.py)
    3:  python mm10db-rejects-accepted-by-other-tools.py

Output:
    - CSV formatted data
'''

import os

mm10dbAccepts = r"mm10dbAccepted.normalised"

mm10dbRejects = r"mm10dbRejected.tsv"

acceptedTargetsDir = r"normalised/exon-only/500k"

toolsToConsider = [
        'casdesigner',
        'casfinder',
        'cctop',
        'chopchop',
        'crispor',
        'crisprdo',
        'crisprera',
        'ctfinder',
        'flashfry',
        'gtscan',
        'guidescan',
        'mm10db',
        'phytocrispex',
        'sgrnacas9',
        'sgrnascorer',
        'ssc',
        'tuscan',
        'wucrispr',
        'vicrispr'
    ]
#toolsToConsider = ['guidescan']

def doesFileNameMatchToolsToConsider(fileName):
    matches = False
    for toolToConsider in toolsToConsider:
        if toolToConsider.lower() in fileName.replace('-', '').lower():
            matches = True
    return matches

mm10RejectGuides = {}
with open(mm10dbRejects, 'r') as fRead:
    for z in fRead.read().split('\n')[:-1]:
        x = z.split('\t')
        mm10RejectGuides[x[0]] = x[1]
        mm10RejectGuides[x[0][:-3]+'...'] = x[1] # special case where PAM was not concatenated so we did it manually in normalisation
        mm10RejectGuides[x[0][:-3]+'NGG'] = x[1] # special case where NGG was concatenated and not actual sequence from FASTA file

mm10AcceptGuides = {}
with open(mm10dbAccepts, 'r') as fRead:
    for z in fRead.read().split('\n')[:-1]:
        x = z.split(',')
        # mm10crisprdatabase,AAAGTCTTATCTGCAGAGAG...,916614,916636,+
        mm10AcceptGuides[x[1]] = None
        mm10AcceptGuides[x[1][:-3]+'...'] = None # special case where PAM was not concatenated so we did it manually in normalisation
        mm10AcceptGuides[x[1][:-3]+'NGG'] = None # special case where NGG was concatenated and not actual sequence from FASTA file
        
headerPrinted = False    

# walk over the directory
for (dirpath, dirnames, filenames) in os.walk(acceptedTargetsDir):
                    
    # for file within the directory
    for filename in filenames:
        if not doesFileNameMatchToolsToConsider(filename):
            continue

        fullpath = os.path.join(acceptedTargetsDir, filename)
        with open(fullpath, 'r') as fNormalisedRead:
            # casdesigner,TGGGGCGGAGCCGCAAACCCAGG,18214,18236,+
            toolAcceptedGuides = [z.split(',')[1][:-3]+'NGG' for z in fNormalisedRead.read().split('\n')[:-1]]
         

        stats = {
            # the number of accepted guides that were rejected by mm10db because...
            'AT%' : 0,
            'Secondary structure or energy' : 0,
            'Multiple matches in exons' : 0,
            'TTTT' : 0,
            'Multiple matches in genome' : 0,
            'Too close to reverse primer' : 0,
            'Off-target score' : 0,
            
            # tool A produced the guide and mm10db accepted it too
            'both accepted' : 0, # this is inaccurate
            
            # tool A produced the guide but mm10db did not accept or reject it
            'not considered by mm10db' : 0,
        }
        
        if not headerPrinted:
            headerPrinted = True
            print ','.join(['toolName', 'acceptedGuides'] + stats.keys())
          
        # find which rejected guides were accepted
        temp = 0
        for toolAcceptedGuide in toolAcceptedGuides:
        
            if toolAcceptedGuide not in mm10RejectGuides and toolAcceptedGuide not in mm10AcceptGuides:
                temp += 1
            if toolAcceptedGuide in mm10RejectGuides:
                # the accepted guide was rejected by mm10db
                stats[mm10RejectGuides[toolAcceptedGuide]] = stats[mm10RejectGuides[toolAcceptedGuide]] + 1
                continue
                
            # mm10db did not even consider the guide
            # ie: guide was not accepted and not rejected by mm10db
            if toolAcceptedGuide not in mm10AcceptGuides:
                stats['not considered by mm10db'] = stats['not considered by mm10db'] + 1
                continue
            stats['both accepted'] += 1
     
        print '%s' % ','.join([filename, str(len(toolAcceptedGuides))] + map(str, stats.values()))