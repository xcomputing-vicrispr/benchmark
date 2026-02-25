'''
A Benchmark of Computational CRISPR-Cas9 Guide Design Methods

Jacob Bradford, Dimitri Perrin. 2019.

This script generates the heatmap plots in the paper. This is the most complex
plotting script to run successfully.

Run:
    python plot-guideConsensusHeatmaps.py
    
Input:
    - normalised data (with scores) in directory specified in TOOLS_NORMALISED_DATA_DIR 
    - (process the raw data using normalise-Xu2015.py)
    
Output:
    - image file
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

PLOT_OUTPUT_DPI = 300

for FOR_EXON_GUIDES_ONLY in (True, False):
    OUTPUT_DIR = '.'

    if FOR_EXON_GUIDES_ONLY:
        tools = [    
            'Cas-Designer',
            'CasFinder', 
            'CCTop',
            'CHOPCHOP',  
            'CRISPOR',
            'CRISPR-DO',
            'CRISPR-ERA', 
            'CT-Finder', 
            'FlashFry',
            'GT-Scan', 
            'GuideScan',
            'mm10db', 
            'PhytoCRISP-Ex',
            'sgRNAcas9',
            'sgRNAScorer2',
            'SSC',  
            'TUSCAN',
            'WU-CRISPR'
        ]

        # "C:\git\IEEE-BIBM-Paper\scripts\compare-one-with-all-others-output\500k_exons_results.csv"
        dataTemp = [map(float, x.split('\t')) for x in '''1	0.994845361	1	0.27867268	1	0.569104381	1	1	1	1	1	0.185244845	0.160760309	0.912210052	0.951514175	1	1	0.17123067
0.880022799	1	1	0.2799943	1	0.572385295	1	1	1	1	1	0.178113423	0.157737247	0.91892277	0.952265603	1	1	0.170703904
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.872855701	0.991422805	1	1	1	0.823915237	1	1	1	1	1	0.147830474	0.157921292	0.894046418	0.974772957	1	1	0.352169526
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.879512074	1	1	0.40652228	1	1	1	1	1	1	1	0.17326363	0.156086632	0.922827981	0.947722181	1	1	0.253422952
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.92	1	1	0.2344	1	0.5568	1	1	1	1	1	1	0.2056	0.996	1	1	1	0.2024
0.896675651	0.994609164	1	0.281221923	1	0.563342318	1	1	1	1	1	0.230907457	1	0.909254268	1	1	1	0.163522013
0.878120639	1	1	0.274771282	1	0.574817801	1	1	1	1	1	0.193053187	0.156923554	1	0.954101411	1	1	0.182508916
0.879017857	0.994494048	1	0.2875	1	0.566517857	1	1	1	1	1	0.186011905	0.165625	0.915625	1	1	1	0.168303571
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.879693921	0.994473572	1	0.280855888	1	0.569222049	1	1	1	1	1	0.177129092	0.157715743	0.91384441	0.952245997	1	1	0.169760521
0.887312187	1	1	0.58263773	1	0.849749583	1	1	1	1	1	0.211185309	0.151919866	0.982470785	0.944073456	1	1	1'''.split('\n')]

        data = np.array(dataTemp)

    else:
        tools = [
            'CasFinder',            
            'CCTop',
            'CHOPCHOP', 
            'CRISPOR', 
            'CRISPR-DO',
            'CRISPR-ERA', 
            'CT-Finder', 
            'FlashFry',
            'GT-Scan', 
            'GuideScan', 
            'PhytoCRISP-Ex',
            'sgRNAcas9',
            'sgRNAScorer2', 
            'SSC',
            'TUSCAN',
            'WU-CRISPR'
        ]

        dataTemp = [map(float, x.split('\t')) for x in 
        '''1	1	0.261674675	1	0.547737947	1	1	1	1	1	0.155104465	0.942018032	0.953665713	1	0.999973007	0.174404794
0.986814416	1	0.261021283	1	0.540515703	1	1	1	1	1	0.154551024	0.929596974	0.953770544	1	0.999973362	0.172105165
0.989284621	1	1	1	0.802020614	1	1	1	1	1	0.143075824	0.926216961	0.970456169	1	1	0.346463925
0.986814416	1	0.261021283	1	0.540515703	1	1	1	1	1	0.154551024	0.929596974	0.953770544	1	0.999973362	0.172105165
1	1	0.387305029	1	1	1	1	1	1	1	0.150875983	0.971539807	0.95182712	1	1	0.270433433
0.986788131	0.999973363	0.261014331	0.999973363	0.540501305	1	0.999973363	0.999973363	0.999973363	0.999973363	0.154546907	0.929572212	0.953745139	0.999973363	0.999946726	0.172100581
0.986814416	1	0.261021283	1	0.540515703	1	1	1	1	1	0.154551024	0.929596974	0.953770544	1	0.999973362	0.172105165
0.986814416	1	0.261021283	1	0.540515703	1	1	1	1	1	0.154551024	0.929596974	0.953770544	1	0.999973362	0.172105165
0.986814416	1	0.261021283	1	0.540515703	1	1	1	1	1	0.154551024	0.929596974	0.953770544	1	0.999973362	0.172105165
0.986814416	1	0.261021283	1	0.540515703	1	1	1	1	1	0.154551024	0.929596974	0.953770544	1	0.999973362	0.172105165
0.990348156	1	0.241640814	1	0.527662875	1	1	1	1	1	1	0.935884178	0.987676663	1	1	0.17847294
1	1	0.26007221	1	0.564903433	1	1	1	1	1	0.155596309	1	0.953550347	1	0.999971345	0.184595106
0.986705953	1	0.265587689	1	0.539414336	1	1	1	1	1	0.160045244	0.929382357	1	1	0.999972071	0.173073969
0.986814416	1	0.261021283	1	0.540515703	1	1	1	1	1	0.154551024	0.929596974	0.953770544	1	0.999973362	0.172105165
0.986814065	1	0.261028237	1	0.540530101	1	1	1	1	1	0.154555141	0.929595099	0.953769313	1	1	0.17210975
1	1	0.525460455	1	0.84932673	1	1	1	1	1	0.160269308	0.997059279	0.959139452	1	1	1.099365423'''
    .split('\n')]

        
        data = np.array(dataTemp)

    fig, ax = plt.subplots(figsize=(11,10))
    im = ax.imshow(data, cmap=cm.summer_r) #bone

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Consensus Rate', rotation=-90, va="bottom")
    ticklbls = ['%i%%' % (z * 100) for z in map(float, cbar.get_ticks())]
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(ticklbls)
    #cbar.set_ticklabels(['92%', '96%', '100%'])

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(tools)))
    ax.set_yticks(np.arange(len(tools)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(tools)
    ax.set_yticklabels(tools)
    ax.set_ylabel('Tool B')
    ax.set_xlabel('Tool A')
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(tools)):
        for j in range(len(tools)):
            text = ax.text(j, i, '%s%%' % (round(data[i, j] * 100, 1)),
                           ha="center", va="center", color="black", size=7)

    

    if FOR_EXON_GUIDES_ONLY:
        #ax.set_title("Consensus of Exon Targeting Guides Across Tools")
        outFigFileName = '%s/guideConsensusExonsIEEE.eps' % (OUTPUT_DIR)
    else:
        #ax.set_title("Consensus of Guides Across Non-Annotating Tools")
        outFigFileName = '%s/guideConsensusIEEE.eps' % (OUTPUT_DIR)
    
    fig.tight_layout()
    outFig = ax.get_figure()    
    outFig.savefig(outFigFileName, format='eps', dpi=PLOT_OUTPUT_DPI)
    print 'Wrote to: %s' % outFigFileName
    plt.clf()