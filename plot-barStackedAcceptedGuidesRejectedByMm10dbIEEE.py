'''
A Benchmark of Computational CRISPR-Cas9 Guide Design Methods

Jacob Bradford, Dimitri Perrin. 2019.

This draws the stacked bar graph used in the paper.
The plot describes what percentage of exon-targeting guides reported by a 
particular tool would have been rejected by mm10db. The data for this is 
calculated by mm10db-rejects-accepted-by-other-tools.py and pasted here. 

Run:
    1:  (normalise the raw data using: normalise.py)
    2:  (extract exon targeting guides using: normalised-extract-exon-guides.py)
    3:  (calculate data for this script using: mm10db-rejects-accepted-by-other-tools.py)
    4:  python plot-barStackedAcceptedGuidesRejectedByMm10dbIEEE.py
    
Input:
    - TSV formatted data pasted into the multiline string named "raw_data".
    
Output:
    - image file
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

OUTPUT_DIR = '.'
PLOT_OUTPUT_DPI = 300

# each column is a tool, each row is a stack category

raw_data = [map(float, x.split('\t')) for x in 
'''0.209568299	0.24707894	0.2489726513	0.239152371	0.248972651	0.232262883	0.248972651	0.248972651	0.248972651	0.248972651	0.248972651	0.50853549	0.252287176	0.258928571	0.248972651	0.248972651	0.242904841	0.255769
0.185244845	0.178113423	0.1771290917	0.147830474	0.177129092	0.17326363	0.177129092	0.177129092	0.177129092	0.177129092	0.177129092	0.112309075	0.193053187	0.186011905	0.177129092	0.177129092	0.211185309	0.353846
0.099065722	0.092049017	0.0918237211	0.052976791	0.091823721	0.093602191	0.091823721	0.091823721	0.091823721	0.091823721	0.091823721	0.001796945	0.088230733	0.048809524	0.091823721	0.091823721	0.101001669	0.046154
0.000322165	0.000284981	0.0002834065	0.000504541	0.000283407	0.000497884	0.000283407	0.000283407	0.000283407	0.000283407	0.000283407	0	0.000310126	0.000297619	0.000283407	0.000283407	0.000834725	0.001923
0.110985825	0.109147905	0.1085447074	0.113017154	0.108544707	0.105302465	0.108544707	0.108544707	0.108544707	0.108544707	0.108544707	0.080862534	0.118778105	0.113988095	0.108544707	0.108544707	0.133555927	0.176923
0.004349227	0.004132231	0.0041093949	0.002018163	0.004109395	0.00373413	0.004109395	0.004109395	0.004109395	0.004109395	0.004109395	0.002695418	0.004496821	0.004315476	0.004109395	0.004109395	0.004173623	0.015385
0.014014175	0.01154175	0.0131784044	0.01160444	0.013178404	0	0.013178404	0.013178404	0.013178404	0.013178404	0.013178404	0.012578616	0	0.013839286	0.013178404	0.013178404	0	0
0.376449742	0.357651753	0.3559586226	0.432896065	0.355958623	0.391336819	0.355958623	0.355958623	0.355958623	0.355958623	0.355958623	0.281221923	0.342843852	0.373809524	0.355958623	0.355958623	0.306343907	0.157692'''.split('\n')]

data = np.array(raw_data)


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
    'PhytoCRISP-Ex',
    'sgRNAcas9',
    'sgRNAScorer2',
    'SSC',
    'TUSCAN',
    'WU-CRISPR',
    'ViCRISPR'
]

reasons = [
    'Not considered by mm10db',
    'Accepted by both',
    '* Multiple exact matches',
    #'* Multiple matches in exons',
    #'* Multiple matches in genome',
    '* Off-target score',
    '* Secondary structure or energy',
    '* Too close to reverse primer',
    '* Poly-thymine',
    '* GC-Content',
]

def stacked_bar(data, series_labels, category_labels=None, 
                show_values=False, value_format="{}", y_label=None, 
                grid=True, reverse=False):
    """Plots a stacked bar chart with the data and labels provided.

    Keyword arguments:
    data            -- 2-dimensional numpy array or nested list
                       containing data for each series in rows
    series_labels   -- list of series labels (these appear in
                       the legend)
    category_labels -- list of category labels (these appear
                       on the x-axis)
    show_values     -- If True then numeric value labels will 
                       be shown on each bar
    value_format    -- Format string for numeric value labels
                       (default is "{}")
    y_label         -- Label for y-axis (str)
    grid            -- If True display grid
    reverse         -- If True reverse the order that the
                       series are displayed (left-to-right
                       or right-to-left)
    """

    ny = len(data[0])
    ind = list(range(ny))

    axes = []
    cum_size = np.zeros(ny)

    data = np.array(data)

    if reverse:
        data = np.flip(data, axis=1)
        category_labels = reversed(category_labels)

    for i, row_data in enumerate(data):
        axes.append(plt.bar(ind, row_data, bottom=cum_size, 
                            label=series_labels[i]))
        cum_size += row_data

    if category_labels:
        plt.xticks(ind, category_labels, rotation=25)

    if y_label:
        plt.ylabel(y_label)

    

    if grid:
        plt.grid()

    if show_values:
        for axis in axes:
            for bar in axis:
                w, h = bar.get_width(), bar.get_height()
                plt.text(bar.get_x() + w/2, bar.get_y() + h/2, 
                         value_format.format(h), ha="center", 
                         va="center")

fig = plt.figure(figsize=(16,12))
fig = fig.add_subplot(111)
#fig.get_yaxis().set_major_formatter(
    #matplotlib.ticker.FuncFormatter(lambda y, p: "{0:.0%}".format(int(y))))
fig.get_yaxis().set_ticklabels(['%s%%' % z for z in xrange(0, 120, 20)])
stacked_bar(
    data, 
    reasons, 
    category_labels=tools, 
    show_values=False, 
    grid=False,
    value_format="{:.1f}",
    y_label="Distribution of Guides"
)
lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.125), ncol=2)

#fig.set_title('Distribution of mm10db Classifications on Exon-Targeting Guides')
# save to file
outFig = fig.get_figure()
outFigFileName = '%s/barStackedAcceptedGuidesRejectedByMm10dbIEEE.eps' % (OUTPUT_DIR)
outFig.savefig(outFigFileName, format='eps', dpi=PLOT_OUTPUT_DPI, bbox_extra_artists=(lgd,), bbox_inches='tight')
print 'Wrote to: %s' % outFigFileName
plt.clf()    