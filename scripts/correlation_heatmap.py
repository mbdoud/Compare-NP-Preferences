'''This script calculates Pearson correlation coefficients between all pairwise comparisons
for a given set of amino-acid preference files. Makes a heatmap colored by 
the correlation coefficients for those pairwise comparisons. The correlation is 
only calculated for the set of sites provided in the infile, so the script can
be run multiple times with different subsets of sites to see how the correlations between
particular amino-acid preference files change for particular sites in the protein.

This script was written by Mike Doud, September 2014.

This script is run with a single argument, the name of an infile specifying the following key/value pairs:

preference_files:  space-separated list of preference files to be used for comparisons
sample_names: space-separated list of labels corresponding, in order, to the preference_file_list
sites_type:  either range or list, as described below.
sites: sites used for comparison.  if site_type range, give two numbers, eg. `2 498`; if site_type list, use comma and space-separated site numbers, eg. `2, 3, 4, 5, 6, 7, 19, 41`
outputprefix:  file prefix for output files
includestop:  True or False
'''


import sys
import os
import mapmuts.io
import mapmuts.plot
import mapmuts.sequtils
import scipy.stats
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy as np
import rmsdtools

def main():
    args = sys.argv[1:]
    if len(args) != 1:
        raise IOError("must call with one argument as input file")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))

    # infile parsing:
    file_list = mapmuts.io.ParseFileList(d, 'preference_files' )
    outputprefix = mapmuts.io.ParseStringValue(d, 'outputprefix')
    samples =  mapmuts.io.ParseStringValue(d, 'sample_names').split()
    sample_names = []
    for sample in samples:
        sample_names.append(sample)
    assert len(file_list) == len(sample_names)

    prefs_dict = {} # holds preference distributions keyed by sample_name
    for i,preffile in enumerate(file_list):
        prefs_dict[sample_names[i]] = rmsdtools.ReadDMSToolsFormattedPrefsToMMFormattedDict(preffile)
    sites_type = mapmuts.io.ParseStringValue(d, 'sites_type')

    if sites_type == 'range':
        sites = mapmuts.io.ParseStringValue(d, 'sites')
        (start,end) = sites.split()
        (start,end) = (int(start),int(end))
        sites = [r for r in range(start,end+1)]

    elif sites_type == 'list':
        sites = mapmuts.io.ParseStringValue(d, 'sites')
        sites = sites.split(', ')
        sites = [int(r) for r in sites]

    #print sites 
    includestop = mapmuts.io.ParseBoolValue(d, 'includestop')

    aas = mapmuts.sequtils.AminoAcids(includestop=includestop)

    correlations = {} # keyed by 'sample1,sample2'

    correlation_data = np.zeros(( len(sample_names), len(sample_names) ))

    outputfile = outputprefix + "-correlations.txt"
    fileout = open(outputfile, 'w')
    fileout.write("Sample\tSample\tR\n")

    for i,sample1 in enumerate(sample_names):
        for j,sample2 in enumerate(sample_names): # in sample_names[i:] or in sample_names[i+1:]
            sample1data = []
            sample2data = []
            for site in sites:
                for aa in aas:
                    sample1data.append(prefs_dict[sample1][site]['PI_%s' % aa])
                    sample2data.append(prefs_dict[sample2][site]['PI_%s' % aa])
            (this_r, this_p) = scipy.stats.pearsonr(sample1data, sample2data)
            correlations[ "%s,%s" % (sample1, sample2)]  = this_r
            fileout.write("%s\t%s\t%f\n" % (sample1, sample2, this_r))
            correlation_data[i][j] = this_r

    fileout.close()

    plotHeatmap(correlation_data,sample_names,outputprefix)


def plotHeatmap(data,labels,outputprefix):
    column_labels = labels
    row_labels = labels
    fig, ax = pylab.subplots()
    heatmap = ax.pcolor(data, cmap=pylab.cm.Greens, vmin=0, vmax=1)
    cbar = pylab.colorbar(heatmap)
    cbar.solids.set_edgecolor("face")
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.set_xticklabels(row_labels, minor=False, rotation='vertical')
    ax.set_yticklabels(column_labels, minor=False)
    pylab.tight_layout()
    pylab.gca().set_aspect(1)
    pylab.show()
    pylab.savefig("%s.pdf" % outputprefix)

main()


