'''This module contains functions for calculating, plotting, and analyzing preference-space RMSD within and between 
groups of replicate amino-acid preference measurements.

Written by Mike Doud 2014'''


import mapmuts.bayesian
import mapmuts.sequtils
import mapmuts.plot
import dms_tools.file_io
import math
import random
import matplotlib
matplotlib.use('pdf')
from matplotlib.ticker import NullFormatter
import pylab
import numpy as np
import scipy
import os


def site_distance(prefs_A, prefs_B, dfunction):
    """Computes some metric representing the difference in amino-acid
    preferences at a single site (position in the protein) 
    between site preference distributions `prefs_A` and `prefs_B`.    
    
    Calling arguments:

    `prefs_A` : a dictionary of preferences *at a single site*, keyed by 
    column labels in the format of a preferences file made by `mapmuts` 
    (only the preferences for a single site should be in this dictionary 
    - see mapmuts.io.ReadEntropyAndEquilFreqs return variable structure - 
    `prefs_A` should be the value of that return variable keyed by the 
    site of interest).

    `prefs_B` : similar to `prefs_A`, but from a different set of preferences,
    eg., from another replicate experiment of the same protein, or from an 
    experiment on a homolog of the protein A.

    `dfunction` : the distance metric to be used. Currently only uses the
    Jensen-Shannon distance metric, which is the square root of the Jensen-Shannon
    divergence.
    """
    
    if dfunction == 'Jensen-Shannon':
        # here we will use the Jensen-Shannon distance metric, which is
        # the square root of the Jensen-Shannon divergence.
        # the divergence can be calculated with mapmuts.bayesian.ShannonJensenDivergence, which 
        # requires each preference distribution to be in a np.ndarray.
        p_A_list = []
        p_B_list = []
        #for aa in mapmuts.sequtils.AminoAcids()+['*']:
        for aa in mapmuts.sequtils.AminoAcids():
            p_A_list.append(prefs_A['PI_%s' % aa])
            p_B_list.append(prefs_B['PI_%s' % aa])
        p_A_array = np.array(p_A_list)
        p_B_array = np.array(p_B_list)
        return math.sqrt( mapmuts.bayesian.ShannonJensenDivergence(p_A_array, p_B_array)  )
    else:
        raise ValueError("not a valid distance metric function")


def RMSD_Within(group_A_prefs, group_B_prefs, sites, dfunction):
    """Returns a dictionary of the within-group root mean square distances (averaged evenly across both groups), keyed by site.

    Calling arguments:

    `group_A_prefs` : A list of preference dictionaries, each in the form returned by mapmuts.io.ReadEntropyAndEquilFreqs,
    corresponding to experimental replicates of a deep mutational scanning experiment.

    `group_B_prefs` : Similar to `group_A_prefs`, but for a different group of replicate experiments (eg, a different homolog).

    `sites` : A list of sites to use. Typically all the sites that were mutagenized in the deep mutational scanning experiment.

    `dfunction` : The distance function to use. Typically 'Jensen-Shannon', which uses the Jensen-Shannon distance metric.
    """

    RMSD_within = {} # keyed by sites

    for site in sites:
        RMSD_sum = 0.0 # to sum the RMSD within both groups, later to be divided by two to return average within-group RMSD
        for prefs_i in [group_A_prefs, group_B_prefs]:
            n = 0 # number of within-group comparisons made within preference group `prefs_i`
            sum_distance_squared = 0 # for within-group comparisons within preference group `prefs_i`
            for i in range(len(prefs_i)):
                for j in range(len(prefs_i)): 
                    if (j>i):
                        sum_distance_squared += ( site_distance(prefs_i[i][site], prefs_i[j][site], dfunction) )**2
                        n += 1
            RMSD_sum += math.sqrt( sum_distance_squared / float(n) )
        
        RMSD_within[site] = 0.5 * (RMSD_sum) # the even average of within-group distances of the two groups

    return RMSD_within


def RMSD_Between(group_A_prefs, group_B_prefs, sites, dfunction):
    """Returns a dictionary of average between-group root mean square distances, keyed by site.
    
    Calling arguments:

    `group_A_prefs` : A list of preference dictionaries, each in the form returned by mapmuts.io.ReadEntropyAndEquilFreqs,
    corresponding to experimental replicates of a deep mutational scanning experiment.

    `group_B_prefs` : Similar to `group_A_prefs`, but for a different group of replicate experiments (eg, a different homolog).

    `sites` : A list of sites to use. Typically all the sites that were mutagenized in the deep mutational scanning experiment.

    `dfunction` : The distance function to use. Typically 'Jensen-Shannon', which uses the Jensen-Shannon distance metric.
    """

    RMSD_between = {} # dictionary keyed by site to hold that site's average between-group squared distance
    for site in sites:
        n = 0 # total number of between-group distances calculated
        sum_distance_squared = 0
        for prefsAi in group_A_prefs:
            for prefsBi in group_B_prefs:
                sum_distance_squared += ( site_distance(prefsAi[site], prefsBi[site], dfunction) )**2
                n += 1
        RMSD_between[site] = math.sqrt( sum_distance_squared/float(n) )
    return RMSD_between


def CalculateAndWriteRMSD(group_A_prefs, group_B_prefs, sites, output, dfunction):
    """Calculates site-specific RMSD values comparing site preferences between replicates of group A and replicates of group B.
    Outputs calculations to file specified by output.

    Calling arguments:

    `group_A_prefs` : A list of preference dictionaries, each in the form returned by mapmuts.io.ReadEntropyAndEquilFreqs,
    corresponding to experimental replicates of a deep mutational scanning experiment

    `group_B_prefs` : Similar to `group_A_prefs`, but for a different group of replicate experiments (eg, a different homolog)

    `sites` : A list of sites to use. Typically all the sites that were mutagenized in the deep mutational scanning experiment.

    `output` : file to write calculations out to.

    `dfunction` : The distance function to use. Typically 'Jensen-Shannon', which uses the Jensen-Shannon distance metric.
    """

    rmsd_within = RMSD_Within(group_A_prefs, group_B_prefs, sites, dfunction)
    rmsd_between = RMSD_Between(group_A_prefs, group_B_prefs, sites, dfunction)
    
    fileout = open(output, 'w')
    fileout.write("SITE\tRMSD_BETWEEN\tRMSD_WITHIN\tRMSD_CORRECTED\n")

    for site in sites:
        fileout.write("%d\t%f\t%f\t%f\n" % (site, rmsd_between[site], rmsd_within[site], rmsd_between[site] - rmsd_within[site]) )
    
    fileout.close()


def PlotHistogram(  vals_list, label_list, plotfile, bin_size, xlabel, ylabel, 
                    rwidth=1.0, legend_loc=0, normalize=False, title='', show_legend=True, 
                    hist_type='bar', xlim=False, show_medians = False, label_medians = False, colors = False):
    '''This function can plot one or more histograms for sets of data stored in vals_list.
    Each histogram is associated with a label in label_list.

    Importantly, vals_list must be a list of lists, so even if only one list is being sent
    as an argument, for example a, send the list as [a]. Also label_list must be a list,
    not a string, even if there is only 1 label.

    Calling arguments:

    `vals_list` : List where each element is a list of values to be plotted as a histogram.

    `label_list` : List where each element is a label corresponding to an element in vals_list.

    `plotfile` : Name of file to save plot.

    `bin_size` : Width of bins to use in histogram.

    `label_list` : label for each corresponding element in vals_list.

    `xlabel` : Label of x-axis.

    `ylabel` : Label of y-axis.

    `rwidth` : Relative width of bins. Ranges from 0 to 1. If 1, bins are touching one another.

    `legend_loc` : Value can range from 0 to 10, corresponding to different locations in
    plot window. The default code 0 is considered the best location.

    `title` : Histogram title.

    `normalize` : Whether to normalize the height of the bars, so that the heights add up
    to 1 like in a true probability distribution. Set to True or False.

    `show_legend` : Whether to show a legend. Set to True or False.

    `hist_type` : Default value is 'bar', but other types are possible in matplotlib.

    `xlim` : Optional x-axis limit.

    `show_medians` : Optionally draw dashed lines to show median of distributions.

    `label_medians` : Optionally label the values of the medians.

    `colors`: Set custom color cycle.
    '''

    pylab.rc('text', usetex=True)

    assert len(vals_list) == len(label_list) # Need a label for each histogram to be plotted

    # Set up bins so they are equivalent width across different histograms
    min_val = min([min(vals) for vals in vals_list])
    max_val = max([max(vals) for vals in vals_list])
    bins = np.linspace(min_val, max_val, (max_val-min_val)/bin_size) 

    # Plot one or more histograms on a common set of axes
    fig = pylab.figure()
    ax = fig.add_subplot(111)

    if not colors:
        colors = [ 'Red','RoyalBlue','DimGray','DodgerBlue','g','b','burlywood','r','m','y']

    for i, (vals, label) in enumerate(zip(vals_list, label_list)):
        n, bins, patches = ax.hist(vals, bins, color=colors[i], rwidth = rwidth, alpha = 0.7, label = label, histtype = hist_type, normed=normalize)
        for patch in patches:
            patch.set_edgecolor('None')
        if show_medians:
            pylab.axvline(np.median(vals), color=colors[i], linestyle='dashed', linewidth=4)
            if label_medians:
                pylab.text(np.median(vals), pylab.gca().get_ylim()[1]-i, "%.3f" % np.median(vals), horizontalalignment='left', verticalalignment='top',  size=17)

    ax.set_xlabel(xlabel, fontsize = 25)
    ax.set_ylabel(ylabel, fontsize = 24)
    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(22)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(22)

    if xlim:
        pylab.gca().set_xlim(xlim)

    if title == 'Aichi/1968 previous study vs Aichi/1968 current study':
        title = 'Aichi/1968 previous study vs\nAichi/1968 current study'
    title.replace('vs','vs.')
    fig.suptitle(title, fontsize = 20)
    if show_legend:
        ax.legend(loc = legend_loc, fontsize=16)

    # remove right and top spines
    pylab.gca().spines['right'].set_color('none')
    pylab.gca().spines['top'].set_color('none')
    pylab.gca().xaxis.set_ticks_position('bottom')
    pylab.gca().yaxis.set_ticks_position('left')

    pylab.tight_layout()
    fig.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotGroupedBoxPlot(list_of_groups_of_data, xlabels, ylabel, plotfileout, title=False, y_lim=False, colors=False):
    '''Plots grouped boxplots. 

    Calling arguments:

    `list_of_groups_of_data` : a nested list. Each element of this list is an ordered group
    of measurements.  For example, the first element of this list could be a list of N experimental measurements,
    the second element could be a list of N simulated measurements, and the third element could be a list of 
    N randomized measurements. N is the number of grouped boxplots that will be drawn, and each item in 
    `list_of_groups_of_data` is drawn with a different color. The length of each item in `list_of_groups_of_data` must be the same.

    `xlabels` : list of tick labels for each group.

    `ylabel` : y-axis label.

    `plotfileout` : file to save plot to.

    `title` : Optional title.

    `y_lim` : optional y-axis limits. 

    `colors` : optional color cycle.
    '''

    pylab.rc('text', usetex=True)

    num_groups = len(list_of_groups_of_data) # number of colors
    groups_of_data = zip(*list_of_groups_of_data)
    fig = pylab.figure()
    ax = pylab.axes()
    
    if not colors:
        colors = ['LightSlateGray', 'MediumOrchid','SlateBlue','ForestGreen','Tomato','Lime', 'DimGray', 'DarkGray', 'PapayaWhip', 'RoyalBlue', 'PeachPuff' 'SlateGray', 'Silver', 'DodgerBlue', 'DeepSkyBlue', 'RoyalBlue']

    width = 0.8

    for i,group in enumerate(groups_of_data):

        pos = [ i*(num_groups+1) + p for p in range(1,num_groups+1) ]
        bp = pylab.boxplot(group, positions = pos, widths = width)
        for i_box,box in enumerate(bp['boxes']):
            pylab.setp(bp['boxes'][i_box], color='Black')
            box = bp['boxes'][i_box]
            boxX = []
            boxY = []
            for j in range(5):
                boxX.append(box.get_xdata()[j])
                boxY.append(box.get_ydata()[j])
            boxCoords = zip(boxX,boxY)
            boxPolygon = pylab.Polygon(boxCoords, facecolor=colors[i_box])
            ax.add_patch(boxPolygon)
            pylab.setp(bp['fliers'][2*i_box], marker='+', alpha=0.9, color=colors[i_box], markersize=8)
            pylab.setp(bp['fliers'][2*i_box+1], marker='+', alpha=0.9, color=colors[i_box], markersize=8)
        pylab.setp(bp['whiskers'], color='Black', linestyle='solid')

        for m in range(num_groups):
            med = bp['medians'][m]
            medianX = []
            medianY = []
            for j in range(2):
                medianX.append(med.get_xdata()[j])
                medianY.append(med.get_ydata()[j])
                pylab.plot(medianX, medianY, 'k') # , linewidth=1.5
                # medians[i] = medianY[0] # if want to label medians at top of plot

    pylab.xlim( 0, (num_groups+1)*len(groups_of_data) )

    if y_lim:
        pylab.ylim(y_lim)
    ax.set_xticklabels(xlabels, rotation=0, fontsize=15)

    tickpos = []
    for i in range(len(groups_of_data)):
        tickpos.append( i*(num_groups+1) + (num_groups+1)/2.0 )
    ax.set_xticks(tickpos)
    pylab.ylabel(ylabel)

    # plot y=0
    # pylab.plot( [0,(num_groups+1)*len(groups_of_data)], [0,0], 'k--', linewidth=0.5)

    if title:
        ax.set_title('%s' % title)

    ax.set_ylabel(ylabel, fontsize = 25)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(18)

    # remove right and top spines
    pylab.gca().spines['right'].set_color('none')
    pylab.gca().spines['top'].set_color('none')
    pylab.gca().xaxis.set_ticks_position('bottom')
    pylab.gca().yaxis.set_ticks_position('left')

    pylab.savefig("%s.pdf" % plotfileout)
    pylab.clf()
    pylab.close()


def PlotDistancesBetweenWithinWithSubgroup(dist_between, dist_within, subgroup_site_lists, total_label, subgroup_label_list, 
    plotfileout, subgroup_colors=False, enforce_limits=False, set_x_label=False, set_y_label=False):
    '''Makes a scatterplot RMSD_within vs. RMSD_between for one comparison between two groups of preferences, 
    and highlights one or more subgroups of sites.

    Calling arguments:

    `dist_between` : dictionary of RMSD_between, keyed by site, as formed by
    the function ReadFstFile.

    `dist_within` : dictionary of RMSD_within, keyed by site, as formed in 
    the function ReadFstFile. 

   `subgroup_site_lists` : A list of site subgroups.  Each element in this list is a list of sites for a different subgroup to highlight.
    Can be "none".

    `total_label` : is the label describing the comparison used to generate
    dist_between and dist_within, and subgroup_label is the label describing that subgroup of sites.
    '''
    

    # use tex
    pylab.rc('text', usetex=True)

    if not subgroup_colors:
        subgroup_colors = ['r','b','g','m','y']

    xs = []
    ys = []
    pylab.plot([0,1], c='k', alpha=0.5) # y=x line for comparison

    # Plot all sites
    for site in dist_within.keys():
        site = str(site)
        xs.append(dist_within[site])
        ys.append(dist_between[site])
    pylab.plot(xs,ys, marker='o', linestyle='None', c='DimGray', markeredgewidth=0.0, label = None, alpha = 0.4)  

    for i,subgroup_sites in enumerate(subgroup_site_lists):
        subgroup_xs = []
        subgroup_ys = []
        for site in subgroup_sites:
            site = str(site)
            subgroup_xs.append(dist_within[site])
            if enforce_limits:
                assert dist_within[site] <= 1.0 # because x and y lim are set 0-1
            subgroup_ys.append(dist_between[site])
            if enforce_limits:
                assert dist_between[site] <= 1.0 # because x and y lim are set 0-1
        pylab.plot(subgroup_xs,subgroup_ys, marker='o', linestyle='None', c=subgroup_colors[i], markeredgewidth=0.8, label = subgroup_label_list[i], alpha = 0.8)


    pylab.xlabel(r"$RMSD_{within}$")
    pylab.ylabel(r"$RMSD_{between}$")
    if set_x_label:
        pylab.xlabel(set_x_label)
    if set_y_label:
        pylab.ylabel(set_y_label)
    
    pylab.gca().yaxis.label.set_size(26)
    pylab.gca().xaxis.label.set_size(26)
    pylab.setp( pylab.gca().get_xticklabels(), fontsize=22)
    pylab.setp( pylab.gca().get_yticklabels(), fontsize=22)

    if total_label == 'Aichi/1968 previous study vs Aichi/1968 current study':
        total_label = 'Aichi/1968 previous study vs\nAichi/1968 current study'
    total_label.replace('vs','vs.')
    pylab.gca().set_title(total_label, y=1.07, fontsize = 20)

    if enforce_limits:
        pylab.gca().set_ylim([0,1])
        pylab.gca().set_xlim([0,1])
    pylab.gca().set_aspect('equal')

    # remove right and top spines
    pylab.gca().spines['right'].set_color('none')
    pylab.gca().spines['top'].set_color('none')
    pylab.gca().xaxis.set_ticks_position('bottom')
    pylab.gca().yaxis.set_ticks_position('left')

    pylab.legend(loc = 'lower right', prop={'size':16}, numpoints=1 )

    pylab.tight_layout()
    pylab.show()
    pylab.savefig("%s.pdf" % plotfileout)
    pylab.close()


def PlotAAPrefCorrelations(prefs_A, prefs_B, sites, xlabel, ylabel, plotfileout, read_from_file=False):
    '''Plots correlation between all amino-acid preferences at the indicated sites between two sets of preferences.

    Calling arguments:

    `prefs_A` : nested dictionary of amino acid preferences, keyed by site and then keyed by PI_X where X = amino acid. Can be
    a file of the format produced by mapmuts when using `read_from_file` = True.

    `prefs_B` : same as `prefs_A`.

    `sites` : a list of sites that contain preferences in prefs_A and prefs_B. Only these sites will be used.

    `read_from_file` : Set this to True to read in preferences from mapmuts-formatted preference files.

    `xlabel` : TeX-formatted label for x-axis (Should describe the preferences given as `prefs_A`).

    `ylabel` : TeX-formatted label for y-axis (Should describe the preferences given as `prefs_B`).

    `plotfileout` : Filename to save the plot to. Doesn't need to include an extension.
    '''


    if read_from_file:
        file_a = prefs_A
        file_b = prefs_B

        prefs_A = ReadDMSToolsFormattedPrefsToMMFormattedDict(file_a)
        prefs_B = ReadDMSToolsFormattedPrefsToMMFormattedDict(file_b)

    pylab.rc('text', usetex=True)
    
    xs = []
    ys = []

    for site in sites:
        for aa in mapmuts.sequtils.AminoAcids():
            xs.append(prefs_A[site]['PI_%s' % aa])
            ys.append(prefs_B[site]['PI_%s' % aa])

    pylab.plot(xs,ys, marker='o', linestyle='None', c='Black', markeredgewidth=0.0, alpha = 0.2, markersize=6)

    (r, p) = scipy.stats.pearsonr(xs, ys)
    r = '$R = %.2f$' % r
    if p < 1e-10:
        p = '$P < 10^{-10}$'
    else:
        p = '$P = %s$' % mapmuts.plot.Base10Formatter(p, 2, 1, 2)

    text = '%s\n%s' % (r, p)
    pylab.text(0.05, 0.96, text, horizontalalignment='left', verticalalignment='top',  size=25) #transform=ax.transAxes

    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)

    pylab.gca().yaxis.label.set_size(26)
    pylab.gca().xaxis.label.set_size(26)
    pylab.setp( pylab.gca().get_xticklabels(), fontsize=23)
    pylab.setp( pylab.gca().get_yticklabels(), fontsize=23)

    pylab.gca().set_ylim([-0.02,1])
    pylab.gca().set_xlim([-0.02,1])
    pylab.gca().set_aspect('equal')

    # remove right and top spines
    pylab.gca().spines['right'].set_color('none')
    pylab.gca().spines['top'].set_color('none')
    pylab.gca().xaxis.set_ticks_position('bottom')
    pylab.gca().yaxis.set_ticks_position('left')

    yticker = matplotlib.ticker.FixedLocator([0, 0.5, 1])
    xticker = matplotlib.ticker.FixedLocator([0, 0.5, 1])
    pylab.gca().yaxis.set_major_locator(yticker)
    pylab.gca().xaxis.set_major_locator(xticker)

    pylab.tight_layout()
    pylab.show()
    pylab.savefig("%s.pdf" % plotfileout)
    pylab.close()


def AveragePairwisePearsonsR(list_of_preference_files, sites, read_from_file = True):
    '''Returns the average Pearson's R between all pairwise comparisons between individual amino-acid 
    preference measurements.

    Calling arguments:

    `list_of_preference_files` : a list of preference files. Optionally, this can be a list of dictionaries of preferences 
    when `read_from_file` is set to False.

    `sites` : a list of sites to be used to calculate correlations.  Typically all the sites should be used here.

    `read_from_file` : set this to False if `list_of_preference_files` is actually a list of nested dictionaries.
    '''
    if read_from_file:
        prefs_dicts = [ ReadDMSToolsFormattedPrefsToMMFormattedDict(f) for f in list_of_preference_files ]
    else:
        prefs_dicts = list_of_preference_files
    aas = mapmuts.sequtils.AminoAcids(includestop=False)
    sum_of_r = 0.0
    num_of_r = 0
    for i in range(len(prefs_dicts)-1):
        for j in range(i+1, len(prefs_dicts)):
            sample1data = []
            sample2data = []
            for r in sites:
                for aa in aas:
                    sample1data.append(prefs_dicts[i][r]['PI_%s' % aa])
                    sample2data.append(prefs_dicts[j][r]['PI_%s' % aa])
            (pearson_r, p) = scipy.stats.pearsonr(sample1data, sample2data)
            sum_of_r += pearson_r
            num_of_r += 1
    average_r = sum_of_r / num_of_r
    return average_r


def ReadRMSDFile(filename):
    '''Reads the RMSD_between, RMSD_within, and RMSD_corrected from an RMSD file
    of the format written by CalculateAndWriteRMSD. Returns a 3-tuple of dictionaries
    which are each keyed by site: (rmsd_between_d, rmsd_within_d, rmsd_corrected_d).
    '''
    rmsd_between_d = {} #rmsd_between
    rmsd_within_d = {} #rmsd_within
    rmsd_corrected_d = {} #corrected_rmsd
    lines = (line.rstrip('\n') for line in open(filename) if line[0:4] != 'SITE')
    for line in lines:
        
        result = line.split('\t')
        site = result[0]
        rmsd_between = result[1] 
        rmsd_within = result[2] 
        rmsd_corrected = result[3]

        rmsd_between_d[site] = float(rmsd_between)
        rmsd_within_d[site] = float(rmsd_within)
        rmsd_corrected_d[site] = float(rmsd_corrected)

    return (rmsd_between_d, rmsd_within_d, rmsd_corrected_d)



def DrawDirichletEquilFreqs(preferences, aas, scale_conc = 1):
    '''This function draws equilibrium preferences for each amino acid 
    at each codon site from a Dirichlet distribution. The alpha parameter 
    for each category (amino acid) at each site is simply set to the inferred 
    equilibrium preference for that amino acid at the site. The function
    returns the set of equilibrium preferences sampled from a Dirichlet.

    `preferences` : This dictionary contains all the amino-acid preferences 
    for a set of codon sites. The  dictionary is keyed by residue numbers, 
    and with the values for each residue number being dictionaries keyed by 
    all the other column labels ('WT_AA', 'SITE_ENTROPY', 'PI_A', 'PI_C', ...)
    and possibly 'PI_*' if that is present in the header.
    
    `aas` : A list of the amino-acids found in the preferences dictionary. This
    list may include stop codon (*).

    `scale_conc` : A scalar that scales the alpha concentration parameters in
    the Dirichlet distribution. This scalar is uniformly applied to all
    the alpha parameters at each codon site. This scalar should be positive.
    The default value is 1.
    
    `alphas` : A dictionary where each key is a residue number for a site, 
    and the value is a list of the alpha parameters for the site. The alpha
    parameters for each amino-acid follow the same ordering as in the 
    aas list. The alpha parameters at a site sum to 1.
    
    `sites` : The list of residue numbers that make up the keys of the 
    preferences dictionary.
    
    `sample` : A new set of equilibrium preferences drawn from a Dirichlet 
    distribution. This is a dictionary with the same structure as the
    preferences dictionary.

    Written by Orr Ashenberg.  

    Extended description
    --------------------

    When inferring amino-acid preferences from deep mutational scanning
    data, there is noise both from the experiments and from the computational
    analyses. This noise results in less than perfect correlation of amino-acid
    preferences between replicates of the same experimental libraries. We
    would like to distinguish biologically significant differences in amino-acid
    preferences, from differences that solely appear due to noise in the
    measurement. To do this, we need a model that describes the noise in the
    measurements. We can then use such a model to simulate the differences 
    we expect to see in replicates of the same libraries. If comparisons between
    libraries from the real, experimental data show differences greater than
    what is seen in these simulations, that will indicate those differences
    have biological meaning.

    The equilibrium preferences are a form of compositional data, as the preferences
    at a codon site must always sum to 1. One way to model sampling variance in
    compositional data is through a Dirichlet distribution. For a given codon site,
    the support of the distribution is the vector of 20 possible equilibrium preferences at 
    the site. The parameters of the distribution are the number of categories (20 amino-acids)
    the random variable can fall into, and the concentration parameter for each category.

    The concentration parameter (alpha) sets the mean and variance of the Dirichlet
    distribution. In particular, the mean of random variable X_i is alpha_i / sum(alpha_i)
    where i ranges over the number of categories. Initially, the alpha parameters at
    each site are set to the experimental amino-acid preferences. For example, if
    alanine occurs at position 5 with frequency 0.07, then alpha_alanine at position
    5 is set to 0.07. Since the amino-acid frequencies at a site sum to 1, this ensures
    that each random variable (category) has a mean that is the experimental amino-acid fequency.

    By multiplying all the alpha parameters at each site by the same scalar, the variance
    of each amino-acid preference can be changed, without changing the mean. If
    alpha = 1, the Dirichlet distribution reduces to a uniform distribution. If alpha<1,
    only a few categories at each site will have nonzero probability. In other words, 
    the mass of the distribution will be highly concentrated on a few categories. For example,
    alanine at a site may be assigned frequency 0.99, while all residues will be assigned
    frequencies < 0.01. If alpha>1, most categories for each site will be sampled. In other words,
    the mass will be dispersed more equally among all the categories. For example, alanine at a site may
    be assigned frequency 0.1, while 10 other amino-acids at the site may be assigned frequency
    0.05. To visualize this behavior, you can plot x^(alpha-1)*(1-x)^(alpha-1) for x in (0,1) as
    you vary alpha from 0 to infinity. This function is a Dirichlet distribution with 2 categories, 
    or more simply a Beta distribution.

    As alpha increases beyond value 1, the sampling for each category stays closer and closer to 
    the mean value of the category, ie the variance decreases. So increasing alpha, by multiplying
    the vector of alpha parameters at a site by a scalar > 1, results in samples that will
    look more similar to the experimental amino-acid preferences, which were used to set up
    the vector of alpha parameters. This also results in higher correlations between replicates
    of the same library.
    '''

    if scale_conc < 0:
        raise ValueError('scale_conc value must be positive')

    alphas = {}
    sample = {}
    sites = sorted(preferences.keys())

    # Fix concentration parameters for amino acid at each site 
    # to the corresponding amino-acid preference
    for r in sites:
        alphas[r] = [preferences[r]['PI_%s' % aa] for aa in aas]
        #print r, alphas[r]

    # Iterate through each site, and sample amino-acid preferences 
    # at the site using a Dirichlet distribution.
    for r in sites:
        # Scale alpha concentration parameters
        alpha = [i*scale_conc for i in alphas[r]]

        # Draw from Dirichlet distribution and save draw in sample[r]
        # threshold enforces individual elements to be greater than a value to prevent 
        # divide by zero errors later:
        threshold = 1e-180
        
        draw = np.random.dirichlet(alpha)
        while( sum( [i<threshold for i in draw] ) > 0):
            draw = np.random.dirichlet(alpha)
        sample[r] = dict( ('PI_%s' % aas[i], draw[i]) for i in range(len(aas)) )
        sample[r]['SITE_ENTROPY'] = mapmuts.bayesian.SiteEntropy(sample[r])
        sample[r]['WT_AA'] = preferences[r]['WT_AA']
            
    return sample


def WriteEquilFreqs(sample, aas, sites, filename):
    '''Write out file containing amino-acid preferences for a 
    set of codon sites.

    `sample` : This dictionary contains all the amino-acid preferences 
    for a set of codon sites. The  dictionary is keyed by residue numbers, 
    and with the values for each residue number being dictionaries keyed by 
    all the other column labels ('WT_AA', 'SITE_ENTROPY', 'PI_A', 'PI_C', ...)
    and possibly 'PI_*' if that is present in the header.
    
    `aas` : A list of the amino-acids found in the preferences dictionary. This
    list may include stop codon (*).
    
    `sites` : The list of residue numbers that make up the keys of the 
    sample dictionary.
    
    `filename` : Name of file for writing equilibrium preferences.

    Written by Orr Ashenberg.
    '''

    with open (filename, 'w') as equilibriumpreferencesfile:
        equilibriumpreferencesfile.write('#SITE\tWT_AA\tSITE_ENTROPY')
        for aa in aas:
            equilibriumpreferencesfile.write('\tPI_%s' % aa)
        equilibriumpreferencesfile.write('\n')
        for r in sites:
            h = sample[r]['SITE_ENTROPY']
            equilibriumpreferencesfile.write('%d\t%s\t%g' % (r, sample[r]['WT_AA'], h))
            for aa in aas:
                equilibriumpreferencesfile.write('\t%g' % sample[r]['PI_%s' % aa])
            equilibriumpreferencesfile.write('\n')

def ReadDMSToolsFormattedPrefsToMMFormattedDict(f):
    '''Reads a preferences file of the format generated by dms_tools (which cannot be read by mapmuts.io.ReadEntropyAndEquilFreqs())
    and returns preferences in a dictionary of the format made by mapmuts.io.ReadEntropyAndEquilFreqs(), where preferences are
    keyed by [site]['PI_X'], wild-type amino-acid is keyed by [site]['WT_AA'], and site entropies are keyed by [site]['SITE_ENTROPY']
    '''
    site_strings, wts, prefs, pi95, h = dms_tools.file_io.ReadPreferences(f)
    newmmprefs_dict = {}
    for site in site_strings:
        newmmprefs_dict[int(site)] = {}
        for aa in mapmuts.sequtils.AminoAcids():
            newmmprefs_dict[int(site)]['PI_%s'%aa] = prefs[site][aa]
        newmmprefs_dict[int(site)]['WT_AA'] = wts[site]
        newmmprefs_dict[int(site)]['SITE_ENTROPY'] = h[site]
    return newmmprefs_dict

def WriteReplicateMeasurementPrefs( list_of_preffiles, sites, aas=mapmuts.sequtils.AminoAcids(includestop=False), outputdir = './', fileprefix='' ):
    '''Creates an amino-acid preference file containing the replicate measurements in `list_of_preffiles` at each site in `sites`.
    The site indices in the output file will correspond to the replicate preference measurements in `list_of_preffiles`.

    This function was used to make preference files that were converted to logoplots for Figure 3 of the paper.
    '''
    
    pref_dists = [ ReadDMSToolsFormattedPrefsToMMFormattedDict(f) for f in list_of_preffiles ]

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    os.chdir(outputdir)

    for site in sites:
        newprefs = {} # will have "site indices" that actually correspond to replicate experimental measurements.
        for i,pref_dist in enumerate(pref_dists): # for each replicate `i` measurement of preferences `pref_dist`
            newprefs[i+1] = pref_dist[site]

        # write out a preferences file containing just this site's preferences from every replicate:
        fake_sites = range(1, len(pref_dists)+1 )
        filename = "%s_site_%s_replicatemeasurements.txt" % (fileprefix,site)
        WriteEquilFreqs(newprefs, aas, fake_sites, filename)


        


