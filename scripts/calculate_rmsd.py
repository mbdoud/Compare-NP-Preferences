"""This script quantifies site-specific differences in amino-acid preferences between two groups of replicate
deep mutational scanning experiments. Each group is a set of replicate measurements of amino-acid preferences
inferred using ``dms_tools``. The script quantifies site-specific noise using replicate experiments, and subtracts
this level of noise from the differences observed between groups.

For example, this script can be used to quantify the differences in amino-acid preferences between the **PR/1934** 
and **Aichi/1968** strains of influenza nucleoprotein. For that purpose, each group will specify a set of amino-acid
preferences measured in biological replicate experiments for each strain.

Since there is noise in the experiment at the level of individual replicates, this analysis quantifies the average differences 
in amino-acid preferences *between* the two groups, as well as *within* each group. The average *within* group differences quantify 
the noise between biological replicates, while the *between* group differences contain this noise as well as any true difference in 
amino-acid preferences between the two groups. This script calculates, for all specified sites, the average difference *between* 
the two groups and *within* the two groups. The differences are calculated using a distance function of the user's choosing, 
as long as the distance function is encoded in rmsdtools.site_distance(). Currently, only the Jensen-Shannon distance metric is coded 
for this purpose.

The average distance *between* groups is calculated over all unique pairwise comparisons of replicates from different 
groups, and the average distance *within* groups is taken from all unique pairwise comparisions of replicates among the same group (for both groups).

The script outputs a tab-delimited file containing, for each site:

    * the root mean square distance between replicates in different groups (RMSD_between)

    * the root mean square distance between replicates within the same groups (RMSD_within)

    * the corrected root mean square distance (RMSD_corrected = RMSD_between - RMSD_within) 

This script was written by Mike Doud.

This script must be run with one argument, the name of an infile consiting of key/value pairs for the following options:

    * `number_of_groups` : the number of groups of preferences to compare. Must be at least 2.

    * `comparisons` : list which group numbers to make comparisons between. Each comparison consists of two group numbers 
    separated by a comma, and multiple comparisons are separated by spaces. If only making comparisons between two groups,
    the value for this key should be `1,2`. To perform all possible comparisons between three specified groups, the value
    for this key could be `1,2 1,3 2,3`.

    * `sites_type` : Desribes how the sites to be used for calculates are listed. Either 'range' or 'list', as described below.

    * `sites` : If providing a range, as declared by `sites_type`, the range of sites to compare should be specified by two numbers, inclusive.
    Otherwise, if providing a list of individual sites, should provide a comma and space between each site, eg. `2, 3, 4, 10, 20, 21`

    * `group_1_files` : a space-separated list of preference files for group 1, with multiple groups. Each file of preferences
    should be in the same form that is written out by ``dms_tools``, with stop codon preferences removed before running this script.

    * `group_2_files` : same as `group_1_files`.

    * `group_3_files` : if comparing more than 2 groups, continue with the pattern by listing each group of preference files.

    * `group_1_name` : a descriptor of group 1 to be used in the file name of the output file.

    * `group_2_name` : a descriptor of group 2 to be used in the file name of the output file.

    * `group_3_name` : if comparing more than 2 groups, continue with the pattern by listing the name for each group.

    * `outfileprefix` : a string used as the prefix for the output file.

    * `distance_function` : 'Jensen-Shannon' is the only supported distance function. Other distance functions can be encoded in rmsdtools.site_distance().
"""

import sys
import os
import mapmuts.io
import rmsdtools

def main(): 
    
    # infile parsing:
    args = sys.argv[1:]
    if len(args) != 1:
        raise IOError("must call with one argument as input file")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    
    number_of_groups = int( mapmuts.io.ParseFloatValue(d, 'number_of_groups') )
    assert number_of_groups > 1  
    outfileprefix = mapmuts.io.ParseStringValue(d, 'outfileprefix')
    distance_function = mapmuts.io.ParseStringValue(d, 'distance_function')

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

    comparisons = mapmuts.io.ParseStringValue(d, 'comparisons')
    comparisons = comparisons.split()
    comps = []
    for comparison in comparisons:
        a,b = comparison.split(',')
        comps.append([int(a),int(b)])
    comparisons = comps # a list of comparisons to make, each comparison designated by two group numbers to compare to each other

    groups_of_files = [] # each element in this list is a file list for a particular group. 
    # eg., groups_of_files[0] is list of files for group 1 replicates, etc...
    for i in range(1, number_of_groups + 1):
        groups_of_files.append( mapmuts.io.ParseFileList(d, 'group_%d_files' % i) )
    group_names = [] # each element in this list is the name of a group, in group 1, 2, 3, 4 order
    for i in range(1, number_of_groups + 1):
        group_names.append( mapmuts.io.ParseStringValue(d, 'group_%d_name' % i) )

    # read in the preferences from each group's file list.
    groups_of_preferences = []  # each element in this list is a list of preference dicts corresponding to a single group
    # eg., groups_of_preferences[0] is a list of preference dicts for group 1 replicates, etc...
    for i in range(number_of_groups):
        groups_of_preferences.append( [ rmsdtools.ReadDMSToolsFormattedPrefsToMMFormattedDict(f) for f in groups_of_files[i] ] )

    # Check to make sure the desired sites have preferences listed and
    # that stop codons have been properly removed from preference files:
    for group in groups_of_preferences:
        for pref_dist in group:
            for site in sites:
                assert site in pref_dist.iterkeys()
            for r in pref_dist.iterkeys():
                assert 'PI_*' not in pref_dist[r].keys()

    # Calculate RMSD for each pair of groups specified:
    for comparison in comparisons:
        i,j = comparison # i, j are group numbers
        label = "%s_vs_%s" % (group_names[i-1], group_names[j-1])
        rmsdfileprefix = "%s_%s_RMSD_calcs.txt"  %  ( outfileprefix, label.replace("/", "").replace(" ", "_") )
        rmsdtools.CalculateAndWriteRMSD(groups_of_preferences[i-1], groups_of_preferences[j-1], sites, rmsdfileprefix, distance_function)

main()
