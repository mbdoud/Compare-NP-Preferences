"""Master Python script that runs the phylogenetic analysis in this directory.

This script uses `phyloExpCM` to fit phylogenetic trees using codon substitution models. The codon substitution models are
either site specific and experimentally informed by nucleoprotein amino-acid preferences, or they are non site-specific and use
the traditional Goldman-Yang 1994 model. After fitting the phylogenetic trees, the script bins sites based on how much 
they improve the likelihood in going from the Aichi1968 model to the Aichi1968_PR1934 model, and checks whether this binning correlates
with per-site RMSD_within or RMSD_corrected. The main inputs to this script are the phylogenetic trees and the amino-acid preferences.


This script should be run with one option at the command line which specifies the path of a configuration file, eg:

python run_phyloanalysis.py phylo_config.txt

Each line of the configuration file should have a key and a value, separated by a space, for the following keys:

base_dir: This is the top directory when pulling the Compare-NP-Preferences Github repository. Important subdirectories in
          this top directory are scripts/, dmstools_output/, phylo_input/, phylo_output/, and compare_prefs_output/.


Here are the required input files that are included in the Compare-NP-Preferences Github repository.
scripts/
   phylo_config.txt

phylo_input/
   Aligned_NPs_Allhosts.fasta  
   Aligned_NPs_Human.fasta   
   Aligned_NPs_Swine.fasta   
   Aligned_NPs_Equine.fasta  
   Aligned_NPs_Avian.fasta     
   codonphyml_Allhosts_tree.newick  
   codonphyml_Human_tree.newick
   codonphyml_Swine_tree.newick
   codonphyml_Equine_tree.newick
   codonphyml_Avian_tree.newick     

dmstools_output/
   mean_Aichi68_both_studies_prefs.txt
   PR8_mean_prefs.txt
   mean_NP_both_homologs_prefs.txt

compare_prefs_output/
    Aichi1968_vs_PR1934/RMSD_Aichi1968_vs_PR1934_RMSD_calcs.txt

Here are required programs that the user must install beforehand.
    python 2.7
    mapmuts (https://github.com/jbloom/mapmuts)
    phyloExpCM (https://github.com/jbloom/phyloExpCM)
       Tested on PhyloExpCM version v0.32 from September 1, 2014 
    sbatch job submission (https://computing.llnl.gov/linux/slurm/sbatch.html)
       sbatch is needed for running phylogenetic calculation as submitted jobs to a computer cluster.
       Uses sbatch and the Python multiprocessing module to run some of
       the analyses in parallel for faster completion. This code has been tested 
       on the FHCRC's computing core. 

Written by Orr Ashenberg, March 2015.
"""


import os
import sys
import re
import string
import time
import math
import operator
import rmsdtools
import pylab
import numpy as np
import multiprocessing
import mapmuts.sequtils
import mapmuts.io
import phyloExpCM.hyphy

def RunScript(rundir, run_name, script_name, commands, use_sbatch, sbatch_cpus, walltime=None):
    """Runs a ``phyloExp_CM`` script.

    *rundir* is the directory in which we run the job. Created if it does
    not exist.

    *run_name* is the name of the run, which should be a string without
    spaces. The input file has this prefix followed by ``_infile.txt``.

    *script_name* is the name of the script that we run.

    *commands* contains the commands written to the input file. It is a list 
    of 2-tuples giving the key / value pairs.
    Both keys and values should be strings.

    *use_sbatch* is a Boolean switch specifying whether we use ``sbatch``
    to run the script. If *False*, the script is just run with the command
    line instruction. If *True*, then ``sbatch`` is used, and the command file
    has the prefix *run_name* followed by the suffix ``.sbatch``.

    *sbatch_cpus* is an option that is only meaningful if *use_sbatch* is 
    *True*. It gives the integer number of CPUs that are claimed via
    ``sbatch`` using the option ``sbatch -c``. 

    *walltime* is an option that is only meaningful if *use_sbatch* is
    *True*. If so, it should be an integer giving the number of hours 
    to allocate for the job. If *walltime* has its default value of 
    *None*, no wall time for the job is specified.

    It is assumed that the script can be run at the command line using::

        script_name infile

    Returns *runfailed*: *True* if run failed, and *False* otherwise.
    """
    print "Running %s for %s in directory %s..." % (script_name, run_name, rundir)
    currdir = os.getcwd()
    if not os.path.isdir(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)
    if (not run_name) or not all([x not in string.whitespace for x in run_name]):
        raise ValueError("Invalid run_name of %s" % run_name)
    infile = '%s_infile.txt' % run_name
    open(infile, 'w').write('# input file for running script %s for %s\n%s' % (script_name, run_name, '\n'.join(['%s %s' % (key, value) for (key, value) in commands])))
    if use_sbatch:
        sbatchfile = '%s.sbatch' % run_name # sbatch command file
        jobidfile = 'sbatch_%s_jobid' % run_name # holds sbatch job id
        jobstatusfile = 'sbatch_%s_jobstatus' % run_name # holds sbatch job status
        joberrorsfile = 'sbatch_%s_errors' % run_name # holds sbatch job errors
        sbatch_f = open(sbatchfile, 'w')
        sbatch_f.write('#!/bin/sh\n#SBATCH\n')
        if walltime:
            sbatch_f.write('#PBS -l walltime=%d:00:00\n' % walltime)
        sbatch_f.write('%s %s' % (script_name, infile))
        sbatch_f.close()
        os.system('sbatch -c %d -e %s %s > %s' % (sbatch_cpus, joberrorsfile, sbatchfile, jobidfile))
        time.sleep(60) # short 1 minute delay
        jobid = int(open(jobidfile).read().split()[-1])
        nslurmfails = 0
        while True:
            time.sleep(60) # delay 1 minute
            returncode = os.system('squeue -j %d > %s' % (jobid, jobstatusfile))
            if returncode != 0:
                nslurmfails += 1
                if nslurmfails > 180: # error if squeue fails at least 180 consecutive times
                    raise ValueError("squeue is continually failing, which means that slurm is not working on your system. Note that although this script has crashed, many of the jobs submitted via slurm may still be running. You'll want to monitor (squeue) or kill them (scancel) -- unfortunately you can't do that until slurm starts working again.")
                continue # we got an error while trying to run squeue
            nslurmfails = 0
            lines = open(jobstatusfile).readlines()
            if len(lines) < 2:
                break # no longer in slurm queue
        errors = open(joberrorsfile).read().strip()
    else:
        errors = os.system('%s %s' % (script_name, infile))
    os.chdir(currdir)
    if errors:
        print "ERROR running %s for %s in directory %s." % (script_name, run_name, rundir)
        return True
    else:
        print "Successfully completed running %s for %s in directory %s." % (script_name, run_name, rundir)
        return False


def RunProcesses(processes, nmultiruns):
    """Runs a list *multiprocessing.Process* processes.

    *processes* is a list of *multiprocessing.Process* objects that
    have not yet been started. If an empty list, then just returns 
    with no action taken.

    *nmultiruns* is an integer >= 1 indicating the number of simultaneous
    processes to run.

    Runs the processes in *processes*, making sure to never have more than
    *nmultiruns* running at a time. If any of the processes fail (return
    an exitcode with a boolean value other than *False*), an exception
    is raised immediately. Otherwise, this function finishes when all
    processes have completed.
    """
    if not processes:
        return
    if not (nmultiruns >= 1 and isinstance(nmultiruns, int)):
        raise ValueError("nmultiruns must be an integer >= 1")
    processes_started = [False] * len(processes)
    processes_running = [False] * len(processes)
    processes_finished = [False] * len(processes)
    while not all(processes_finished):
        if (processes_running.count(True) < nmultiruns) and not all(processes_started):
            i = processes_started.index(False)
            processes[i].start()
            processes_started[i] = True
            processes_running[i] = True
        for i in range(len(processes)):
            if processes_running[i]:
                if not processes[i].is_alive():
                    processes_running[i] = False
                    processes_finished[i] = True
                    if processes[i].exitcode:
                        raise IOError("One of the processes failed to complete.")
        time.sleep(60)

def FormatModelName(name):
    """Returns a LaTex formatted model name."""
    formatted_names = {
            'Aichi1968_fitbetaHalpernBruno':'Aichi/1968',
            'PR1934_fitbetaHalpernBruno':'PR/1934',
            'Aichi1968_PR1934_fitbetaHalpernBruno':'Aichi/1968 + PR/1934',
            'GY94_CF3x4_omega-global-gamma4_rates-gamma4':'GY94, gamma $\omega$, gamma rates'
                      }
    if name in formatted_names:
        return formatted_names[name]
    else:
        return name.replace('_', ' ')


def ReadSiteLikelihoodFile(infile):
    '''Read an `infile` containing delta per-site log likelihoods generated 
    by  phyloExpCM_SiteLikelihoodComparison. This file contains all the delta
    per-site likelihoods across all sites in a gene. The residue site and 
    likelihood value are stored as key/value pairs in the `sitelikelihood` 
    dictionary.

    `infile`:  Name of per-site likelihood file to read.
    `sitelikelihood`: Dictionary where key is residue site number, and value 
       is the delta of log-likelihood values in the second column of the file.
       The delta log-likelihood is calculated between two different substitution 
       models.
    '''
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile %s" % infile)

    lines = open(infile).readlines()[1 : ]
    sitelikelihood = dict([(int(line.split()[0]), float(line.split()[1])) for line in lines])
    return sitelikelihood


def BinValuesByPercentileBoxPlot(xs, ys, plotfile, xlabel, ylabel, percent):
    '''This function takes in a list of (x, y) pairs. It bins the x-values into percentiles, and 
    for each percentile bin, makes a box and whisker plot of the corresponding y-values 
    found in that bin. 
    For instance, the x-values could be sorted and then divided into 10 bins of equal size.Then a box 
    and whisker plot of the corresponding y-values for each bin would be plotted. 

    `xs`: List of x-values
    `ys`: List of y-values corresponding to x-values. Same length as xs.
    `plotfile`: Name of file to save plot
    `xlabel`: Label of x-axis
    `ylabel`: Label of y-axis
    `percent`: Size of bins, in percent, to divide xs data into. The total number of bins
       is 100/percent. For example, if percent = 10, then there will be 10 bins corresponding
       to the first 10% of sorted xs, next 10% etc... to the last 10%.
    '''
    if len(xs) != len(ys):
        raise ValueError('xs and ys list must have same length')
    if (100 % percent):
        raise ValueError('Requested percentile %s in BinValuesByPercentile does not evenly divide into 100' % percent)

    # Find the set of indices in xs list (percentile_index) that correspond to dividing the 
    # sorted xs data into percentile bins. Example: If dividing 498 sites into quintiles,
    # we would use the following indices [0, 99, 198, 297, 396, 496] to designate the quintiles.
    percentiles = range(0, 101, percent)
    index = range(len(xs))
    percentile_index = [ int(math.floor(np.percentile(index, p))) for p in percentiles]

    # Sort paired xs and ys data by the x value
    sortedxy = sorted(zip(xs, ys), key = operator.itemgetter(0))

    # Get corresponding ys in each bin of xs
    ybins = []
    # xs and ys in first percentile bin
    percentile = sortedxy[:percentile_index[1]] 
    ybins.append( [y for (x, y) in percentile] )
    # xs and ys in middle bins
    for i in range(1, len(percentile_index) - 2):
        percentile = sortedxy[percentile_index[i]:percentile_index[i + 1]]
        ybins.append( [y for (x, y) in percentile] )
    # xs and ys in last percentile bin
    percentile = sortedxy[percentile_index[-2]:] 
    ybins.append( [y for (x, y) in percentile] )

    # Set up positions for percentiles
    xpos = percentiles[1:]
    
    # Set up box and whisker plot
    list_of_groups_of_data = [ybins]
    num_groups = len(list_of_groups_of_data) # number of colors
    groups_of_data = zip(*list_of_groups_of_data)
    colors = ['DodgerBlue','Tomato','Lime', 'DimGray', 'DarkGray', 'PapayaWhip', 'RoyalBlue', 'PeachPuff' 'SlateGray', 'Silver', 'DodgerBlue', 'DeepSkyBlue', 'RoyalBlue']
    width = 0.9
    y_lim=False
    pylab.clf()
    pylab.rc('text', usetex=True)
    fig = pylab.figure(tight_layout=True)
    ax = pylab.axes()

    # Set color of each box in box and whisker plot
    for i,group in enumerate(groups_of_data):
        pos = [ i*(num_groups+1) + p for p in range(1,num_groups+1) ]
        bp = pylab.boxplot(group, positions = pos, widths = width, sym = '')
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
        pylab.setp(bp['whiskers'], color='Black', linestyle='solid')

        # Change the median line of each box in the box and whisker plot from red to black
        for m in range(num_groups):
            med = bp['medians'][m]
            medianX = []
            medianY = []
            for j in range(2):
                medianX.append(med.get_xdata()[j])
                medianY.append(med.get_ydata()[j])
                pylab.plot(medianX, medianY, 'k') # , linewidth=1.5

    # Plot overall median value for y-values ( y = median(ys))
    pylab.plot( [0,(num_groups+1)*len(groups_of_data)], [np.median(ys), np.median(ys)], 'k--', linewidth=0.5)

    # Set x-limits and y-limits
    pylab.xlim( 0, (num_groups+1)*len(groups_of_data) )
    if y_lim:
        pylab.ylim(y_lim)

    # Set tick positions, tick labels, x-labels, y-labels
    tickpos = []
    for i in range(len(groups_of_data)):
        tickpos.append( i*(num_groups+1) + (num_groups+1)/2.0 )
    ax.set_xticks(tickpos)
    ax.set_xticklabels(['Bottom','Second','Middle','Fourth','Top'], rotation=0, fontsize=15)
    pylab.xlabel(xlabel, fontsize = 30)
    pylab.ylabel(ylabel, fontsize = 30)
    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(20)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(20)

    # Remove right and top spines
    pylab.gca().spines['right'].set_color('none')
    pylab.gca().spines['top'].set_color('none')
    pylab.gca().xaxis.set_ticks_position('bottom')
    pylab.gca().yaxis.set_ticks_position('left')

    pylab.savefig('%s' % plotfile)
    pylab.clf()
    pylab.close()


def main():
    """Main body of script.
    """

    # Parse the configuration file:
    configfilename = sys.argv[1]
    if not os.path.isfile(configfilename):
        raise IOError("Failed to find configuration file %s" % configfilename)
    d = mapmuts.io.ParseInfile(open(configfilename))
    basedir = mapmuts.io.ParseStringValue(d, 'base_dir')
    
    # If True we don't overwrite existing output; if False we regenerate everything
    use_existing_output = True

    # Do we use sbatch to submit some of the jobs?
    use_sbatch = True

    # density for JPGs created from PDFs
    jpg_density = 150

    # Names of host species for which phylogenetic trees should already exist in the directory phylo_input/
    hosts = ['Allhosts', 'Human', 'Swine', 'Equine', 'Avian']

    # Prefix for sequence alignment files. This prefix and the suffix should make the entire file name for the sequence alignment files.
    # These alignments should already exist in the directory phylo_input/
    infileprefix = 'Aligned_NPs' 

    # List of amino-acid preference files, which should already exist in the directory dmstools_output/. 
    preferencesfile = [ \
    '%s/dmstools_output/mean_Aichi68_both_studies_prefs.txt' % basedir, \
    '%s/dmstools_output/PR8_mean_prefs.txt' % basedir, \
    '%s/dmstools_output/mean_NP_both_homologs_prefs.txt' % basedir]

    # File containing RMSD values for Aichi1968 vs PR1934 should already exist in 
    # the directory compare_prefs_output/Aichi1968_vs_PR1934/
    RMSD_file = '%s/compare_prefs_output/Aichi1968_vs_PR1934/RMSD_Aichi1968_vs_PR1934_RMSD_calcs.txt' % basedir

    # List of string identifiers for subsitution models based on different experimental amino-acid preferences.
    # The order of strings corresponds to the order of preference files in preferencesfile list.
    preferencesstring = ['Aichi1968_', 'PR1934_', 'Aichi1968_PR1934_']

    # Suffices to add to end of file names (for reading FASTA alignment files and phylogenetic trees).
    # Alignment file names MUST match this format. 
    suffices = ['_%s' % i for i in hosts]

    print "Beginning execution of script at %s.\n" % time.asctime()
    print "Using directory %s" % basedir
    if use_existing_output:
        print "Existing output will be used when possible. Note that you want to set use_existing_output to False if you want to regenerate output, such as after changing input data or analysis settings."
    else:
        print "All existing output will be deleted, and new output regenerated."

    # separator to break sections of output
    separator = '*******************************************************************' 

    #########################################
    #########################################
    # Load each nucleoprotein sequence alignment
    #########################################
    #########################################
    print "\n%s\nLoading a set of aligned NP sequences..." % separator
    alignmentfiles = {}  # keyed by model abbreviation, value is corresponding alignmentfile
    for suffix in suffices:
        alignmentfile = '%s/phylo_input/%s%s.fasta' % (basedir, infileprefix, suffix)
        alignmentfiles['GY94%s' % suffix] = alignmentfile
        if os.path.isfile(alignmentfile):
            print "The existing alignment file of %s will be used." % alignmentfile
        else:
            raise ValueError("Failed to find alignment file %s" % alignmentfile)

    #########################################
    #########################################
    # Load the codonPhyML trees that were built using the above nucleoprotein sequence alignments
    #########################################
    #########################################
    print "\n%s\nLoading a tree built by codonPhyML from the aligned sequences..." % separator

    # codonphyml_trees is keyed by model abbreviation, value is phylogenetic tree file
    codonphyml_trees = {}
    for suffix in suffices:
        codonphyml_trees['GY94%s' % suffix] = '%s/phylo_input/codonphyml%s_tree.newick' % (basedir, suffix)
    for (model, tree) in codonphyml_trees.iteritems():
        if os.path.isfile(tree):
            print "The existing tree file %s for model %s will be used." % (tree, model)
        else:
            raise ValueError("Failed to find tree file %s for model %s" % (tree, model))
    print "\n%s\n\n" % separator

    #########################################
    #########################################
    # Optimize the trees for the different substitution models
    #########################################
    #########################################
    print "\n%s\nOptimizing the tree branch lengths for different substitution models and computing the resulting likelihoods..." % separator

    fixationmodels = ['HalpernBruno']
    hyphy_results = {} # keyed by (treemodel, substitionmodel) to give hyphyoutfile
    sbatch_cpus = 4 # get four CPUs as script uses lots of memory, and having more CPUs typically means less other jobs on the node
    walltime =  9 * 24 # give up to nine days for optimization
    processes = []
    optimizedtreedir = "%s/phylo_output/codonmodel_optimized_trees" % basedir
    if not os.path.isdir(optimizedtreedir):
        os.mkdir(optimizedtreedir)

    # Prepare processes to run with experimentally determined substitution models (site-specific)
    commands = {
                'hyphypath':'HYPHYMP CPU=4',
                'mutationrates':'freeparameters', 
                'scalefactor':10000.0,
                'siteslist':preferencesfile[0],
                'keeptempfiles':'False',
                'outfileprefix':'None',
                'persitelikelihoods':'True',
                }
    experimentalmodels = []
    for fixationmodel in fixationmodels:
        commands['fixationmodel'] = fixationmodel
        for (treemodel, tree) in codonphyml_trees.iteritems():
            commands['treefile'] = tree
            commands['fastafile'] = alignmentfiles[treemodel]
            for (fitbeta, fitbetastring) in [('freeparameter', 'fitbeta')]:
                commands['fitbeta'] = fitbeta
                for (preffile, prefstring) in zip(preferencesfile, preferencesstring):
                    commands['aapreferences'] = preffile
                    for (randomizepreferences, randstring) in [('False', '')]:
                        commands['randomizepreferences'] = randomizepreferences
                        istring = '%s%s%s%s' % (prefstring, fitbetastring, fixationmodel, randstring)
                        subdir = '%s/Tree-%s_Model-%s' % (optimizedtreedir, treemodel, istring)
                        hyphyoutfile = '%s/optimizedtree_results.txt' % subdir
                        if os.path.isfile(hyphyoutfile) and use_existing_output:
                            print "The output for the HYPHY analysis of the %s tree with the %s fixation model already exists in %s" % (treemodel, istring, hyphyoutfile)
                        else:
                            print "Using HYPHY to analyze the %s tree with the %s fixation model to create output file %s" % (treemodel, istring, hyphyoutfile)
                            processes.append(multiprocessing.Process(target=RunScript,\
                            args=(subdir, 'phyloExpCM_ExpModelOptimizeHyphyTree', 'phyloExpCM_ExpModelOptimizeHyphyTree.py', list(commands.items()), use_sbatch, sbatch_cpus), kwargs={'walltime':walltime}))
                        experimentalmodels.append(istring)
                        hyphy_results[(treemodel, istring)] = hyphyoutfile

    # Prepare processes to run for non site-specific substitution models
    substitutionmodels = ['GY94_CF3x4_omega-global-gamma4_rates-gamma4']
    substitutionmodels = dict([(x, x) for x in substitutionmodels])
    commands = {
                'hyphypath':'HYPHYMP CPU=1',
                'hyphycmdfile':'hyphy_cmds.bf',
                'hyphyoutfile':'hyphy_output.txt',
                'hyphytreefile':'hyphy_tree.newick',
                'hyphydistancesfile':'None',
                'siteslist':preferencesfile[0],
                'persitelikelihoods':'sitelikelihoods.txt'
               }
    for (treemodel, tree) in codonphyml_trees.iteritems():
        commands['fastafile'] = alignmentfiles[treemodel]
        commands['treefile'] = tree
        for (substitutionmodel, substitutionmodelspecs) in substitutionmodels.iteritems():
            commands['model'] = substitutionmodelspecs
            subdir = '%s/Tree-%s_Model-%s' % (optimizedtreedir, treemodel, substitutionmodel)
            hyphyoutfile = '%s/hyphy_output.txt' % subdir
            if os.path.isfile(hyphyoutfile) and use_existing_output:
                print "The output for the HYPHY analysis of the %s tree with the %s substitution model already exists in %s" % (treemodel, substitutionmodel, hyphyoutfile)
            else:
                print "Using HYPHY to analyze the %s tree with the %s substitution model to create output file %s" % (treemodel, substitutionmodel, hyphyoutfile)
                processes.append(multiprocessing.Process(target=RunScript,\
                    args=(subdir, 'phyloExpCM_optimizeHyphyTree', 'phyloExpCM_optimizeHyphyTree.py', list(commands.items()), use_sbatch, 1), kwargs={'walltime':walltime}))
            hyphy_results[(treemodel, substitutionmodel)] = hyphyoutfile

    # Run phylogenetic fitting for each substitution model
    RunProcesses(processes, nmultiruns=len(processes))
    for ((treemodel, substitutionmodel), hyphyoutfile) in hyphy_results.iteritems():
        if not os.path.isfile(hyphyoutfile):
            raise ValueError("Failed to find expected HYPHY output for the %s tree with the %s substitution model, which should be in the file %s" % (treemodel, substitutionmodel, hyphyoutfile))
    for x in experimentalmodels:
        substitutionmodels[x] = x
    print "\nAll HYPHY analyses have been completed."
    print "\n%s\n\n" % separator

    #########################################
    #########################################
    # Make summary output of log likelihood and number of parameters
    #########################################
    #########################################
    print "\n%s\nCreating summaries of log likelihoods and parameter counts.\n" % separator

    empiricalparameters = [ # (unique re matching model name, number empirical parameters)
            (re.compile('^Aichi1968_fit'), 0),
            (re.compile('^PR1934_fit'), 0),
            (re.compile('^Aichi1968_PR1934_fit'), 0),
            (re.compile('GY94_CF3x4'), 9),
            ]
    # keyed by parameter strings in HYPHY output, value is printed parameter name
    parameters = {
                'beta':'$\\beta$' ,
                'R_AG':'$R_{A \\rightarrow G}$',
                'R_AT':'$R_{A \\rightarrow T}$',
                'R_CA':'$R_{C \\rightarrow A}$',
                'R_CG':'$R_{C \\rightarrow G}$',
                'global rate_alpha':'rate shape',
                'global omega':'$\omega$',
                'global omega_mean':'mean $\omega$',
                'global kappa':'$\kappa$',
                'global omega_alpha':'$\omega$ shape',
                }

    for treemodel in codonphyml_trees:
        fname = '%s/phylo_output/%s_summary.csv' % (basedir, treemodel)
        flatexname = '%s/phylo_output/%s_summary.tex' % (basedir, treemodel)
        print "Writing summary for tree %s to %s and %s..." % (treemodel, fname, flatexname)
        linelist = []
        latexlinelist = []
        for substitutionmodel in substitutionmodels:
            eparameters = None
            for (m, n) in empiricalparameters:
                if m.search(substitutionmodel):
                    if eparameters != None:
                        raise ValueError("%s matches two empiricalparameters" % substitutionmodel)
                    else:
                        eparameters = n
            if eparameters == None:
                raise ValueError("empiricalparameters not specified for %s" % substitutionmodel)
            resultfile = hyphy_results[(treemodel, substitutionmodel)]
            treeparameters = dict([(parameter, None) for parameter in parameters.iterkeys()])
            if 'fitbeta' not in substitutionmodel:
                del treeparameters['beta']
            phyloExpCM.hyphy.ExtractValues(resultfile, treeparameters, allowmissing=True)
            (loglikelihood, mlparameters) = phyloExpCM.hyphy.ExtractLikelihoodAndNParameters(resultfile)
            totparameters = mlparameters + eparameters
            aic = 2.0 * totparameters - 2.0 * loglikelihood
            textstring = "%s, %%g, %g, %d, %d, %d" % (substitutionmodel, loglikelihood, totparameters, mlparameters, eparameters)
            latexsubstitutionmodel = FormatModelName(substitutionmodel)

            print 'SUBSTITUTION MODEL IS %s latex is %s' % (substitutionmodel, latexsubstitutionmodel)

            treeparameters = ['%s = %.1f' % (parameters[parameter], value) for (parameter, value) in treeparameters.iteritems() if value]
            treeparameters.sort()
            treeparameters = ', '.join(treeparameters)
            if not treeparameters:
                treeparameters = 'none'
            latexstring = "%s & %%.1f & %.1f & %d (%d + %d) & %s" % (latexsubstitutionmodel, loglikelihood, totparameters, mlparameters, eparameters, treeparameters)
            linelist.append((aic, textstring))
            latexlinelist.append((aic, latexstring))
        linelist.sort()
        latexlinelist.sort()
        aics = [tup[0] for tup in linelist]
        daics = [aic - min(aics) for aic in aics]
        linelist = [line % daic for ((aic, line), daic) in zip(linelist, daics)]
        latexlinelist = [line % daic for ((aic, line), daic) in zip(latexlinelist, daics)]
        f = open(fname, 'w')
        flatex = open(flatexname, 'w')
        f.write('#Summary for tree %s.\n#\n#SUBSTITUTION_MODEL, dAIC, LOG_LIKELIHOOD, FREE_PARAMETERS, MAXIMUM_LIKELIHOOD_PARAMETERS, EMPIRICAL_PARAMETERS\n%s' % (treemodel, '\n'.join([line for line in linelist])))
        flatex.write('{\\scriptsize\n\\begin{tabular}{ccccc}\nmodel & $\Delta$AIC & \parbox[b]{0.53in}{\center log likelihood} & \parbox[b]{0.7in}{\center parameters (optimized + empirical)} & optimized parameters \\\\ \hline\n%s\n\end{tabular}}' % '\\\\ \n'.join([line for line in latexlinelist]))
        f.close()
        flatex.close()
    print "\nAll summaries created.\n%s\n\n" % separator
    
    #########################################
    #########################################
    # Compare site likelihoods for two experimental site-specific substitution models fitting the Allhosts phylogenetic tree
    #########################################
    #########################################
    print "\n%s\nComparing site likelihoods for two experimental site-specific substitution models.\n\n" % separator

    host = 'Allhosts'
    (model1name, model1file) = ('Aichi1968_PR1934', 'codonmodel_optimized_trees/Tree-GY94_%s_Model-Aichi1968_PR1934_fitbetaHalpernBruno/sitelikelihoods.txt' % host)
    (model2name, model2file) = ('Aichi1968', 'codonmodel_optimized_trees/Tree-GY94_%s_Model-Aichi1968_fitbetaHalpernBruno/sitelikelihoods.txt' % host)

    prefix = model1name + '_minus_' + model2name + '_' + host + '_'
    subdir = '%s/phylo_output' % basedir
    print "The two models being compared are %s (%s) versus %s (%s)" % (model1name, model1file, model2name, model2file)
    commands = [('sitelikelihoodfiles', '%s %s' % (model1file, model2file)),
                ('modelnames', '%s %s' % (model1name, model2name)),
                ('dsspfile', 'None'),
                ('dsspchain', 'None'),      
                ('outfileprefix', prefix)]
    RunScript(subdir, 'phyloExpCM_SiteLikelihoodComparison', 'phyloExpCM_SiteLikelihoodComparison.py', commands, False, 1)

    print "\nCompleted the site likelihood comparison.\n%s\n\n" % separator

    #########################################
    #########################################
    # Bin nucleoprotein sites by how much they improve the per-site likelihood in going from Aichi1968 model to Aichi1968_PR1934 model.
    # Within each bin of sites, display the corresponding RMSD_within or RMSD_corrected values for those sites.
    #########################################
    #########################################
    print "\n%s\nBinning nucleoprotein sites by improvement in per-site likelihood and correlating with RMSD_within or RMSD_corrected.\n\n" % separator

    host = 'Allhosts'

    # Read file containing delta per-site log likelihood values of Aichi1968_PR1934 minus Aichi1968
    subdir = '%s/phylo_output' % basedir
    likelihoodfile = '%s/Aichi1968_PR1934_minus_Aichi1968_%s_sitelikelihoods.txt' % (subdir, host)
    sitelikelihood = ReadSiteLikelihoodFile(likelihoodfile)
    prefix = re.sub('_[A-Za-z0-9]*_sitelikelihoods.txt', '', os.path.basename(likelihoodfile))

    # Read RMSD file with RMSD values for Aichi/1968-PR/1934 calculated by scripts/compare_preferences.py
    RMSD = {}
    rmsd_between, rmsd_within, rmsd_corrected = rmsdtools.ReadRMSDFile(RMSD_file) 
    RMSD['rmsd_between'] = rmsd_between
    RMSD['rmsd_within'] = rmsd_within
    RMSD['rmsd_corrected'] = rmsd_corrected

    # Plot binned RMSDwithin vs delta log likelihood for protein sites shared between two data sets
    xlabel = '$\Delta$ (log-likelihood) quintile'
    ylabel = '$RMSD_{within}$'
    plotfile = '%s/%s_%s_binnedlikelihood_RMSDwithin.pdf' % (subdir, prefix, host)
    sites = list( set(sitelikelihood.keys()) & set([int(site) for site in RMSD['rmsd_within']]) )
    xs, ys = zip( *[(sitelikelihood[site], RMSD['rmsd_within'][str(site)]) for site in sites] )
    BinValuesByPercentileBoxPlot(xs, ys, plotfile, xlabel, ylabel, 20)

    # Plot binned RMSDcorrected vs delta log likelihood for protein sites shared between two data sets
    xlabel = '$\Delta$ (log-likelihood) quintile'
    ylabel = '$RMSD_{corrected}$'
    plotfile = '%s/%s_%s_binnedlikelihood_RMSDcorrected.pdf' % (subdir, prefix, host)
    sites = list( set(sitelikelihood.keys()) & set([int(site) for site in RMSD['rmsd_corrected']]) )
    xs, ys = zip( *[(sitelikelihood[site], RMSD['rmsd_corrected'][str(site)]) for site in sites] )
    BinValuesByPercentileBoxPlot(xs, ys, plotfile, xlabel, ylabel, 20)

    print "\nCompleted execution of script at %s.\n" % time.asctime()


if __name__ == '__main__':
    main() # run the script
