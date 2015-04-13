'''This script plots amino-acid preference correlations between samples and makes heatmaps of correlations. This was used
to make Figure 2 in the paper.

The correlations to plot are hardcoded in the script; however, a configuration file is needed to specify some paths.

This script should be run with one option at the command line which specifies a configuration file, eg:

python make_correlation_plots.py configfile.txt

The configuration file specifies paths to various directories containing scripts and preference files.
Each line of the configuration file should have a pathname and a path, separated by a space, for the following pathnames:

heatmap_script_path: path to the correlation_heatmap.py script.
preferencefiles_dir: path to the preference files made by `dms_tools`. likely ends in '/dmstools_output/'
plot_output_dir: path to a directory to save the correlation plots and heatmaps.

example configuration file:

heatmap_script_path /home/user/NP_analysis/scripts/correlation_heatmap.py
preferencefiles_dir /home/user/NP_analysis/dmstools_output/
plot_output_dir /home/user/NP_analysis/correlation_plots/
'''

import os
import sys
import time
import string
import rmsdtools
import mapmuts.io
import numpy as np
import mapmuts.sequtils
import itertools
import sitegroups
import scipy


def RunScript(rundir, run_name, script_name, commands, use_sbatch, sbatch_cpus, walltime=None):
    """Runs a `python` script with a single argument being an infile composed of key/value pairs
    specified in the dictionary *commands* like so:

    python script_name infile

    *rundir* is the directory in which we run the job. Created if it does
    not exist.

    *run_name* is the name of the run, which should be a string without
    spaces. The input file has this prefix followed by ``_infile.txt``.

    *script_name* is the name of the script that we run. This should contain
    the full path to the script, which should be hardcoded once in the main function below.

    *commands* contains the commands written to the input file. Itis a list 
    of 2-tuples giving the key / value pairs.
    Both keys and values should be strings.

    *use_sbatch* is a Boolean switch specifying whether we use ``sbatch``
    to run the script. If *False*, the script is just run with the command
    line instruction. If *True*, then ``sbatch`` is used, and the command file
    has the prefix *run_name* followed by the suffix ``.sbatch``.

    *sbatch_cpus* is an option that is only meaningful if *use_sbatch* is 
    *True*. It gives the integer number of CPUs that are claimed via
    ``sbatch`` using the option ``sbatch -c``. 

    *waltime* is an option that is only meaningful if *use_sbatch* is
    *True*. If so, it should be an integer giving the number of hours 
    to allocate for the job. If *walltime* has its default value of 
    *None*, no wall time for the job is specified.

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
        sbatch_f.write('python %s %s' % (script_name, infile))
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
                if nslurmfails > 180: # error is squeue fails at least 180 consecutive times
                    raise ValueError("squeue is continually failing, which means that slurm is not working on your system. Note that although this script has crashed, many of the jobs submitted via slurm may still be running. You'll want to monitor (squeue) or kill them (scancel) -- unfortunately you can't do that until slurm starts working again.")
                continue # we got an error while trying to run squeue
            nslurmfails = 0
            lines = open(jobstatusfile).readlines()
            if len(lines) < 2:
                break # no longer in slurm queue
        errors = open(joberrorsfile).read().strip()
    else:
        errors = os.system('python %s %s' % (script_name, infile))
    os.chdir(currdir)
    if errors:
        print "ERROR running %s for %s in directory %s." % (script_name, run_name, rundir)
        return True
    else:
        print "Successfully completed running %s for %s in directory %s." % (script_name, run_name, rundir)
        return False



def main():
    
    # Parse config file:
    configfilename = sys.argv[1]
    if not os.path.isfile(configfilename):
        raise IOError("Failed to find configuration file %s" % configfilename)
    d = mapmuts.io.ParseInfile(open(configfilename))
    heatmap_script_path = mapmuts.io.ParseStringValue(d, 'heatmap_script_path')
    print "Using heatmap script %s" % heatmap_script_path
    preferencefiles_dir = mapmuts.io.ParseStringValue(d, 'preferencefiles_dir')
    print "Using preferencefiles_dir %s" % preferencefiles_dir
    plot_output_dir = mapmuts.io.ParseStringValue(d, 'plot_output_dir')
    print "Using plot_output_dir %s" % plot_output_dir

    sites = range(2,499)
    aas = mapmuts.sequtils.AminoAcids(includestop=False)

    all_sites = sites
    variable_sites = sitegroups.VariableSites()
    rna_sites = sitegroups.RNABindingGrooveSites()
    conserved_contact_sites, all_contact_sites = sitegroups.SitesContactingVariableSites()

    site_groups = [all_sites, rna_sites, variable_sites, conserved_contact_sites]
    site_group_labels = ["All sites", "RNA-binding", "Variable", "Contacting variable"]

    if not os.path.isdir(plot_output_dir):
        os.mkdir(plot_output_dir)
    os.chdir(plot_output_dir)

    # Specify the names of the preference files for individual replicates to be used for correlation heatmaps
    # All these files should be in the mapmuts_output/preferencefiles directory after running the `run_mapmuts.py` script.
    PR8_preference_files = ['PR8_replicate_1_prefs.txt',
                            'PR8_replicate_2_prefs.txt',
                            'PR8_replicate_3_prefs.txt']
    Aichi68A_preference_files = [   'Aichi68A_replicate_1_prefs.txt',
                                    'Aichi68A_replicate_2_prefs.txt',
                                    'Aichi68A_replicate_3_prefs.txt',
                                    'Aichi68A_replicate_4_prefs.txt'
                                    ]
    Aichi68B_preference_files = [   'Aichi68B_replicate_1_prefs.txt',
                                    'Aichi68B_replicate_2_prefs.txt',
                                    'Aichi68B_replicate_3_prefs.txt',
                                    'Aichi68B_replicate_4_prefs.txt'
                                    ]
    Aichi68C_preference_files = [   'Aichi68C_replicate_1_prefs.txt',
                                    'Aichi68C_replicate_2_prefs.txt']
    WSN_HA_files = [ 'WSN_HA_rep1_prefs.txt', 'WSN_HA_rep2_prefs.txt', 'WSN_HA_rep3_prefs.txt' ]

    # Assemble the list of preference files to use in heatmaps 
    all_pref_files = PR8_preference_files + Aichi68A_preference_files + Aichi68B_preference_files + Aichi68C_preference_files + WSN_HA_files
    for i,f in enumerate(all_pref_files):
        all_pref_files[i] = preferencefiles_dir + '/' + f
    labels = [  'PR8-1', 'PR8-2', 'PR8-3', 'Aichi68-A1', 'Aichi68-A2', 'Aichi68-A3', 'Aichi68-A4', 
                'Aichi68-B1', 'Aichi68-B2', 'Aichi68-B3', 'Aichi68-B4',  
                'Aichi68-C1', 'Aichi68-C2', 'WSN-HA-1', 'WSN-HA-2', 'WSN-HA-3']
    assert len(labels) == len(all_pref_files)

    # Make heatmaps of correlation coefficients for all sites and for subgroups of sites.
    for site_group, site_group_label in zip(site_groups, site_group_labels):
        command_d = {   'includestop' : 'False',
                        'outputprefix' : '%s_correlation_heatmap' % site_group_label.replace(" ",""),
                        'sites_type' : 'list',
                        'sites' : str(site_group).rstrip(']').lstrip('['),
                        'preference_files' : ' '.join(all_pref_files),
                        'sample_names' : ' '.join(labels)
                         }
        RunScript(plot_output_dir, '%s_generate_heatmap' % site_group_label.replace(" ", ""), heatmap_script_path, 
                    list(command_d.items()), False, 1)
    
    # Specify location of mean preferences for correlation plots
    mean_HA = "%s/WSN_HA_mean_prefs.txt" % preferencefiles_dir
    mean_PR8 = "%s/PR8_mean_prefs.txt" % preferencefiles_dir
    mean_AichiOld = "%s/mean_Aichi68_previous_study_prefs.txt" % preferencefiles_dir
    mean_AichiNew = "%s/Aichi68C_mean_prefs.txt" % preferencefiles_dir
    mean_Aichi = "%s/mean_Aichi68_both_studies_prefs.txt" % preferencefiles_dir

    # Make correlation plots for all sites comparing mean preferences from various experiments:
    rmsdtools.PlotAAPrefCorrelations(mean_PR8, mean_Aichi, sites, "PR/1934 preference", "Aichi/1968 preference", "pr8_aichi_correlation", read_from_file=True)
    rmsdtools.PlotAAPrefCorrelations(mean_AichiOld, mean_AichiNew, sites, "Aichi/1968 preference\n(previous study)", "Aichi/1968 preference\n(current study)", "aichi_aichi_correlation", read_from_file=True)
    rmsdtools.PlotAAPrefCorrelations(mean_PR8, mean_HA, sites, "PR/1934 preference", "HA preference", "pr8_ha_correlation", read_from_file=True)


main()
