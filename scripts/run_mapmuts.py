"""This script runs a suite of ``mapmuts`` scripts to align deep sequencing reads and parse mutation counts for many samples.

This script should be run with one option at the command line which specifies the path of a configuration file, eg:

python run_mapmuts.py configfile.txt

The configuration file specifies paths to various directories containing FASTQ files, previously published amino-acid preferences, etc.
Each line of the configuration file should have a pathname and a path, separated by a space, for the following pathnames:

adapter_dir: contains FASTA files with adapter sequences
refseq_dir: contains reference NP sequences for PR/1934, Aichi/1968
mapmuts_output_dir: output directory for mapmuts analysis
FASTQ_directory: contains FASTQ files in the appropriate directory structure (/PR8_replicate_X/DNA, /Aichi68C_replicate_Y/virus, etc)

example configuration file:

adapter_dir /home/user/NP_analysis/adapter_seqs/
refseq_dir /home/user/NP_analysis/ref_seqs/
mapmuts_output_dir /home/user/NP_analysis/mapmuts_output/
FASTQ_directory /home/user/NP_analysis/FASTQ_files/

---

This is a master Python script that runs the ``mapmuts`` scripts to analyze deep sequencing data. 
Can run multiple scripts at the same time, using either direction submission
of the jobs on the current CPU or submission to a queue via ``sbatch``.

The script is currently set up to submit jobs using sbatch. To use the current
CPU for all jobs, set `use_sbatch` in the main() function to False.

The locations and names of various files are hard-coded into the script.
Because this is a custom script for running jobs on the FHCRC computing core, 
the code is not fully documented elsewhere, and you will have to look at the 
script source code to understand exactly how it works. 
Note however that the ``mapmuts`` scripts run by this code are fully documented 
with the ``mapmuts`` package, and do not requiring understanding the source code.
This script is just running the ``mapmuts`` scripts.
This script was written by Jesse Bloom and Mike Doud.
"""


import os
import sys
import string
import time
import multiprocessing
import mapmuts.io


def RunScript(rundir, run_name, script_name, commands, use_sbatch, sbatch_cpus, walltime=None):
    """Runs a ``mapmuts`` script.

    *rundir* is the directory in which we run the job. Created if it does
    not exist.

    *run_name* is the name of the run, which should be a string without
    spaces. The input file has this prefix followed by ``_infile.txt``.

    *script_name* is the name of the script that we run.

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
                if nslurmfails > 180: # error is squeue fails at least 180 consecutive times
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
    have not yet been started.

    *nmultiruns* is an integer >= 1 indicating the number of simultaneous
    processes to run.

    Runs the processes in *processes*, making sure to never have more than
    *nmultiruns* running at a time. If any of the processes fail (return
    an exitcode with a boolean value other than *False*), an exception
    is raised immediately. Otherwise, this function finishes when all
    processes have completed.
    """
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
        

def main():
    
    # parse the configuration file:
    configfilename = sys.argv[1]
    if not os.path.isfile(configfilename):
        raise IOError("Failed to find configuration file %s" % configfilename)
    d = mapmuts.io.ParseInfile(open(configfilename))
    adapter_dir = mapmuts.io.ParseStringValue(d, 'adapter_dir')
    print "Using adapter sequence directory %s" % adapter_dir
    refseq_dir = mapmuts.io.ParseStringValue(d, 'refseq_dir')
    print "Using reference sequence directory %s" % refseq_dir
    mapmuts_output_dir = mapmuts.io.ParseStringValue(d, 'mapmuts_output_dir')
    print "Using output directory for mapmuts analysis %s" % mapmuts_output_dir
    FASTQ_files_directory = mapmuts.io.ParseStringValue(d, 'FASTQ_directory')

    plotout_dir = mapmuts_output_dir + '/mapmuts_plots/' # for mapmuts plots that don't fall under directories for each replicate
    
    if not os.path.isdir(mapmuts_output_dir):
        os.mkdir(mapmuts_output_dir)
    if not os.path.isdir(plotout_dir):
        os.mkdir(plotout_dir)

    # Do we use sbatch when appropriate? Some simple operations will still
    # be run without sbatch even when True
    use_sbatch = True
    # Maximum number of CPUs to try to use at once. If not using sbatch, 
    # don't make this bigger than the number of available cores.
    max_cpus = 100 
    # Density for JPGs converted from PDFs
    jpg_density = 150

    # Specify the replicates and amplicons.
    replicates = ['PR8_replicate_1', 'PR8_replicate_2', 'PR8_replicate_3', 'Aichi68C_replicate_1', 'Aichi68C_replicate_2']
    amplicons = ['DNA', 'mutDNA', 'virus', 'mutvirus']

    # here you can toggle which steps to run or skip. However, this script is not smart enough to 
    # know if the output required from each step exists or is correct, so use with caution! 
    # These sections should generally be run in order.
    run_makealignments = True
    run_alignmentsummaryplot = True
    run_parsecounts = True
    run_parsesummaryplots = True
    run_countparsedmuts = True

    #########################################
    #########################################
    # run mapmuts_makealignments.py. 
    # Each sample is run its own subdirectory.
    #########################################
    #########################################
    if run_makealignments:
        print "entering makealignments section of masterscript"

        processes = []
        convert_to_jpgs = []

        # dictionary to convert replicate/amplicon to appropriate barcoded adapter sequence for trimming
        R1_trims_d = {  'PR8_replicate_1_DNA' : 'R1_adapter_AR001.fa',
                        'PR8_replicate_1_mutDNA' : 'R1_adapter_AR008.fa',
                        'PR8_replicate_1_virus' : 'R1_adapter_AR010.fa',
                        'PR8_replicate_1_mutvirus' : 'R1_adapter_AR011.fa',
                        'PR8_replicate_2_DNA' : 'R1_adapter_AR003.fa', 
                        'PR8_replicate_2_mutDNA' : 'R1_adapter_AR009.fa',
                        'PR8_replicate_2_virus' : 'R1_adapter_AR025.fa',
                        'PR8_replicate_2_mutvirus' : 'R1_adapter_AR022.fa',
                        'PR8_replicate_3_DNA' : 'R1_adapter_AR002.fa',
                        'PR8_replicate_3_mutDNA' : 'R1_adapter_AR004.fa',
                        'PR8_replicate_3_virus' : 'R1_adapter_AR007.fa',
                        'PR8_replicate_3_mutvirus' : 'R1_adapter_AR016.fa', 
                        'Aichi68C_replicate_1_DNA' : 'R1_adapter_AR003.fa',
                        'Aichi68C_replicate_1_mutDNA' : 'R1_adapter_AR009.fa',
                        'Aichi68C_replicate_1_virus' : 'R1_adapter_AR025.fa',
                        'Aichi68C_replicate_1_mutvirus' : 'R1_adapter_AR022.fa',
                        'Aichi68C_replicate_2_DNA' : 'R1_adapter_AR003.fa',
                        'Aichi68C_replicate_2_mutDNA' : 'R1_adapter_AR009.fa',
                        'Aichi68C_replicate_2_virus' : 'R1_adapter_AR025.fa',
                        'Aichi68C_replicate_2_mutvirus' : 'R1_adapter_AR022.fa'
                        } 

        # Commands for makealignments which are the same for all samples:
        command_d = {'gzipped':'True',
                     'applyfilter':'True',
                     'minq':'25',
                     'a2file':'%s/R2_adapterUniversal_RC.fa' % adapter_dir,
                     'maxn':'5',
                     'minoverlap':'100',
                     'maxrm':'1',
                     'maxa1m':'1',
                     'maxa2m':'1',
                     'maxgenem':'10',
                     'upcase':'test',
                     'write_unaligned':'True',
                    }
          
        # set all makealignments commands that vary across replicates and amplicons:    
        for replicate in replicates:
            if not os.path.isdir(mapmuts_output_dir + '/' + replicate):
                os.mkdir(mapmuts_output_dir + '/' + replicate)

            if 'PR8' in replicate:
                command_d['generange'] = '46 1539'
                command_d['fullgenefile'] = '%s/PR8_NP.fa' % refseq_dir
            elif 'Aichi68' in replicate:
                command_d['generange'] = '62 1555'
                command_d['fullgenefile'] = '%s/Aichi68-NP_amplicon.fa' % refseq_dir
            else:
                raise ValueError("Invalid replicate name. Was expecting all replicates to be PR8 or Aichi68.")
 
            for amplicon in amplicons:
                if not os.path.isdir("%s/%s/%s" % (mapmuts_output_dir, replicate, amplicon)):
                    os.mkdir("%s/%s/%s" % (mapmuts_output_dir, replicate, amplicon))
                
                # specify the a1file for this sample using R1_trims_d.  Convert sample name to barcode-specific adapter sequence file.
                command_d['a1file'] = '%s/%s' % (adapter_dir, R1_trims_d[ replicate + '_' + amplicon ])
                
                command_d['r1files'] = '%s/%s/%s/*R1*.gz' % (FASTQ_files_directory, replicate, amplicon)
                command_d['r2files'] = '%s/%s/%s/*R2*.gz' % (FASTQ_files_directory, replicate, amplicon)
                          
                subdir = '%s/%s/%s' % (mapmuts_output_dir, replicate, amplicon)
                
                # outfileprefixes are specified as "replicate_amplicon"
                command_d['outfileprefix'] = replicate+'_'+amplicon
                # samplenames are specified by "replicate, amplicon"
                command_d['samplename'] = replicate+', '+amplicon
                
                # Add this makealignments.py call to the list of processes
                processes.append(multiprocessing.Process(target=RunScript,\
                    args=(subdir, "makealignments", 'mapmuts_makealignments.py', list(command_d.items()), use_sbatch, 1))) # use *list* to make a copy of *commands* sinces lists are mutable on subsequent runs
                
                convert_to_jpgs.append('%s/%s_alignmentstatistics.pdf' % (subdir, command_d['outfileprefix']))
                convert_to_jpgs.append('%s/%s_insertlengths.pdf' % (subdir, command_d['outfileprefix']))
                
        # run all the makealignments processes
        RunProcesses(processes, nmultiruns=max_cpus)
        
        # convert all the pdfs
        for f in convert_to_jpgs:
            assert os.path.isfile(f), "Cannot find file %s" % f
            os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))

        # end makealignments block of masterscript

    else:
        print "Skipping makealignments section of masterscript"



    #########################################
    #########################################
    # run mapmuts_alignmentsummaryplot.py to make one alignment summary plot
    #########################################
    #########################################
    if run_alignmentsummaryplot:
        print "entering alignmentsummaryplot section of masterscript"

        convert_to_jpgs = []
        commands = []
        for replicate in replicates:
            for amplicon in amplicons:
                commands.append(("%s/%s/%s/%s_%s_alignmentstatistics.txt" % (mapmuts_output_dir, replicate, amplicon, replicate, amplicon), "%s %s" % (replicate, amplicon)))
            
        plotfile = "%s/alignmentsummaryplot.pdf" % (plotout_dir)
        convert_to_jpgs.append(plotfile)
        commands.append(('plotfile', plotfile))
        RunScript(plotout_dir, 'alignmentsummaryplot', 'mapmuts_alignmentsummaryplot.py', commands, False, 1)
        
        for f in convert_to_jpgs:
            assert os.path.isfile(f), "Cannot find file %s" % f
            os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))
    else:
        print "skipping alignmentsummaryplot section of masterscript"


    #########################################
    #########################################
    # run mapmuts_parsecounts.py to parse counts from alignments.
    #########################################
    #########################################
    if run_parsecounts:
        print "entering parsecounts section of masterscript"
        processes = []
        convert_to_jpgs = []
        command_d = {'upcase':'test',
                     'r1exclude':'1 2 3 4 5 6 7 8 9 10 11 12 13 14 15',
                     'r2exclude':'1 2 3 4 5 6 7 8 9 10 11 12 13 14 15',
                    }
        for replicate in replicates:
            if 'PR8' in replicate:
                command_d['generange'] = '46 1539'
                command_d['fullgenefile'] = '%s/PR8_NP.fa' % refseq_dir
            elif 'Aichi68' in replicate:
                command_d['generange'] = '62 1555'
                command_d['fullgenefile'] = '%s/Aichi68-NP_amplicon.fa' % refseq_dir
            else:
                raise ValueError("Invalid type of replicate - strain not specified")
            for amplicon in amplicons:
                subdir = '%s/%s' % (replicate, amplicon)
                command_d['outfileprefix'] = replicate+'_'+amplicon
                command_d['samplename'] = replicate+', '+amplicon
                command_d['alignmentfile'] = '%s/%s/%s_alignments.txt.gz' % (mapmuts_output_dir, subdir, command_d['outfileprefix'])
                processes.append(multiprocessing.Process(target=RunScript,args=(mapmuts_output_dir+'/'+subdir, 'parsecounts', 'mapmuts_parsecounts.py', list(command_d.items()), use_sbatch, 1)))
                convert_to_jpgs.append('%s/%s/%s_codondepth.pdf' % (mapmuts_output_dir, subdir, command_d['outfileprefix']))    
        RunProcesses(processes, nmultiruns=max_cpus)
        for f in convert_to_jpgs:
            assert os.path.isfile(f), "Cannot find file %s" % f
            os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))
    else:
        print "skipping parsecounts section of masterscript"



    #########################################
    #########################################
    # run mapmuts_parsesummaryplots.py for each replicate
    #########################################
    #########################################
    if run_parsesummaryplots:
        print "entering parsesummaryplots section of masterscript"
        convert_to_jpgs = []
        for replicate in replicates:
            commands = []
            for amplicon in amplicons:
                commands.append(("%s/%s/%s/%s_%s" % (mapmuts_output_dir, replicate, amplicon, replicate, amplicon), "%s %s" % (replicate, amplicon)))
            plotfileprefix = "%s/%s/parsesummary" % (mapmuts_output_dir, replicate)
            commands.append(('plotfileprefix', plotfileprefix))
            commands.append(('writefracs', 'True'))
            commands.append(('textwritefracs', plotfileprefix))
            commands.append(('pairedcodonplot', 'True'))
            convert_to_jpgs.append("%s_codon_types_and_nmuts.pdf" % plotfileprefix)
            RunScript("%s/%s/" % (mapmuts_output_dir, replicate), 'parsesummaryplots', 'mapmuts_parsesummaryplots.py', commands, False, 1)
        for f in convert_to_jpgs:
            assert os.path.isfile(f), "Cannot find file %s" % f
            os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))
    else:
        print "skipping parsesummaryplots section of masterscript"



    #########################################
    #########################################
    # run mapmuts_countparsedmuts.py for each replicate, and for pooled samples across replicates within both strains.
    #########################################
    #########################################
    if run_countparsedmuts:
        print "entering countparsedmuts section of masterscript"

        convert_to_jpgs = []

        # First run individually for each replicate.
        for replicate in replicates:
            commands = [('plotfileprefix', 'countparsedmuts'), ('maxn', '50'), ('legendloc', 'right'), ('writecounts', 'False'), ('sites','2 498')]
            for amplicon in amplicons:
                files = ['%s/%s/%s/%s_%s_codoncounts.txt' % (mapmuts_output_dir, replicate, amplicon, replicate, amplicon)]
                commands.append((amplicon, ' '.join(files)))
            RunScript(mapmuts_output_dir + '/' + replicate, 'countparsedmuts', 'mapmuts_countparsedmuts.py', commands, False, 1)
            convert_to_jpgs.append('%s/%s/countparsedmuts_multi-nt-codonmutcounts.pdf' % (mapmuts_output_dir, replicate))

        # Now pool data across replicates for each strain
        # first, for PR8:
        PR8_replicates = ['PR8_replicate_1', 'PR8_replicate_2', 'PR8_replicate_3']
        commands = [('plotfileprefix', 'countparsedmuts_total_PR8'), ('maxn', '50'), ('legendloc', 'right'), ('writecounts', 'True'), ('sites', '2 498')]
        for amplicon in amplicons:
            files = ['%s/%s/%s/%s_%s_codoncounts.txt' % (mapmuts_output_dir, replicate, amplicon, replicate, amplicon) for replicate in PR8_replicates]
            commands.append((amplicon, ' '.join(files)))
        RunScript(plotout_dir, 'countparsedmuts_total_PR8', 'mapmuts_countparsedmuts.py', commands, False, 1)
        convert_to_jpgs.append('%s/countparsedmuts_total_PR8_multi-nt-codonmutcounts.pdf' % plotout_dir)

        # second, for Aichi68 current study:
        Aichi68C_replicates = ['Aichi68C_replicate_1', 'Aichi68C_replicate_2']
        commands = [('plotfileprefix', 'countparsedmuts_total_Aichi68C'), ('maxn', '50'), ('legendloc', 'right'), ('writecounts', 'True'), ('sites', '2 498')]
        for amplicon in amplicons:
            files = ['%s/%s/%s/%s_%s_codoncounts.txt' % (mapmuts_output_dir, replicate, amplicon, replicate, amplicon) for replicate in Aichi68C_replicates]
            commands.append((amplicon, ' '.join(files)))
        RunScript(plotout_dir, 'countparsedmuts_total_Aichi68C', 'mapmuts_countparsedmuts.py', commands, False, 1)
        convert_to_jpgs.append('%s/countparsedmuts_total_Aichi68C_multi-nt-codonmutcounts.pdf' % plotout_dir)

    else:
        print "skipping countparsedmuts section of masterscript"



if __name__ == '__main__':
    main() # run the script
