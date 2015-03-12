======================
Compare-NP-Preferences
======================

Summary
-------

These scripts perform a comparative analysis of site-specific amino-acid preferences inferred from deep mutational scanning data for two homologs of influenza nucleoprotein (NP). First, `mapmuts`_, a package for analyzing deep-sequencing data to identify mutation frequencies in mutant libraries of a gene, is used to process Illumina sequencing data (SRA Accession #SRP056028) for biological replicate mutant libraries of PR/1934(H1N1) and Aichi/1968(H3N2) NP. Next, `dms_tools`_, a package for analyzing deep mutational scanning data, is used to infer site-specific amino-acid preferences for each replicate mutant library. To quantify the extent that amino-acid preferences differ between NP homologs, ``compare_preferences.py`` calculates distance-based statistics that account for experimental noise by describing how different the amino-acid preferences are at each site *between* homologs relative to the experimental noise in amino-acid preferences across replicate experiments *within* each homolog. Finally, `phyloExpCM`_ is used to compare the performance of various codon substituion models based on the amino-acid preferences measured in each homolog to describe sequence evolution within human, swine, avian, and equine influenza NP phylogenetic trees. These scripts were written by Mike Doud and Orr Ashenberg and extensively use `mapmuts`_, `dms_tools`_, and `phyloExpCM`_, which were written by Jesse Bloom.

Running the analysis scripts
----------------------------

The entire analysis can be completed by sequentially running five master scripts: ``run_mapmuts.py``, ``run_dmstools.py``, ``make_correlation_plots.py``, ``compare_preferences.py``, and ``run_phyloanalysis.py``. Each of these master scripts are run at the command line with one argument, the name of a configuration file, which specifies paths to various directories that are used by that script (eg. ``python run_mapmuts.py mapmuts_config_file.txt``). Within each configuration file, each line must contain a directory name and the corresponding path, separated by a space (eg. ``FASTQ_directory /home/user/NP_homologs/FASTQ_files`` )

Running `mapmuts`_
------------------

The configuration file for ``run_mapmuts.py`` must contain entries for four paths:

  * ``adapter_dir`` specifies the directory containing adapter sequences used for trimming during sequence alignment (``R1_adapter_AR0XX.fa``).
  * ``refseq_dir`` specifies the directory containing the reference sequences used for alignment (``PR8_NP.fa`` and ``Aichi68-NP_amplicon.fa``).
  * ``FASTQ_directory`` specifies the directory containing subdirectories storing the FASTQ files. These subdirectories are structured as: ``FASTQ_directory/(replicate)/(amplicon)/sequences.fastq``, with possible replicates including ``PR8_replicate_1``, ``PR8_replicate_2``, ``PR8_replicate_3``, ``Aichi68C_replicate_1``, and ``Aichi68C_replicate_2``, and possible amplicons including ``DNA``, ``virus``, ``mutDNA``, and ``mutvirus``.
  *  ``mapmuts_output_dir`` specifies a directory in which the alignments and mutation parsing done by `mapmuts`_ will be saved.

Running `dms_tools`_
--------------------

The configuration file for ``run_dmstools.py`` ...

Running ``make_correlation_plots.py``
-------------------------------------

The configuration file for ``make_correlation_plots.py`` ...

Running ``compare_preferences.py``
----------------------------------

The configuration file for ``compare_preferences.py`` ...

Running ``run_phyloanalysis.py``
--------------------------------

The configuration file for ``run_phyloanalysis.py`` ...





.. _`mapmuts`: https://github.com/jbloom/mapmuts
.. _`dms_tools`: https://github.com/jbloom/dms_tools
.. _`Python`: http://www.python.org/
.. _`phyloExpCM`: https://github.com/jbloom/phyloExpCM
.. _`previously published`: http://dx.doi.org/10.1093/molbev/msu173
.. _`WSN-HA`: https://github.com/jbloom/mapmuts/tree/master/examples/WSN_HA_2014Analysis
