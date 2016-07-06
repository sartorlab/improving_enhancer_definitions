# Overview
This repository allows one to repeat much of the analysis that I did in my rotation in the Sartor lab in May/June 2016.

It includes the following:

- Create new enhancers lists using ENCODE ChromHMM/DNase-seq data
- Process [ENCODE ChIA-PET](https://www.encodeproject.org/search/?type=Experiment&assay_term_name=ChIA-PET) data using [Mango](https://github.com/dphansti/mango) / download other processed ChIA-PET data from other sources
- Use processed ChIA-PET data to link enhancers to genes
- Evaluate how well different enhancer lists capture ENCODE ChIP-seq peaks
- Run chipenrich using different enhancer lists and compare the resulting p-values

I've set it up so that everything can be run using a series of ```make``` commands.  The make commands are running the commands listed in commands files in the respective directories.  I've also added --help messages to most or all of the R scripts that are used (in the ```scripts``` directory), so hopefully by looking at the commands files and the options for the R scripts it will be clear how any desired modifications can be made. I've tried to scatter comments throughout the commands files and R scripts as well.

# Relevant file formats
The instructions below mention a number of file formats that I used.  They are described here. There also files named ```commands``` or ```*.commands``` throughout the directories; these are the files called by the Makefiles, listing the commands to be run.

## Pairs file

A pairs file ends in *.pairs

## Interactions file

An interactions file ends in *.interactions and lists interacting genomic regions based on e.g. a ChIA-PET experiment.  Each line in the file symbolizes an interaction between two genomic regions.  The first three columns are the chromosome, start, and end of the first of the two regions, and the second three columns are the chromosome, start, and end of the second of the two interacting regions.  This can be used to associate genes with enhancers.

My scripts assume that there are no headers in these files.


# Instructions
## Download public data
You'll need to start by downloading some public datasets (ENCODE ChIP-seq, ChromHMM tracks, and master DNase).  You can do this from this directory by running:

	make chipseq # To fetch ENCODE chip-seq data
	make dnase # To fetch ENCODE master DNase data
	make chromhmm # To fetch the ~9 cell types with ChromHMM results on ENCODE

## Generate enhancer lists
Once these have been downloaded, you can generate enhancer lists:

	cd new_enhancer_lists
	make enhancer_lists
	cd ..

You'll likely need to change the path at the top of the ```commands``` file in order for this to work, as I've set it for my own home path.  This will be the case in every commands file in other directories as well.

As it is now, this will generate four different enhancer lists:

- One list based on the ChromHMM enhancers alone
- One list based on the ChromHMM enhancers alone, extending enhancers less than 2kb in length to 2kb
- One list based on the ChromHMM enhancers unioned with the DNase regions (requiring that a DNase region shows up in at least 2 samples)
- One list based on the ChromHMM enhancers unioned with the DNase regions (requiring that a DNase region shows up in at least 2 samples), and extending enhancers less than 2kb in length to 2kb

Enhancers overlapping with the gene list (gene list is in current_definitions/genes.txt) are kicked out of the enhancer list.

## Gather ChIA-PET interaction data
Next you need to make sure that you have ChIA-PET interaction lists (simple text files listing the significantly interacting regions based on ChIA-PET data).  I've included all of the ones I've used in the ```interaction_lists``` directory (all files ending with \*interactions), so if you don't want to change anything then you can skip this step.  If you'd like to add new interactions files, just add them to this directory (with the file name following the format 'experiment_name.interactions') and they'll be included in the downstream processing (also make sure to add the experiment(s) to the ```experiment_info.txt``` file in the ```new_locus_lists``` directory -- explained later).

The files in the ```interaction_lists``` directory are from three sources:
- The files starting with "EN" (e.g., ```ENCSR000CAD.interactions```) are from ENCODE ChIA-PET data and have been processed with Mango (v. 1.1.7).  To recreate the processing pipelines for these files, go into the ```mango_pipelines``` directory and run:

		make prepare_workspace
		make ENCODE_pipelines

then you can run the resulting ```*.pipeline``` files individually with bash.
- The files ```ruan_2015*.interactions``` are from [this](http://www.ncbi.nlm.nih.gov/pubmed/26686651) publication.  I've simply downloaded their supplementary data; you can download the data again using the script ```interaction_lists/ruan_2015_data.sh```.  This data was __not__ processed using Mango, but rather with another method that (by the look of it) isn't as statistically stringent as Mango would be.
- The files ```naive_hesc.interactions``` and ```primed_hesc.interactions``` are from [this](doi:10.1016/j.stem.2015.11.007) publication.  They ran an earlier version of Mango on their data and the results are included in their supplementary files (as Excel files).  I downloaded these, converted them to text files, and included them in this repository.  It should be noted that although these were processed by the authors using Mango, it was an earlier version that I believe wasn't yet optimized for the protocol that they used; therefore, processing them with a newer version of Mango will give different results.

The number of interactions per experiment are in the spreadsheet [here](https://docs.google.com/a/umich.edu/spreadsheets/d/13h4WubwexHCBToSYIgumvLkR6S0gIuvcyubQYzyZSm0/edit?usp=sharing)

## Link enhancers to genes

Once the ChIA-PET data is all there, and you've generated the enhancer lists, you can link enhancers to genes.  To do this, ```cd``` into the ```new_locus_lists``` directory and run:

	make pipelines

This will create a number of *.pipeline files which you can simply run using bash.  When these are finished there will be files ending with "*.all.E.pairs" or "*.all.P2P.pairs" (described below).

This part of the pipeline links enhancers to genes in two different ways.  One is the "point-to-point" method (P2P).  This just means that if there are any enhancer(s) at one end of the interaction, and gene(s) at the other end, then these enhancers are assigned to these genes.  The second is the "encompassing" method (E).  This is something that I tried based on the CTCF-looping model of genome organization; if CTCF molecules interact to form loops, and genes and enhancers in these loops interact, then in CTCF/cohesin ChIA-PET data we may wish to take not only the _ends_ of the interactions but everything in between.  At the moment I have it set up so that the within a CTCF loop, if there is only one gene then this gene will be linked to all the enhancers in the loop.  If there is more than one gene present, then no assignments are made based on that loop (this is the --max_genes_per_region argument for the link_genes_and_enhancers.R script, as indicated by the --help menu for that script).

In order to use the 'E' method on just the CTCF/cohesin experiments, the scripts obviously need to know which proteins were pulled down in which experiments.  The experiment_info.txt file in this directory lists this information.  The first column lists the experiments (this should correspond to the experiment's interactions file name (in the interaction_lists directory), with the ".interactions" suffix removed); the second column lists the protein that was pulled down in the experiment.  Naturally this assumes that it's a ChIA-PET experiment.  If it's a HiC experiment or something like that and you'd like to use the experiment to create the locus definitions, just treat it like a ChIA-PET experiment and put "CTCF" in the second column of the experiment_info.txt file.

The output files at this step (using bash to run the *.pipeline files) are *.pairs files (one for each enhancer list and experiment, corresponding to the "E" and the "P2P" methods.  If it's not a CTCF/cohesin ChIA-PET experiment, the "E" and "P2P" files will actually be the same; I just make one of each because I concatenate the *pairs files for each enhancer list and each method (P2P/E) together to create the final locus definitions for that enhancer list and that method, and creating an "E" file even for the non-CTCF/cohesin experiments (equivalent to a P2P file for these experiments) allows me to have a consistent naming scheme for concatenating the files that belong together).  The *.pairs files for a given enhancer list can be concatenated together to get the locus definition corresponding to that enhancer list; the *.pipeline files do this concatenation to create the "*.all.E.pairs" and "*.all.P2P.pairs" files mentioned earlier in this section.

## Evaluation of locus definitions

Next, there are a variety of ways in which the new locus definitions can be compared to the old locus definitions (note: the old locus definitions are included in the ```current_definitions``` file; these are the "TSS +/- 5kb and 2kb enhancers based on FANTOM5" definitions.  If you wish to run my analyses using a different set of locus definitions, you can do so.  Just take the locus list of interest, make sure it's in the same format as shown in the ```current.ldef``` file, and change the name to ```current.ldef``` so that it's compatible with all my scripts.  Also, you'll need to include a list of the regions in this locus definition that correspond to genes (format: chrom, start, end, geneid; tab-separated, no header; see ```genes.txt```) and to enhancers (format: chrom, start, end; tab-separated, no header; see ```enhancers.txt```).  For the gene file and the enhancer file, change the names to ```genes.txt``` and ```enhancers.txt```; essentially, you're just replacing the files currently in the ```current_definitions``` directory).

### Plot the overlap between the old and the new enhancer lists

The first comparison that can be done is checking how many of the old enhancers overlap with the new enhancers (i.e., what fraction of the old enhancers overlap with at least one enhancer from the new enhancer list, and what fraction of base pairs from the old enhancers are covered in the new enhacers).  You can make these plots by changing into the ```overlap_old_enhancers_with_new``` directory and running:

	make plot

### Overlap all the ENCODE ChIP-seq peaks with the locus lists
To compare the fraction of ChIP-seq peaks caught by the new lists with the fraction caught by the old lists, ```cd``` into the ```peak_catching``` directory and run:

	make pipeline
	make run

This will overlap each ChIP-seq experiment's peaks with the locus definitions and print out a few basic statistics.  You can run

	make master_list

to concatenate the results into a single file that can be used to plot things in R

### Compare chipenrich p values using different locus lists
You can additionally run chipenrich using the different locus lists and plot negative log10 p values against each other.  To do this, ```cd``` into the ```chipenrich``` directory and run:

	make pipeline
	make run
	make pairs_plots

The ```make run``` can take a long time to complete since it's running chipenrich on many different ChIP-seq experiments and locus lists; therefore you'll probably want to do this step with ```nohup```.
