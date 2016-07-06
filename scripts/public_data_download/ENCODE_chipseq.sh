#!/bin/bash

################################################################################
# The purpose of this script is to download ENCODE Uniform TFBS .narrowPeak and
# ENCODE histone .broadPeak files
################################################################################

# Create the folder and go to there
# NOTE: ~/espresso is a specific symlink in my (rcavalca) home directory
# if you're not me, you'll likely need to change the paths in mkdir and cd
# but everything else should be fine.
mkdir ENCODE_chipseq
cd ENCODE_chipseq

# Get the TFBS file list from ENCODE (hopefully this URL is stable for a while)
# and rename it to something clearer.
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/files.txt
mv files.txt tfbs_files.txt

# Get the histone file list from ENCODE (hopefully this URL is stable for a while)
# and rename it to something clearer.
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/files.txt
mv files.txt histone_files.txt

# Use awk to create a bunch of wgets in a bash script for the narrowPeak
# and broadPeak files only. The directory contains raw-er data, but we
# don't need it.
awk '{print "wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/"$1}' <(cut -f 1 tfbs_files.txt | grep narrowPeak) > tfbs_wget.sh
awk '{print "wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/"$1}' <(cut -f 1 histone_files.txt | grep broadPeak) > histone_wget.sh

# Go and get them
bash tfbs_wget.sh
bash histone_wget.sh

# Create an annotation file with columns for filename, cell line, treatment, and antibody

# TFBS
	# Get the file name
	cut -f 1 tfbs_files.txt > tmp_files.txt
	# Get cell line
	cut -f 2 tfbs_files.txt | cut -d ';' -f 6 | cut -d '=' -f 2 > tmp_lines.txt
	# Get treatment
	cut -f 2 tfbs_files.txt | cut -d ';' -f 7 | cut -d '=' -f 2 > tmp_treatments.txt
	# Get antibody / tf
	cut -f 2 tfbs_files.txt | cut -d ';' -f 8 | cut -d '=' -f 2 > tmp_antibodies.txt
	cut -f 2 tfbs_files.txt | cut -d ';' -f 8 | cut -d '=' -f 2 | cut -d '_' -f 1 > tmp_tf.txt

	# Stitch it together with paste
	paste tmp_files.txt tmp_lines.txt tmp_treatments.txt tmp_tf.txt tmp_antibodies.txt > tmp_annotated.txt
	grep narrowPeak tmp_annotated.txt > tfbs_annotated_files.txt
	rm tmp_*

# Histones
	# Get the file name
	cut -f 1 histone_files.txt > tmp_files.txt
	# Get cell line
	cut -f 2 histone_files.txt | cut -d ';' -f 7 | cut -d '=' -f 2 > tmp_lines.txt
	# Get treatment
	cut -f 2 histone_files.txt | cut -d ';' -f 8 | cut -d '=' -f 2 > tmp_treatments.txt
	# Get antibody / tf
	cut -f 2 histone_files.txt | cut -d ';' -f 9 | cut -d '=' -f 2 > tmp_antibodies.txt
	cut -f 2 histone_files.txt | cut -d ';' -f 9 | cut -d '=' -f 2 | cut -d '_' -f 1 > tmp_tf.txt

	# Stitch it together with paste
	paste tmp_files.txt tmp_lines.txt tmp_treatments.txt tmp_tf.txt tmp_antibodies.txt > tmp_annotated.txt
	grep broadPeak tmp_annotated.txt > histone_annotated_files.txt
	rm tmp_*
