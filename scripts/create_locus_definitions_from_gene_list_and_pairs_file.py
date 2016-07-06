#!/bin/python

from optparse import OptionParser
import re

### Handle command line options
parser = OptionParser()
parser.add_option("--gene_list", action = "store", dest = "genes", type = str, help = "[Required] The list of genes; chromosome, start, end, geneid (header optional)")
parser.add_option("--pairs_list", action = "store", dest = "pairs", type = str, help = "[Required] The list of gene-enhancer pairs; geneid, enhancerid (enhancerid = chromosome:start:end) (header optional)")


(options, args) = parser.parse_args()

### Read in the gene list
f = open(options.genes)

genes = {}

for line in f:
	if "chrom" in line:
		# it's a header; skip
		continue

	line = line.rstrip()
	(chromosome, start, end, geneid) = line.split("\t")

	genes[geneid] = chromosome + ":" + start + ":" + end

f.close()

# Output the header for the locus definitions:
print "chrom\tstart\tend\tgeneid"

# Output the gene loci for the locus definitions:
for (geneid, val) in genes.items():
	(chromosome, start, end) = val.split(":")
	print chromosome + "\t" + start + "\t" + end + "\t" + geneid


### Read in the pairs list
f = open(options.pairs)

for line in f:
	if "gene" in line:
		# it's a header
		continue

	line = line.rstrip()
	(geneid, enhancerid) = line.split("\t")
	(chromosome, start, end) = enhancerid.split(":")
	print chromosome + "\t" + start + "\t" + end + "\t" + geneid

f.close()
