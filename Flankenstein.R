#!/usr/bin/env Rscript
# This is a script for extracting flanking regions of a fixed length
# from a collecting of genes
# Author : Patrick Munk

# Setup
bedtoolspath = "/services/tools/bedtools/2.28.0/bin/bedtools"
kmapath = "/home/projects/cge/people/pmun/bin/kma-1.3.3/kma"

# R libraries
library(optparse)

# Make argument list
option_list = list(
  make_option(c("-b", "--blastfile"), action="store", default=NA, type='character',
              help="Tabular BLAST output indicating which genes are found where on which contig"),
  make_option(c("-c", "--contig"), action="store", default=NA, type="character",
              help="Contigs in fasta format that have corresponding BLAST output"),
  make_option(c("-n", "--namedFeature"), action="store", default="ALL", type='character',
              help="Exact name of feature to extract intervals for [default %default]"),
  make_option(c("-l", "--minLen"), action="store", default=1000,
              help="Exact name of feature to extract intervals for [default %default]"),
  make_option(c("-m", "--maxLen"), action="store", default=1000,
              help="What should be the maximum flanking region to pad with? [default %default]"),
  make_option(c("-t", "--threads"), action="store", default=1,
              help="How many threads should be used at once? [default %default]"),
  make_option(c("-o", "--outputName"), action="store", default="FlankedGenes", type="character",
              help="Output file prefix name [default %default]"))

# Parge arguments to R
opt <- parse_args(OptionParser(option_list=option_list))

# Make output directory if it does not already exist
fullpathoutdir = file.path(getwd(), opt$outputName)
dir.create(fullpathoutdir, showWarnings = FALSE)
print(paste("Saving results to:", fullpathoutdir))

# Read in table with BLAST results
blastdata = read.delim(file=opt$blastfile, h = F)

# Subset dataset if user has specified a specific feature of interest
if (opt$namedFeature != "ALL") {
  blastdata = blastdata[blastdata$V1 == opt$namedFeature,]
}

# blastdata = blastdata[!duplicated(blastdata[,c(1:6)]),]

# Order the features according to priority for easier filtering on first occurence
# Select first based on percent identity, then by common name, then by accession number
blastdata = blastdata[order(blastdata$V2, blastdata$V1, blastdata$V9, decreasing = T),]

# Extract features for BED file
contigpos = blastdata$V7
startpos = sapply(strsplit(as.character(contigpos), "..", fixed=T), `[`, 1)
endpos = sapply(strsplit(as.character(contigpos), "..", fixed=T), `[`, 2)
contigname = blastdata$V6
featurename = blastdata$V1
contigInfo = strsplit(as.character(contigname), "_NODE_", fixed=T)
sampleName = sapply(contigInfo, `[`, 1)
contigLen = sub(".*_length_ *(.*?) *_cov_.*", "\\1", blastdata$V6)

# Make a BED file for the selected feature
bedfile = data.frame(contigname, startpos, endpos, featurename)

# Filtering of features
# Disallow that any two features share start or end positions on same contig. Choose feature
bedfile = bedfile[!duplicated(bedfile[,c("contigname","startpos")]),]
bedfile = bedfile[!duplicated(bedfile[,c("contigname","endpos")]),]

bedfileName = paste(file.path(fullpathoutdir, opt$outputName), ".bed", sep = "")
write.table(bedfile, file = bedfileName, quote = F, row.names = F, col.names = F, sep = "\t")

# Make a genome file for the selected feature
genomefile = data.frame(contigname, contigLen)
genomefileName = paste(file.path(fullpathoutdir, opt$outputName), ".genome.txt", sep = "")
write.table(genomefile, file = genomefileName, quote = F, row.names = F, col.names = F, sep = "\t")

# Extend the intervals with user-specified lengths
slopbedfileName = paste(file.path(fullpathoutdir, opt$outputName), ".slop.bed", sep = "")
system(paste(bedtoolspath, "slop -i", bedfileName, "-g", genomefileName, "-b", opt$maxLen, ">", slopbedfileName, sep = " "))

# Read in extended and original intervals for compairson
slop = read.delim(file = slopbedfileName, h = F)
genes = read.delim(file = bedfileName, h = F)

# Only retain BED intervals that got extended by at least the minimum length up- and down-stream
geneWithFlank = (slop$V2 <= genes$V2 - opt$minLen & slop$V3 >= genes$V3 + opt$minLen)
fullflanks = slop[geneWithFlank,]

# Overwrite BED file with gene intervals with new version only containing genes with proper flanks
bedfileGene = bedfile[geneWithFlank,]
write.table(bedfileGene, file = "genes.temp.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# Export the flank file to the user
print(paste("Retaining", nrow(fullflanks), "genes out of", nrow(slop), "genes.", sep = " "))
write.table(fullflanks, file = "fullflanks.temp.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# Export fasta with genes and flanks, but the flanks masked
maskfilename = paste(file.path(fullpathoutdir, opt$outputName), ".masked.fa", sep="")
system(paste(bedtoolspath, "maskfasta -fullHeader -fi", opt$contig, "-bed", bedfileName, "-fo ", maskfilename, sep = " "))

# Export fasta with masked genes and fixed flank lengths
maskFlankfilename = paste(file.path(fullpathoutdir, opt$outputName), ".GeneMaskedFlank.fa", sep="")
system(paste(bedtoolspath, "getfasta -bed fullflanks.temp.bed -fi" , maskfilename, "-fo", maskFlankfilename, sep = " "))

# Export fasta with fixed flank length and no masking done
flankfilename = paste(file.path(fullpathoutdir, opt$outputName), ".GeneFlank.fa", sep="")
system(paste(bedtoolspath, "getfasta -bed fullflanks.temp.bed -fi" , opt$contig, "-fo", flankfilename, sep = " "))

# Export fasta of genes that have corresponding flank files
genefilename = paste(file.path(fullpathoutdir, opt$outputName), ".Gene.fa", sep="")
system(paste(bedtoolspath, "getfasta -bed genes.temp.bed -fi" , opt$contig, "-fo", genefilename, sep = " "))

# Clean up temp files
system("rm fullflanks.temp.bed")
system("rm genes.temp.bed")
system(paste("rm",maskfilename))

# File suffixes to remove
files2remove = c(".bed", ".genome.txt", ".slop.bed", ".masked.fa.fai")
files2remove = file.path("rm" ,paste(fullpathoutdir, files2remove))
# TBD
# system(paste("echo" ,paste(opt$outputName, files2remove, sep="")))

# Make KMA indeces for the gene and masked flank fasta files
kmaindexcommand = paste(kmapath, "index", "-Sparse -", "-i", maskFlankfilename, "-o", maskFlankfilename)
system(kmaindexcommand)

kmaindexcommand2 = paste(kmapath, "index", "-Sparse -", "-i", genefilename, "-o", genefilename)
system(kmaindexcommand2)

# Use KMA to calculate kmer overlap coefficients for the flanks and genes
kmadistcommand = paste(kmapath, "dist", "-t", 1, "-d", 2048, "-t_db", maskFlankfilename) 
system(kmadistcommand)

kmadistcommand2 = paste(kmapath, "dist", "-t", 1, "-d", 2048, "-t_db", genefilename) 
system(kmadistcommand2)

# Example
# ./Flankenstein.R -b contigs_resfinder.txt -c resfinder_contigs_allsamples.fa -n "tet(X)" -l 100 -o ThisIsATest
# ./Flankenstein.R -b contigs_resfinder.txt -c resfinder_contigs_allsamples.fa -l 100 -o TestAll