#!/usr/bin/env Rscript
#This script is used to merge a collection of mapstat files into a catmapstat file

# Generates python-like options used when the Rscript from a terminal. Requires the R-package 'optparse'
library("optparse")
option_list = list(
  make_option(c("-m", "--mapstatdir"), type="character", default=NULL,
              help="Directory with *.mapstat files", metavar="character"),
        make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="Output catmapstat file name prefix [default= %default]", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Determine the prefix to be used for output files. Date or user-defined
if (is.null(opt$prefix)){
        prefix = format(Sys.time(), "%Y%m%d%H%M")
        } else {
                prefix <- opt$prefix
        }

# Check if mandatory dir is supplied
if (is.null(opt$mapstatdir)){
  stop("A directory with mapstat file needs to be supplied (-m)", call.=FALSE)
}

#setwd(opt$mapstatdir)


################# PROGRAM START #################

# Read in file names from wd
#file_list <- list.files()
file_list = list.files(opt$mapstatdir)

# Select just the mapstat files in the working directory
file_list = file_list[grepl(".mapstat",file_list, fixed=T)]
file_paths = file.path(opt$mapstatdir, file_list) 

# Initialize output 
catmapstat = NULL

# Iteratively read in every file, add sample name and concatenate results
#for (i in file_paths) {
#	mapstat = read.delim(i, h=T, sep="\t", skip=6, check.names=F)
#  print(nrow(mapstat))
#	sampleName = gsub(".mapstat", "", i)
#	mapstat = data.frame(sample = sampleName, mapstat)
#	catmapstat = rbind(catmapstat, mapstat)
#}

for (i in file_paths) {
  mapstat = tryCatch({
    read.delim(i, h=T, sep="\t", skip=6, check.names=F)
    }, error = function(err) {
      print(paste("Could not read a table. Empty? Check output!",err))
      })
      sampleName = gsub(".mapstat", "", i)
      mapstat = data.frame(sample = sampleName, mapstat)
      catmapstat = rbind(catmapstat, mapstat)
}

# Fix the bad column name
colnames(catmapstat) = gsub("X..refSequence", "refSequence", colnames(catmapstat))

# Summarize output
print(paste("Writing a catmapstat file with", 
  length(unique(catmapstat$refSequence)), 
  "features in", 
  length(unique(catmapstat$sample)), 
  "samples."))

# Determine output name and output file
outname = paste(opt$prefix, ".cms", sep="") 
write.table(catmapstat, file=outname, sep="\t", quote=F, row.names = F)
