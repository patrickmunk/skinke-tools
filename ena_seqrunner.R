# ENA Sequence runner
# For submitting sequence run to ENA

library(dplyr)

# Some settings that are likely to stay static for all the samples/runs in a project
projectSettings = list(STUDY = "PRJEB70907", INSTRUMENT = "Illumina NovaSeq 6000",
                       LIBRARY_SOURCE = "METAGENOMIC", INSERT_SIZE = 300,
                       LIBRARY_SELECTION = "RANDOM", LIBRARY_STRATEGY = "WGS")


ENAreadManifestWriter = function(ena_sample_files_tbl, out_dir, project_settings) {
  # input: a df with sample and full paths to reads
  # output: A manifest file for each samples' reads to submit
  num_samples = nrow(ena_sample_files_tbl)
  col1 = c("STUDY", "SAMPLE", "NAME", "INSTRUMENT", "INSERT_SIZE", "LIBRARY_SOURCE",
           "LIBRARY_SELECTION", "LIBRARY_STRATEGY", "FASTQ", "FASTQ")
  for (i in 1:num_samples) {
     col2 = c(project_settings$STUDY, 
            ena_sample_files_tbl$ACCESSION[i], 
            ena_sample_files_tbl$ALIAS[i],  
            project_settings$INSTRUMENT, 
            project_settings$INSERT_SIZE, 
            project_settings$LIBRARY_SOURCE, 
            project_settings$LIBRARY_SELECTION, 
            project_settings$LIBRARY_STRATEGY,
            ena_sample_files_tbl$read1[i], 
            ena_sample_files_tbl$read2[i])
    result = cbind(col1, col2)
    # Now write the file
    outfile_path_name = file.path(out_dir, paste(ena_sample_files_tbl[i,"ALIAS"], "manifest", sep = "."))
    write.table(result, outfile_path_name, sep = " ", col.names = F, row.names = F, quote = F)
  }
}


# Read ENA-provided sample accessions

sampleAccFile = "C:/Users/pmun/Downloads/Webin-accessions-2023-12-15T16_59_41.688Z.txt"

sampleAccTbl = read.delim(sampleAccFile)
sampleAccTbl = sampleAccTbl[1:nrow(sampleAccTbl)-1,]

# Read in the information on file names that should be used for matching
readFilesFofn = read.delim("C:/Users/pmun/Downloads/fonfn.txt", h = F) %>%
  rename(fullpath = V1) %>%
  arrange(fullpath) %>%
  mutate(filename = basename(fullpath)) %>%
  mutate(runName = gsub(".gz", "", filename)) %>% 
  mutate(runName = tools::file_path_sans_ext(runName))

# Find which read files match what alias and join
MatchAliasToFiles = function(sample_acc_tbl, reads_fofn) {
  result = data.frame(sample_acc_tbl, read1 = NA, read2 = NA)
  for (s in 1:nrow(sample_acc_tbl)) {
    file_hits = grepl(sample_acc_tbl[s,3], reads_fofn$runName)
    result[s,"read1"] =  reads_fofn[file_hits,"fullpath"][1]
    result[s,"read2"] =  reads_fofn[file_hits,"fullpath"][2]
    #print(file_hits)
  }
  return(result)
}

sampleFilesJoined = MatchAliasToFiles(sampleAccTbl, readFilesFofn)

# Subset to just the ena_sample_files_tbl format
sampleFilesJoined = sampleFilesJoined %>%
  dplyr::select(ACCESSION, read1, read2, ALIAS)

# Now write the manifest files
ENAreadManifestWriter(sampleFilesJoined, "C:/Users/pmun/Downloads/manifests/", projectSettings)
