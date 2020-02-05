## The script visualize the BED files by mapping the reads on the human genome scheme.
## This script doesn't perform any alignments, this is merely the visualization step!
## Required input:
##  1. Text file containing chromosome length info (GRCh38 and hg19 are available in the github repo)
##  2. Text file containing centromeres coordinates (also can be found in the repo)
##  3. Sorted BED file - reads aligned to the genome
## Output:
##  1. PDF files containing the plots, showing how the alignments are located across the human genome
##  2. Report file with logs.
## To make the script work you need to set:
##  1. The working directory (line 17)
##  2. A folder containing the BED files (line 18).
##     If left empty or NA, the script will search for BED files in the working directory or any downstream folder
##  3. Put the chromosome length and centromeres coordinates files into the working directory or any downstream folder

# Directories -------------------------------------------------------------
##################################################################

#####         LAM #####
#input
#folder.with.sorted.bed <- "/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/Mishas_Demultiplexed_LAM/bed_files"
#output
#setwd("/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/Mishas_Demultiplexed_LAM/bed_plots")

#################       Tn5  ###################
# input
folder.with.sorted.bed <- "/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/P_Eric4_Tn5/bed_files"
# output
setwd("/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/P_Eric4_Tn5/bed_plots")



# Logging -----------------------------------------------------------------
report.file <- paste("lam_htgts_report_", format(Sys.time(), "%d%m%y_%H%M"), ".txt", sep = "")
write(paste("LAM-HTGTS Mapping script report ", format(Sys.time(), "%d%m%y_%H%M"), sep = ""),
      file = report.file, append = TRUE)
write(paste("Working directory selected: ", getwd(), sep = ""), file = report.file, append = TRUE)


# Libraries ---------------------------------------------------------------
libs <- c('tidyr', 'ggplot2', 'ggrepel', 'reshape2', 'stringr', 'Biostrings', 'reader', 'lazyeval')
library(dplyr)
lapply(libs, require, character.only = TRUE) %>% unlist() %>%
  sum() -> libs.installed
if (libs.installed == length(libs)) {
  print("All libraries installed successfully")
}else{
  stop("Some libraries were not installed properly")
}
rm(libs, libs.installed)

# Functions ---------------------------------------------------------------

### what does this do
search_files <- function(folder = NA, file.or.pattern, log.file = report.file){
  if(folder == "" | is.na(folder)){
    write(paste("No folder specified for", file.or.pattern, "file(s), looking for the files in the working directory.", sep = " "),
          file = log.file, append = TRUE)
    files <- list.files(path = getwd(), pattern = file.or.pattern)
    if(length(files) == 0){
      write(paste("No", file.or.pattern, "files found in the working directory, proceeding to recursive search.", sep = " "),
            file = log.file, append = TRUE)
      files <- list.files(path = getwd(), pattern = file.or.pattern, recursive = TRUE)
      if(length(files) == 0){
        write("Recursive search did not give any result, execution stopped.")
        stop(paste("No", file.or.pattern, "files found, execution aborted", sep = " "))
      }else{
        write(paste(length(files), file.or.pattern, "file(s) found in downstream folders.\n", sep = " "),
              file = log.file, append = TRUE)
      }

    }else{
      write(paste(length(files), file.or.pattern, "file(s) found in the working directory.\n", sep = " "),
            file = log.file, append = TRUE)
    }
  }else{
    if(dir.exists(folder)){
      files <- list.files(path = folder, pattern = file.or.pattern)
      if(length(files) == 0){
        write(paste("No", file.or.pattern, "files found in the specified folder, proceeding to recursive search.", sep = " "),
              file = log.file, append = TRUE)
        files <- list.files(path = getwd(), pattern = file.or.pattern, recursive = TRUE)
        if(length(files) == 0){
          write("Recursive search did not give any result, execution stopped.\n")
          stop(paste("No", file.or.pattern, "files found, execution aborted", sep = " "))
        }else{
          write(paste(length(files), file.or.pattern, "file(s) found in downstream folders.\n", sep = " "),
                file = log.file, append = TRUE)
        }
      }else{
        write(paste(length(files), file.or.pattern, "file(s) found in the specified folder.\n", sep = " "),
              file = log.file, append = TRUE)
      }
    }else{
      write("The specified folder doesn't exist.\n",
            file = log.file, append = TRUE)
    }
  }
  return(files)
}

lookup <- function(input, table){
  return(table$output[which(str_detect(input, table$input))])
}

### what does this do
calculate_offsets <- function(df, method = 'middle'){
  if (method == 'middle'){
    bins_center <- bins + binwidth/2
    df$AX <- bins_center[df$binX]
  }

  df <- unite(df, 'binXchr', c('binX', 'chromosome'), sep = ";")

  df <- left_join(df, df %>% group_by(binXchr) %>% summarise(cnt = n()), by = 'binXchr')

  if (method == 'mean'){
    df <- left_join(df, df %>% group_by(binXchr) %>% summarise(meanAX = mean(AX)), by = 'binXchr')
    df$AX <- df$meanAX
    df <- select(df, -'meanAX')
  }

  df <- df %>% group_by(binXchr) %>% mutate(elementrank = min_rank(ins.id))

  df$nlevel <- ceiling(-0.5 + 0.5*sqrt(8*df$cnt+1))
  df$level <- -0.5 + 0.5*sqrt(8*df$elementrank+1)
  df$levelr <- ceiling(df$level)

  df <- df %>% group_by(binXchr, levelr) %>% mutate(elementposinlevel = min_rank(level))

  df <- unite(df, 'binXchrLevelR', c('binXchr','levelr'), sep = ":")
  dfg <- df %>%
    group_by(binXchrLevelR) %>%
    summarise(elementsonlevel = max(elementposinlevel))
  df <- left_join(df, dfg, by = 'binXchrLevelR')
  df <- separate(df, 'binXchrLevelR', c('binXchr','levelr'), sep = ":")
  df$levelr <- as.numeric(df$levelr)

  df$offsetY <- df$levelr - 1
  df$offsetX <- 2*df$elementposinlevel - df$elementsonlevel - 1
  df <- select(df, -c('levelr','nlevel','level','elementsonlevel','elementposinlevel',
                      'elementrank','cnt'))
  df <- separate(df, 'binXchr', c('binX', 'chromosome'), sep = ";")
  return(df)
}

count.overlap <- function(subject.starts, subject.ends, query.start, query.end){
  subject.starts <- as.numeric(subject.starts)
  subject.ends <- as.numeric(subject.ends)
  query.start <- as.numeric(query.start)
  query.end <- as.numeric(query.end)
  if(sum(subject.starts > subject.ends, na.rm = TRUE)){
    stop("Some of the subject ranges you provided are flipped")
  }
  if(query.start > query.end){
    stop("Query range you provided is flipped")
  }
  count.overlap <- sum(query.start >= subject.starts & query.start <= subject.ends|
                         query.end >= subject.starts & query.end <= subject.ends|
                         query.start <= subject.starts & query.end >= subject.ends)

  return(count.overlap)
}

# File search -------------------------------------------------------------
#chrlength.file <- search_files(folder = getwd(), file.or.pattern = "chrlength")

chrlength.file <- read.table("/home/edanner/workspace/uditas/Pipeline2Plotting/GRCh38_chrlength.txt",
                             header = TRUE, stringsAsFactors = FALSE)
# ERIC PC or LEBMIH ASUS
if (Sys.info()["nodename"] != "LEBMIHASUS") {
  print("Hi, sexy stranger!")
  offtargets <- read.table("/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/offtargets_hg38-unknownLoc.txt",
                           header = TRUE, stringsAsFactors = FALSE)
} else {
  print("Hi, misha!")
  offtargets <- read.table("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics/eric_faa/offtargets_hg38-unknownLoc.txt",
                           header = TRUE, stringsAsFactors = FALSE)
  setwd("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics")
}


## The program uses just one file with chromosome length, it chooses the most recent one
# if(length(chrlength.file) > 1){
#   file.info.table <- bind_cols(data.frame(filename = chrlength.file), file.info(chrlength.file))
#   file.info.table <- arrange(file.info.table, desc(atime))
#   chrlength.file <- file.info.table$filename[1]
#   warning(paste("More than one file with chromosome length found, using the most recent one, ", chrlength.file, ".", sep = ""))
#   write(paste("More than one file with chromosome length found, using the most recent one, ", chrlength.file, ".\n", sep = ""),
#         file = report.file, append = TRUE)
#   rm(file.info.table)
# }else{
#   write(paste("Chromosome length file used: ", chrlength.file, ".\n", sep = ""),
#         file = report.file, append = TRUE)
# }

centromeres <- read.table("/home/edanner/workspace/uditas/Pipeline2Plotting/hg38_centromeres.txt", header = TRUE, stringsAsFactors = FALSE)
## The program uses just one centromere file, it chooses the most recent one
# if(length(centromere.file) > 1){
#   file.info.table <- bind_cols(data.frame(filename = centromere.file), file.info(centromere.file))
#   file.info.table %>% mutate(filename = as.character(filename)) -> file.info.table
#   file.info.table <- arrange(file.info.table, desc(atime))
#   centromere.file <- file.info.table$filename[1]
#   warning(paste("More than one file with centromere coordinates found, using the most recent one, ", centromere.file, ".", sep = ""))
#   write(paste("More than one file with centromere coordinates found, using the most recent one, ", centromere.file, ".\n", sep = ""),
#         file = report.file, append = TRUE)
#   rm(file.info.table)
# }else{
#   write(paste("Centromere coordinate file used: ", centromere.file, ".\n", sep = ""),
#         file = report.file, append = TRUE)
# }

# loading and concationating all the bed files for processing
# we then remove the extra comments from the bowtie2 as they were filtered in python for primary reads, high quality
bed.files <- search_files(folder = folder.with.sorted.bed, file.or.pattern = ".sorted.bed")
bed <- data.frame()
for (i in bed.files) {
  filename <- paste0(folder.with.sorted.bed, "/", i)
  bed.temp <- read.table(filename, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  bed.temp$filename <- i
  bed <- bind_rows(bed, bed.temp)
}
bed.files <- bed %>% select(V1:V12, filename)
rm(bed, bed.temp)
#
# bed.files <- read.table(paste0(folder.with.sorted.bed, "/N707_N505.sorted.bed"), stringsAsFactors = FALSE,
#                         fill = TRUE)
#
colnames(bed.files) <- c("chrom", "start", "end", "read.id", "quality", "strand",
                         "flag", "cigar", "asterisk", "zero1", "zero2", "seq",
                         "filename")
mc <- bed.files %>% filter(chrom == "pE038_MC")
bed.files <- bed.files %>% filter(!str_detect(chrom, "_") & chrom != "pE038_MC")
# Quality based filtering
threshold <- 25
bed.fileslowQ <- bed.files %>% dplyr::filter(quality < threshold)
bed.files <- bed.files %>% dplyr::filter(quality >= threshold)

# -------------------- Counting the overlap of the Cas9-offtarget sites with the reads (Cas9 driven events)  -------------------

# this measures for all the off targets in our table
offtarget.overlap <- c()
for(i in seq(nrow(offtargets))){
  bed.files.tmp <- bed.files %>% dplyr::filter(chrom == offtargets$chrom[i])
  bed.files.tmp <- bed.files.tmp %>%
    mutate(center = (start + end)/2,
           dist.to.off = center - (offtargets$start[i] + offtargets$end[i])/2,
           dist.to.off = abs(dist.to.off))
  overlaps <- sum(bed.files.tmp$dist.to.off <= 500)
  # num <- sum((bed.files.tmp$start >= offtargets$start[i] & bed.files.tmp$end <= offtargets$end[i])|
  #              (bed.files.tmp$start <= offtargets$start[i] & bed.files.tmp$end >= offtargets$end[i])|
  #              (bed.files.tmp$start >= offtargets$start[i] & bed.files.tmp$start <= offtargets$end[i])|
  #              (bed.files.tmp$end >= offtargets$start[i] & bed.files.tmp$end <= offtargets$end[i]))
  offtarget.overlap <- c(offtarget.overlap, overlaps)
}
offtargets$overlap <- offtarget.overlap
offtargets <- offtargets %>%
  arrange(desc(overlap))
# need to check off target table manually by opening dataframe and scroll to right and look at "overlap" column count

# Filtering of the off-target table to show the top_n on the plot later
top_n <- 10
offtargets <- offtargets %>%
  arrange(desc(mitOfftargetScore)) %>%
  dplyr::slice(1:top_n)



# Library of values for check step ----------------------------------------

## All names for the columns of interest should be listed here
chromosome.columns <- c("Chromosome", "chromosome", "Chrom", "chrom", "chr", "Chr")
length.columns <- c("Length.bp", "length.bp", "Length", "length", "Len.bp", "len.bp", "Len", "len")
start.columns <- c("Start", "start", "chromStart", "chrStart", "chromstart", "chrstart", "ChromStart",
                   "Chromstart", "Chrstart", "ChrStart")
end.columns <- c("End", "end", "chromEnd", "chrEnd", "chromend", "chrend", "ChromEnd",
                 "Chromend", "Chrend", "ChrEnd")
## Look-up table for chromosome names standartization
chrom.names.lookuptable <- data.frame(input = c(paste(as.character(c(1:22)),"$", sep = ""), "[Xx]$", "[Yy]$", "[mM][tT]{0,1}$"),
                                      output = paste("chr", c(as.character(c(1:22)), "X", "Y", "M"), sep = ""),
                                      stringsAsFactors = FALSE)
chrom.names.lookuptable$input <- paste("^", chrom.names.lookuptable$input, "|^[a-zA-Z]{3,10}",
                                       chrom.names.lookuptable$input, sep = "")



#---------------------------------- Data import and check ---------------------------------------------------

write("\n#### DATA CHECK ####\n", file = report.file, append = TRUE)
# Detection of the assembly in the chromosome length file
# chrlength.file.last.piece <- str_split(chrlength.file, "/") %>% unlist()
# chrlength.file.last.piece <- tail(chrlength.file.last.piece, n = 1)
# if(any(str_detect(chrlength.file.last.piece, c("hg19", "GRCh37", "Grch37", "grch37", "GRCH37")))){
#   chrlength.file.assembly <- "hg19"
#   write("hg19 assembly used for chromosome length data", file = report.file, append = TRUE)
# }else if(any(str_detect(chrlength.file.last.piece, c("hg38", "GRCh38", "Grch38", "grch38", "GRCH38")))){
#   chrlength.file.assembly <- "hg38"
#   write("hg38 assembly used for chromosome length data", file = report.file, append = TRUE)
# }else{
#   warning("Unable to determine the assembly used for chromosome length data")
#   write("Unable to determine the assembly used for chromosome length data", file = report.file, append = TRUE)
# }
# rm(chrlength.file.last.piece)
# # Detection of the assembly in the centromere coordinates file
# centromere.file.last.piece <- str_split(centromere.file, "/") %>% unlist()
# centromere.file.last.piece <- tail(centromere.file.last.piece, n = 1)
# if(any(str_detect(centromere.file.last.piece, c("hg19", "GRCh37", "Grch37", "grch37", "GRCH37")))){
#   centromere.file.assembly <- "hg19"
#   write("hg19 assembly used for centromere data", file = report.file, append = TRUE)
# }else if(any(str_detect(centromere.file.last.piece, c("hg38", "GRCh38", "Grch38", "grch38", "GRCH38")))){
#   centromere.file.assembly <- "hg38"
#   write("hg38 assembly used for centromere data", file = report.file, append = TRUE)
# }else{
#   warning("Unable to determine the assembly used for centromere data")
#   write("Unable to determine the assembly used for centromere data", file = report.file, append = TRUE)
# }
# rm(centromere.file.last.piece)
# # Checking whether the assembly used is the same in both files
# if(chrlength.file.assembly != centromere.file.assembly){
#   warning("Chromosome length data and centromeres coordinates originate from different assemblies!")
#   write("Chromosome length data and centromeres coordinates originate from different assemblies!",
#         file = report.file, append = TRUE)
# }
# rm(chrlength.file.assembly, centromere.file.assembly)

# Reading the chromosome length data
chrlen <- chrlength.file
# Checking the presence of chromosome and length columns
if(!any(colnames(chrlen) %in% chromosome.columns)){
  stop("No chromosome column found in chromosome length file!")
  write("No chromosome column found in chromosome length file!", file = report.file, append = TRUE)
}
if(!any(colnames(chrlen) %in% length.columns)){
  stop("No length column found in chromosome length file!")
  write("No length column found in chromosome length file!", file = report.file, append = TRUE)
}
# Checking the number of observations
if(!nrow(chrlen) %in% c(24,25)){
  stop("The number of observations in chromosome length data is not equal to number of chromosomes in human genome!")
  write("The number of observations in chromosome length data is not equal to number of chromosomes in human genome!",
        file = report.file, append = TRUE)
}

# Reading the centromeres coordinates data
#centromeres <- read.table(file = centromere.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# Checking the presence of chromosome and coordingates columns
if(!any(colnames(centromeres) %in% chromosome.columns)){
  stop("No chromosome column found in centromere file!")
  write("No chromosome column found in centromere file!", file = report.file, append = TRUE)
}
if(!any(colnames(centromeres) %in% start.columns)){
  stop("No start coordinate column found in centromere file!")
  write("No start coordinate column found in centromere file!", file = report.file, append = TRUE)
}
if(!any(colnames(centromeres) %in% end.columns)){
  stop("No end coordinate column found in centromere file!")
  write("No end coordinate column found in centromere file!", file = report.file, append = TRUE)
}
# Checking the number of observations
if(nrow(centromeres) > 24){
  warning("There are more observations in centromeres coordinates data than chromosomes in human nuclear genome!
          Perhaps, you used unfiltered gaps table.")
  write("There are more observations in centromeres coordinates data than chromosomes in human nuclear genome!
        Perhaps, you used unfiltered gaps table.",
        file = report.file, append = TRUE)
}else if(nrow(centromeres) < 24){
  stop("There are less observations in centromeres coordinates data than chromosomes in human nuclear genome!")
  write("There are less observations in centromeres coordinates data than chromosomes in human nuclear genome!",
        file = report.file, append = TRUE)
}



#--------------------------------- Data preprocessing ------------------------------------------------------
## Column names to standard
chrlen$chrom <- chrlen[,which(colnames(chrlen) %in% chromosome.columns)]
chrlen$length <- chrlen[,which(colnames(chrlen) %in% length.columns)]

centromeres$chrom <- centromeres[,which(colnames(centromeres) %in% chromosome.columns)]
centromeres$start <- centromeres[,which(colnames(centromeres) %in% start.columns)]
centromeres$end <- centromeres[,which(colnames(centromeres) %in% end.columns)]

# Chromosome length data: filtering out all columns except chromosome name and length
chrlen <- select(chrlen, chrom, length)
# Centromeres data: filtering out all cols except chrom, start and end
centromeres <- select(centromeres, chrom, start, end)

## Variables to standard
# In case length or coordinate values contain extra characters:
chrlen$length <- as.integer(str_replace_all(chrlen$length, "\\D", ""))
centromeres$start <- as.integer(str_replace_all(centromeres$start, "\\D", ""))
centromeres$end <- as.integer(str_replace_all(centromeres$end, "\\D", ""))
# Chromosome names standartization
chrom.names.standard <- sapply(chrlen$chrom, lookup, table = chrom.names.lookuptable)
if(is.list(chrom.names.standard)){
  stop("Ambiguous chromosome names, check the chromosome length data")
  write("Ambiguous chromosome names, check the chromosome length data", file = report.file, append = TRUE)
}else{
  chrlen$chrom <- chrom.names.standard
}
rm(chrom.names.standard)
# The same for centromeres
chrom.names.standard <- sapply(centromeres$chrom, lookup, table = chrom.names.lookuptable)
if(is.list(chrom.names.standard)){
  stop("Ambiguous chromosome names, check the chromosome length data")
  write("Ambiguous chromosome names, check the chromosome length data", file = report.file, append = TRUE)
}else{
  centromeres$chrom <- chrom.names.standard
}
rm(chrom.names.standard)

## For centromeres file
# In case this file contains centromeres of different types
if("type" %in% colnames(centromeres)){
  centromeres <- filter(centromeres, type == 'centromere')
}
if(nrow(centromeres) > 24){
  centromeres <- centromeres %>% group_by(chrom) %>% summarise(start = min(start),
                                                                   end = max(end))
}

centromeres$center <- (centromeres$start + centromeres$end)/2
# For plotting I won't need start and end coords, just the center coordinate
centromeres <- select(centromeres, chrom, center)

# Merging chrlen and centromeres tables to get rid of redundant tables
chrlen <- left_join(chrlen, centromeres, by = "chrom", copy = FALSE)
rm(centromeres)

# For sake of plotting in downstream code:
chrlen$relative.length <- chrlen$length/max(chrlen$length)



# -----------------------------Plot parameters specifying ----------------------------------------------
## Input Parameters
plot.height <- 1500
plot.width <- 1000
x.margin <- 50
y.margin <- 10
# Chromosomes rectangles thickness
thick <- 0
# Triangles parameters
tw <- 6
th <- 10

# Plot characteristics
#plot.params <- list(height = 1500, width = 1000, x.margin = 50,
#                    y.margin = 10, thick = 0, trianglewidth = 6,
#                    triangleheight = 10)

## Dependent Parameters
x.inner <- plot.width - x.margin
y.inner <- plot.height - y.margin
inner.width <- plot.width - 2*x.margin
inner.height <- plot.height - 2*y.margin
N <- nrow(chrlen)
spacing <- inner.height/(N-1)
chrlen$element.i <- seq(N)

#--------------------------------- Plotting coordinates calculation ----------------------------------------
## Chromosomes coordinates
chrlen$AX <- x.margin
chrlen$AY <- y.inner - spacing*(chrlen$element.i - 1) - thick/2
chrlen$BX <- x.margin + inner.width*chrlen$relative.length
chrlen$BY <- y.inner - spacing*(chrlen$element.i - 1) + thick/2

## Centromeres coordinates
chrlen$center.relative <- chrlen$center/max(chrlen$length)
chrlen$center.Y <- chrlen$AY
chrlen$center.X <- x.margin + inner.width*chrlen$center.relative
chrlen$length.mbp <- chrlen$length / 1000000

## The mitochondrial DNA is not plotted
chrlen <- filter(chrlen, chrom != "chrM")

# length.mbp will be used as the text labels on the plot
chrlen$length.mbp <- round(chrlen$length.mbp, 0)
chrlen$length.mbp <- as.character(chrlen$length.mbp)
chrlen$length.mbp <- paste(chrlen$length.mbp, 'Mbp')
color_of_lines <- 'black'

#---------------------------------- Read mapping import -----------------------------------------------------
#jsut calling it now lam table (this mapping import was done above )
lam.table <- bed.files


# -----------------------Mappings plotting -------------------------------------------------------
lam.table <- left_join(lam.table, select(chrlen, chrom, AY, BY), by = 'chrom')
max_chr_len <- max(chrlen$length)
lam.table %>% mutate(rel.start = start/max_chr_len) %>%
  mutate(rel.end = end/max_chr_len) -> lam.table

lam.table %>% mutate(AX = x.margin + inner.width*rel.start) %>%
  mutate(BX = x.margin + inner.width*rel.end) -> lam.table

binwidth = 2
bins <- seq(0, inner.width + x.margin, binwidth)
if (bins[length(bins)] != inner.width + x.margin){
  bins <- c(bins, inner.width + x.margin)
}
lam.table$binX <- .bincode(lam.table$AX, bins)

# Calculating the coordinates for the offtargets
offtargets <- offtargets %>%
  filter(chrom != "chrM") %>%
  mutate(center = (start+end)/2,
    AX = x.margin + inner.width*center/max(chrlen$length),
         BX = AX-2.5, CX = AX+2.5) %>%                          # change the numbers here to adjust the widht of the triangles
  left_join(select(chrlen, chrom, AY), by = c("chrom")) %>%
  mutate(AY = AY, BY = AY - 15, CY = BY)                  # Change this number to adjust the height of the triangles
offtargets <- offtargets %>%
  unite(A, AX,AY, sep = ";") %>%
  unite(B, BX,BY, sep = ";") %>%
  unite(C, CX,CY, sep = ";")
offtargets <- offtargets %>%
  gather(key = "dot", value = "coord", c("A", "B", "C"))
offtargets <- offtargets %>%
  separate(coord, c("dotX", "dotY"), sep = ";", remove = TRUE) %>%
  mutate(dotX = as.numeric(dotX),
         dotY = as.numeric(dotY))

#lam.table$sample <- "bpA"
#lam.table <- lam.table %>% filter(read.dir == "rev")

#how to make things that have more than 3 reads

for(i in unique(lam.table$filename)){
  lam.table.temp <- filter(lam.table, filename == i & quality >= 20)

  lam.table.temp %>% group_by(chrom, binX) %>% summarise(AY = nth(AY, 1), reads.count = n()) -> lam

  lam$reads.count.log <- log(lam$reads.count + 1, 7)
  lam$AX <- lam$binX * binwidth
  barheight <- 5
  lam$BY <- lam$AY + lam$reads.count.log * barheight
  #lam$BY <- lam$AY + lam$reads.count/2000
  lam$BX <- lam$AX + binwidth


#this contains the information on read count
  lamplot <- ggplot(data = chrlen)+
    coord_cartesian(xlim = c(0, plot.width), ylim = c(0, plot.height))+
    geom_segment(aes(x = AX, y = AY, xend = BX, yend = BY), alpha = 0.2, color = "black", size = 1)+
    geom_segment(data = lam.table.temp,
                 aes(x = AX, y = AY, xend = BX, yend = BY),
                 alpha = 0.25, size = 1, color = 'red')+
    geom_point(aes(x = center.X, y = center.Y), size = 2, shape = 21, fill = 'deeppink2')+
    geom_text(aes(x = AX-35, y = AY+3, label = chrom), colour = "black")+
    geom_text(aes(x = BX + 50, y = AY+3, label = length.mbp), size = 3, colour = "black")+
    # the read.counts here is the minimum bin size!!!!
    geom_segment(data = lam %>% filter(reads.count > 4), aes(x = AX-1, y = AY+2.5, xend = AX-1, yend = BY+2.5),
                 colour = "blue", size = 1)+
    # geom_rect(data = offtargets,
    #           aes(xmin = AX, ymin = AY, xmax = BX, ymax = BY),
    #           fill = "red")+
    
    #this is the color of the triangles for off target locations
    geom_polygon(data = offtargets,
                 aes(x = dotX, y = dotY, group = locusDesc), fill = 'red')+
    theme_void()


  ggsave(paste('htgts_', i,
               format(Sys.time(), "%d%m%y_%H%M"),
               '.pdf', sep = ''),
         width = 15, height = 10,
         device = "pdf",
         plot = lamplot)
  rm(lam.table.temp)
}

lam.table %>% filter(filename == unique(lam.table$filename)[2] & quality >= 25) %>% nrow()


# Plotting the minicircle mappings ----------------------------------------
# Specify the minicircle length here
mc.length <- 1232
# Binning the minicircle
binwidth <- 5 # a width of the single bin in bp
bins <- seq(0, mc.length, binwidth)
if (bins[length(bins)] != mc.length){
  bins <- c(bins, mc.length)
}
for (i in unique(mc$filename)) {
  mc.tmp <- mc %>% dplyr::filter(filename == i)
  coverage <- c()
  for (j in seq(length(bins)-1)) {
    coverage <- c(coverage, count.overlap(mc.tmp$start, mc.tmp$end, bins[j], bins[j+1]))
  }
  covdf <- data.frame(bp = bins[-length(bins)], cov = coverage)

  ggplot(covdf)+
    geom_bar(aes(x = bp, y = cov), stat = "identity")+
    theme_classic()+
    scale_x_continuous(breaks = seq(0, mc.length, 50))

  ggsave(paste('mc_', i,
               format(Sys.time(), "%d%m%y_%H%M"),
               '.pdf', sep = ''),
         width = 10, height = 3,
         device = "pdf")
  rm(mc.tmp, covdf)
}
