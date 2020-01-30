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
folder.with.sorted.bed <- "/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/Mishas_Demultiplexed_LAM/bed_files"
#output
setwd("/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/Mishas_Demultiplexed_LAM/bed_plots")

#################       Tn5  ###################
# input
#folder.with.sorted.bed <- "/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/P_Eric4_Tn5/bed_files"
# output
#setwd("/media/edanner/NewUbuntuSpace/Workspace/LinearAmp/Sequence2_191129_MN00157_0047_A000H2GWGF/P_Eric4_Tn5/bed_plots")



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
if(libs.installed == length(libs)){
  print("All libraries installed successfully")
}else{
  stop("Some libraries were not installed properly")
}
rm(libs, libs.installed)

# Functions ---------------------------------------------------------------
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



# File search -------------------------------------------------------------
chrlength.file <- search_files(folder = getwd(), file.or.pattern = "chrlength")

# chrlength.file <- read.table("/home/edanner/workspace/uditas/Pipeline2Plotting/GRCh38_chrlength.txt",
#                              header = TRUE, stringsAsFactors = FALSE)
## The program uses just one file with chromosome length, it chooses the most recent one
if(length(chrlength.file) > 1){
  file.info.table <- bind_cols(data.frame(filename = chrlength.file), file.info(chrlength.file))
  file.info.table <- arrange(file.info.table, desc(atime))
  chrlength.file <- file.info.table$filename[1]
  warning(paste("More than one file with chromosome length found, using the most recent one, ", chrlength.file, ".", sep = ""))
  write(paste("More than one file with chromosome length found, using the most recent one, ", chrlength.file, ".\n", sep = ""),
        file = report.file, append = TRUE)
  rm(file.info.table)
}else{
  write(paste("Chromosome length file used: ", chrlength.file, ".\n", sep = ""),
        file = report.file, append = TRUE)
}


centromere <- read.table("/home/edanner/workspace/uditas/Pipeline2Plotting/hg38_centromeres.txt",
                              header = TRUE, stringsAsFactors = FALSE)
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

bed.files <- search_files(folder = folder.with.sorted.bed, file.or.pattern = ".sorted.bed")
bed <- data.frame()
for(i in bed.files){
  filename <- paste0(folder.with.sorted.bed, "/", i)
  bed.temp <- read.table(filename, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  bed.temp$filename <- i
  bed <- bind_rows(bed, bed.temp)
}
  bed.files <- bed %>% select(V1:V12, filename)

# 
# bed.files <- read.table(paste0(folder.with.sorted.bed, "/N707_N505.sorted.bed"), stringsAsFactors = FALSE,
#                         fill = TRUE)
# 
colnames(bed.files) <- c("chrom", "start", "end", "read.id", "quality", "strand", "flag", "cigar", "asterisk", "zero1", "zero2", "seq", "filename")

bed.files <- bed.files %>% filter(!str_detect(chrom, "_") | chrom == "pE038_MC")
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




# Data import and check ---------------------------------------------------
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
chrlen <- read.table(file = as.character(chrlength.file), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
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

centromeres <- centromere
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


# Data preprocessing ------------------------------------------------------
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


# Plot parameters specifying ----------------------------------------------
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


# Plotting coordinates calculation ----------------------------------------
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



# Read mapping import -----------------------------------------------------
#jsut calling it now lam table
lam.table <- bed.files
# lam.table <- read.table(file = paste(folder.with.sorted.bed, "/", bed.files[1], sep = ""),
#                         sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#lam.table$sample <- unlist(str_split(bed.files[1], ".sorted.bed"))[1]
# if(length(bed.files)>1){
#   for(i in seq(2,length(bed.files))){
#     lam.table.temp <- read.table(file = paste(folder.with.sorted.bed, "/", i, sep = ""),
#                                  sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#     lam.table.temp$sample <- unlist(str_split(i, ".sorted.bed"))[1]
#     lam.table <- bind_rows(lam.table, lam.table.temp)
#   }
#   rm(i)
# }

# Filtering out non-used columns and standartization of the rest of the columns
# lam.table %>% select(c("V1", "V2", "V3", "V4", "sample")) -> lam.table

#colnames(lam.table) <- c("chrom", "start", "end", "read.id", "sample")
# Chromosome mappings can contain undefined chromosome scaffolds, filtering them out
# if(any(str_detect(unique(lam.table$chrom), "_"))){
#   warning("Mappings contain the undefined chromosome scaffolds alignments, filtering them out.")
#   write("Mappings contain the undefined chromosome scaffolds alignments, filtering them out.",
#         file = report.file, append = TRUE)
#   n.observations.before <- nrow(lam.table)
#   lam.table %>% filter(!str_detect(chrom, "_")) -> lam.table
#   n.observations.after <- nrow(lam.table)
#   n.obs.filtered <- n.observations.before - n.observations.after
#   the.message <- paste(n.obs.filtered, " reads (", round(n.obs.filtered*100/n.observations.before,1)," %)",
#                    " were filtered out on this step (out of ",
#                    n.observations.before, " reads, ", n.observations.after, " reads left).", sep = "")
#   message(the.message)
#   write(the.message, file = report.file, append = TRUE)
#   rm(the.message)
#   lam.table %>% filter(!str_detect(chrom, "chrM")) -> lam.table
#   n.observations.after.chrM <- nrow(lam.table)
#   the.message <- paste("Number of reads filtered out due to mapping to the chrM: ",
#                        n.observations.after - n.observations.after.chrM, " (",
#                        round((n.observations.after - n.observations.after.chrM)*100/n.observations.after, 1),
#                        " %).", sep = "")
#   message(the.message)
#   write(the.message, file = report.file, append = TRUE)
#   rm(the.message)
#   rm(n.observations.before, n.observations.after, n.observations.after.chrM, n.obs.filtered)
# }
# 

# Mappings plotting -------------------------------------------------------
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
    geom_segment(data = lam %>% filter(reads.count > 1), aes(x = AX-1, y = AY+2.5, xend = AX-1, yend = BY+2.5),
                 colour = "blue", size = 1)+
    theme_void()

  ggsave(paste('htgts_', i,
               format(Sys.time(), "%d%m%y_%H%M"),
               '.pdf', sep = ''),
         width = 15, height = 10,
         device = "pdf",
         plot = lamplot)
  rm(lam.table.temp)
}

  lam.table %>% filter(filename == unique(lam.table$filename)[1] & quality >= 20) %>% nrow()
lam.table %>% filter(filename == unique(lam.table$filename)[1]) %>% nrow()
lam.table$chrom %>% unique()

lamplot

ggplot(contigs) +
  geom_histogram(aes(x = length, fill = sample), position = position_identity(), alpha = 0.3, bins = 100) +
  scale_x_continuous(limits = c(0,5000)) + geom_vline(xintercept = 224) +
  theme_classic() + scale_y_continuous(expand = c(0,0)) +
  xlab("Length, bp") +
  ylab("Mappings count")

