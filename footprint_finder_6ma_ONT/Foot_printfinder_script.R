library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

#Loading in Pileup data 

Pileup_data <- "/penrose/projects/J.Roberts12/haplotype_specific_6ma/differentiation_foxa2_diff_meth/110326_h1_S4_hia5_hg38_pileup.bed.gz"

extraCols_pileup_data <- c(
  Nvalid_cov       = "integer",
  percent_modified = "numeric",
  Nmod             = "integer",
  Ncanonical       = "integer",
  Nother_mod       = "integer",
  Ndelete          = "integer",
  Nfail            = "integer",
  Ndiff            = "integer",
  Nnocall          = "integer"
  )

#importing the pileup data as a granges object

gr_pileup_data <- import(Pileup_data, format = "BED", extraCols = extraCols_pileup_data)

#obtaining the fraction modified by dividing the percent by 100

gr_pileup_data$Methylation_ratio <- gr_pileup_data$percent_modified / 100

#Loading in the FIRE peaks and making them into a GRanges - not sure if better to do it this way or with import and defining the extra columns

fire <- read.delim(gzfile("/penrose/projects/J.Roberts12/data/Stage_4_hia5_hg38_FDR_peaks/FDR-FIRE-peaks.bed.gz"), header = TRUE)

fire_gr <- GRanges(seqnames = fire$X.chrom, ranges = IRanges(start = fire$peak_start + 1, end = fire$peak_end))


#------------- Parameters for finding footprint --------------------

LOW_THRESH   <- 0.25   # below this = footprint candidate
HIGH_THRESH  <- 0.50   # flanking regions must reach this
MIN_FP_BASES <- 8     # minimum number of consecutive A-bases
MIN_FP_BP    <- 10     # minimum genomic span in bp

#------------- Defining the function to search for footprints ---------------


find_footprints <- function(fire_peak, peak_chromosome, peak_start, peak_end) { 
  
  if (length(fire_peak) < MIN_FP_BASES) return(NULL) #if the peak length is less than the minimum number of adenines that we're looking for then dont check

  below <- fire_peak$Methylation_ratio < LOW_THRESH # Creating a logical vector for below and above the threshold for a footprint
  
  r1 <- rle(below) # Compresses consecutive TRUE/FALSE into runs - 2,1,13,4,5,6 T,F,T,F,T,F
  
  positions <- start(fire_peak) #storing the positions within the fire peak so that the indexes above can be used to map the footprints found
  
  end_idx <- cumsum(r1$lengths) #Calculating the index of the end of each run from the rle 

  start_idx <- c(1, end_idx[-length(end_idx)] + 1) # Calculating the start of each run from the end of the previous run - this first run starts at 1 and the second one starts at the length of the end idx +1 

  #Creating a data.frame of all the runs - whether they are a footprint or not - finding the length of the run - and then getting the positions of the runs with the indexes created above 
  
  runs <- data.frame(is_footprint = r1$values, n_bases = r1$lengths, g_start = positions[start_idx], g_end = positions[end_idx], g_span = positions[end_idx] - positions[start_idx])
  
  
  #Filtering the runs for footprints that have a span of MIN_FP_BP and MIN_FP_bases  
  candidates <- runs[runs$is_footprint & runs$n_bases >= MIN_FP_BASES & runs$g_span  >= MIN_FP_BP, ]
    
  #Function to check the flanks of the footprint are high 
  check_flanks <- function(start, end, positions, values, flank_size = 5) {
  
  #Getting the values from the positions before and after the candidate footprint
  left_side <- values[positions < start]
  right_side <- values[positions > end]
  
  #Only taking the desired flank size of values 
  left_flank  <- tail(left_side, flank_size)
  right_flank <- head(right_side, flank_size)
  
  #Checking that the length of values above the flank threshold is higher than zero and have an average over the High threshold - same for left flank
  left_ok  <- length(left_flank) > 0 && mean(left_flank) >= HIGH_THRESH
  right_ok <- length(right_flank) > 0 && mean(right_flank) >= HIGH_THRESH
  
  return(left_ok & right_ok)
  }

#Applying check to the candidates 

  candidates$flank_ok <- sapply(seq_len(nrow(candidates)), function(i) {
  check_flanks(
    candidates$g_start[i],
    candidates$g_end[i],
    positions,
    fire_peak$Methylation_ratio
  )
  })
  
  #Keeping the information about the peak in which the candidate footprint was found
  candidates$peak_chromosome <- peak_chromosome
  candidates$peak_start <- peak_start
  candidates$peak_end   <- peak_end
  
  candidates

}


#Running the function to discover footprints - 

footprints <- vector("list", length(fire_gr))

#To find which pileup data falls within peaks - shows which pileup position belongs to which peak

hits <- findOverlaps(fire_gr, gr_pileup_data)

#Groups the pileup data indexes by peak indexes so now there is a list of peaks that contain the corresponding pileup data indexes - if we used subsetByoverlaps inside the loop it would likely take longer as it is 
#Doing a search of all the positions with every fire peak 

hits_by_peak <- split(subjectHits(hits), queryHits(hits))

for (i in seq_along(hits_by_peak)) {

    idx <- hits_by_peak[[i]] #For each peak it extracts the indexes of the Adenines within it 

    peak_methylation_data <- gr_pileup_data[idx] #using those indexes to get the data out of gr_pile_up

    peak_coord <- fire_gr[i]

    footprints[[i]] <- find_footprints(
        fire_peak = peak_methylation_data,
        peak_chromosome = as.character(seqnames(peak_coord)),
        peak_start = start(peak_coord),
        peak_end = end(peak_coord)
    )
}

#Footprints is now a list of elements containing information about footprints - inlcuding if it is NULL - this uses do.call to apply rbind to all of those that make it into the Filter list

all_footprints_df <- do.call(rbind, Filter(Negate(is.null), footprints))

confirmed <- all_footprints_df[all_footprints_df$flank_ok, ]

write_tsv(confirmed,"/penrose/projects/J.Roberts12/haplotype_specific_6ma/differentiation_foxa2_diff_meth/footprint_finding/110326_h1_S4_hia5_hg38_pileup_large_footprints_wholegen.tsv")



footprint_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, confirmed_gr))

writeXStringSet(footprint_sequences, filepath = "/penrose/projects/J.Roberts12/footprint_finder/hg38_Stage4_H1_largefootprints_fasta.fa")
