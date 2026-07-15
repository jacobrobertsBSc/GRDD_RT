library("TFBSTools")
library("JASPAR2022")
library(BiocFileCache)
library(dbplyr)
library("JASPAR2022")
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)

#Getting the JASPAR motifs with code from Michaelis 

# setting options for JASPAR motif db
opts <- list()
opts[["species"]] <- 'Homo sapiens'
opts[["all_versions"]] <- FALSE
opts[["matrixtype"]] <- 'PWM'

# extract motifs corresponding to human transcription factors
pwm_library <- getMatrixSet(JASPAR2022,opts)

pwm_library %>% head
pwm_library[[1]]

#Using lapply() to go through each PWM in pwm library and take out the ID and name 
pwm_library_list <- lapply(pwm_library, function(x){
  data.frame(ID = ID(x), name = name(x))
})

# combine the list into one data frame
pwm_library_df <- dplyr::bind_rows(pwm_library_list)
pwm_library_df %>% nrow
pwm_library_df %>% head

Footprints_tsv <- read.table(file = '/Users/j.roberts12/Desktop/footprint_finder/110326_h1_S4_hia5_hg38_pileup_large_footprints_wholegen.tsv', sep = '\t', header = TRUE)

Footprints_larger_than_20bp <- Footprints_tsv[Footprints_tsv$g_span > 20, ]

Footprints_larger_than_20bp_gr <- GRanges(seqnames = Footprints_larger_than_20bp$peak_chromosome, ranges = IRanges(start = Footprints_larger_than_20bp$g_start, end = Footprints_larger_than_20bp$g_end))

#https://greenleaflab.github.io/motifmatchr/articles/motifmatchr.html
motif_ix <- matchMotifs(pwm_library, Footprints_larger_than_20bp_gr, genome = "BSgenome.Hsapiens.UCSC.hg38", out = "scores")

scores <- motifScores(motif_ix)

#Needed to use colMeans from the Matrix library as the matrix is a sparseMatrix - where the majority of elements are zero it onky stores non-zero values and coordinates - though it can still print them to the screen
mean_scores <- Matrix::colMeans(scores)

#Dataframe for the IDs and mean enrichment scores for each ID
results <- data.frame(ID = names(mean_scores), mean_score = as.numeric(mean_scores))

results <- merge(results, pwm_library_df, by = "ID")

results <- results[order(results$mean_score, decreasing = TRUE), ]


 matches <- motifMatches(motif_ix)

#finding out how many rows (footprints) have a true (1) value in matches - sum only counts 1 - true - not 0 false 

  sum(matches[, "MA0139.1"])

#because CTCF has multiple JASPAR motif IDs - i made a list of the matches with each ID to then sum them 

ctfc_any <- matches[, "MA0139.1"], matches[, "MA1102.2"], matches[, "MA1102.1"], matches[, "MA1102.2"], matches[, "MA1102.3"], matches[, "MA1929.1"], matches[, "MA1929.2"], matches[, "MA1930.1"], matches[, "MA1930.2"]

#but this threw errors as not all of these motifs were in matches so i made a list of just the IDs and checked if they were in matches 

ctcf_ids <- c(
   "MA0139.1",
   "MA1102.1",
   "MA1102.2",
   "MA1102.3",
   "MA1929.1",
   "MA1929.2",
   "MA1930.1",
   "MA1930.2"
)

ctcf_ids %in% colnames(matches)

ctcf_ids_present <- ctcf_ids[ctcf_ids %in% colnames(matches)]

matches_all_ctcf <- matches[, ctcf_ids_present]

#https://stackoverflow.com/questions/51025979/equivalent-of-rowsum-function-for-matrix-class-dgcmatrix

#For every row - is it true that in matched_all_ctcf columns - there are more than 0 matches 

num_ctcf_footprints <- Matrix::rowSums(matches_all_ctcf) > 0

#Gives a table shows number of trues and falses
table(num_ctcf_footprints)

#Taking out all of the footprint rows that are not true in num_ctcf_footprints 
 
non_ctcf_matches <- matches[!num_ctcf_footprints, ]

#Using colSums from the Matrix library to access the sparse matrix and sum the columns for true and false for each footprint
 
non_ctcf_counts <- Matrix::colSums(non_ctcf_matches)

#creating a data frame from the named vector non_ctcf_counts
 
results_non_ctcf <- data.frame(ID = names(non_ctcf_counts), footprints = non_ctcf_counts)

results_non_ctcf <- merge(results_non_ctcf, pwm_library_df, by = "ID")

#ordering the merged data frame by decreasing number of footprints that the TF motif is in
results_non_ctcf <- results_non_ctcf[order(results_non_ctcf$footprints, decreasing = TRUE),]

head(results_non_ctcf, 40)

pdf("/Users/j.roberts12/Desktop/footprint_finder/non_ctcf_footprint_counts_stage_4_hg38_h1.pdf", width = 16, height = 8)

plot <- ggplot(head(results_non_ctcf,30), aes(x = reorder(name, -footprints), y = footprints)) +
	geom_col(stat = "identity", fill = "red") +  
	labs(x = "Transcription factor", y = "Number of footprints") +
	scale_y_continuous(breaks = seq (0,45, by = 5)) + 
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)

dev.off()
