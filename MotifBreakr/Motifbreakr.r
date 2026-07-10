---
title: "R Notebook"
output: html_notebook
---



```{r}
#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("motifbreakR")

```
```{r}
library(motifbreakR)
library(BSgenome)
```
```{r}
reduced<- import("/Users/j.roberts12/Desktop/Differential_methylation_modkit/hglft_genome_hg38_heatmap_regions_foxa2_control.bed", format = "BED")

reduced_df <- as.data.frame(reduced)

reduced_df$region <- paste0(
  reduced_df$seqnames, ":",
  reduced_df$start, "-",
  reduced_df$end
)

regions_list <- as.list(reduced_df$region)

regions_list

#defines a function with the argument being names region_string
#splits the first string element of the list [[1]] at : and - and stores is in p -- 
#the [[ ]] is important as it brings out the actual character string instead of just a subsection of a list which is still a list and not the character vector 
#https://adv-r.hadley.nz/subsetting.html
#The function then uses each index of that character vector and assigns it to the seqnames etc to generate a Granges range 
#https://www.r-bloggers.com/2024/04/exploring-strsplit-with-multiple-delimiters-in-r/#google_vignette -- for p <-strsplit

regions_function <- function(region_string){
  p <- strsplit(region_string, "[:-]")[[1]]
  
  GRanges(seqnames = p[1], 
          ranges = IRanges(start = as.integer(p[2]),
                           end = as.integer(p[3])
          )
  )
}

#applying the function above to the list of regions from regions list

regions <- lapply(X = regions_list, FUN = regions_function )

#because lapply creates a list of regions -- this function basically makes sure that every element is a granges (as assurance) within the list - if they aren't then it will error
#with this it unlists all of the grnages objects and concatenates them into grnages --> much like the do.callc() function but forcing the granges method 

gr_regions <- unlist(
  GRangesList(regions),
  use.names = FALSE
)
```

```{r}
significant_hits_0.05 <- read_tsv("/Users/j.roberts12/Desktop/Differential_methylation_modkit/endoc60_hia5_hg38_original_190326_dmr_filtered_sorted_chr20_22334362_22554362_control_region_significant_hits_0.05_onlyinheatmapregions.bed", col_names = FALSE)



colnames(significant_hits_0.05)[1:3]  <- c("chr", "start", "end")
colnames(significant_hits_0.05)[5]  <- "score"
colnames(significant_hits_0.05)[15] <- "MAP-based p-value"
colnames(significant_hits_0.05)[16] <- "effect_size"

print(significant_hits_0.05)

gr <- makeGRangesFromDataFrame(significant_hits_0.05, 
                               keep.extra.columns = TRUE, 
                               ignore.strand = FALSE)
gr

significant_hits_in_heatmap_regions <- subsetByOverlaps(gr, gr_regions)

bed_gr <- makeGRangesFromDataFrame(significant_hits_in_heatmap_regions, 
                               keep.extra.columns = TRUE, 
                               ignore.strand = FALSE)

length(significant_hits_in_heatmap_regions)

export.bed(bed_gr, con = "/Users/j.roberts12/Desktop/Differential_methylation_modkit/endoc60_adenine_filtered_sorted_chr20_22315000_22535000_hg19_significant_calls_0.05_onlyinheatmapregions.bed")
```
```{r}
bed_df <-as.data.frame(bed_gr)

colnames(bed_df)[1] <- "chr"

bed_df

motifbreakr_ready <- bed_df %>%
  mutate(
    name = paste(chr, start, "A", "N", sep = ":"), 
    score = 0,                                     
  ) %>%
  select(chr, start, end, name, score, X6)

write_tsv(motifbreakr_ready, 
          "/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_input_endoc60_adenine_filtered_sorted_chr20_22315000_22535000_hg19_significant_calls_0.05_onlyinheatmapregions.bed", 
          col_names = FALSE)

motifbreakr_ready
```

```{r}
snps_file <- "/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_input_endoc60_hia5_hg38_original_190326_dmr_filtered_sorted_chr20_22334362_22554362_control_region_significant_hits_0.05_onlyinheatmapregions.bed"

input_data <- read.delim(snps_file, header = FALSE)

print(input_data)
```
```{r}
genomes <- available.genomes()
length(genomes)

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
```
```{r}
my_6ma_sites <- snps.from.file(file = snps_file, search.genome = BSgenome.Hsapiens.UCSC.hg38, format = "bed")

my_6ma_sites

as.data.frame(my_6ma_sites)
```
```{r}
library(MotifDb)

MotifDb
```
```{r}
data(motifbreakR_motif)

motifbreakR_motif
```
```{r}
mcols(my_6ma_sites)$REF <- getSeq(BSgenome.Hsapiens.UCSC.hg38, my_6ma_sites)

mcols(my_6ma_sites)$REF

my_6ma_sites

all_motif_results <- motifbreakR(snpList = my_6ma_sites, 
                                 filterp = TRUE,
                                 pwmList = motifbreakR_motif, 
                                 threshold = 1e-4,
                                 method = "ic", 
                                 bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                 BPPARAM = BiocParallel::SerialParam())

all_motif_results

head(all_motif_results, 70)

motif_df <- as.data.frame(all_motif_results, row.names = NULL)

motif_df

write_tsv(
  motif_df,
  "/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_input_endoc60_adenine_filtered_sorted_chr20_22315000_22535000_hg19_significant_calls_0.05_onlyinheatmapregions_all_motif_results.tsv"
)

tsv <- "/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_input_endoc60_adenine_filtered_sorted_chr20_22315000_22535000_hg19_significant_calls_0.05_onlyinheatmapregions_all_motif_results.tsv"

input_data <- read.delim(tsv)

strong_effect_hits <- input_data %>% filter(effect == "strong")

write_tsv(
  strong_effect_hits,
  "/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_input_endoc60_adenine_filtered_sorted_chr20_22315000_22535000_hg19_significant_calls_0.05_onlyinheatmapregions_all_motif_results_strong_effects.tsv"
)

strong_effect_hits


```
```{r}
top_hits_alleleeffectsize <- all_motif_results[head(order(abs(all_motif_results$alleleEffectSize), decreasing = TRUE), 20),]


top_hits_alleleeffectsize_p_value <- calculatePvalue(top_hits_alleleeffectsize)

top_hits_alleleeffectsize_p_value

top_hits_alleleeffectsize_p_value_df <- as.data.frame(top_hits_alleleeffectsize_p_value, row.names = NULL)

write_tsv(
  top_hits_alleleeffectsize_p_value_df,
  "/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_input_endoc60_hia5_hg38_original_190326_dmr_filtered_sorted_chr20_22334362_22554362_control_region_significant_hits_0.05_onlyinheatmapregions_all_motif_results_strong_effects_top_hits_alleleeffectsize.tsv"
)

```


```{r}
as.data.frame(my_6ma_sites)
```

```{r}

#hg19_endoc60_ecog11_ftnucs_hplotagged regions for poster

site_22471691 <- all_motif_results[all_motif_results$SNP_id == "chr20:22471691:A:N"]

# Sites for Hia5 endoc data 
#site_22402823 <- all_motif_results[all_motif_results$SNP_id == "chr20:22402823:A:N"]
#site_22470506 <- all_motif_results[all_motif_results$SNP_id == "chr20:22470506:A:N"]
#site_22484337 <- all_motif_results[all_motif_results$SNP_id == "chr20:22484337:A:N"]
#site_22484426 <- all_motif_results[all_motif_results$SNP_id == "chr20:22484426:A:N"]
#site_22503819 <- all_motif_results[all_motif_results$SNP_id == "chr20:22503819:A:N"]
#site_22533981 <- all_motif_results[all_motif_results$SNP_id == "chr20:22533981:A:N"]
#site_22418853 <- all_motif_results[all_motif_results$SNP_id == "chr20:22418853:A:N"]
#site_22418853

#site_22448990 <- all_motif_results[all_motif_results$SNP_id == "chr20:22448990:A:N"]
#site_22448990

#site_22449256 <- all_motif_results[all_motif_results$SNP_id == "chr20:22449256:A:N"]
#site_22449256

#site_22449844 <- all_motif_results[all_motif_results$SNP_id == "chr20:22449844:A:N"]
#site_22449844

#site_22460280 <- all_motif_results[all_motif_results$SNP_id == "chr20:22460280:A:N"]
#site_22460280

#site_22463532 <- all_motif_results[all_motif_results$SNP_id == "chr20:22463532:A:N"]
#site_22463532

#site_22463534 <- all_motif_results[all_motif_results$SNP_id == "chr20:22463534:A:N"]
#site_22463534
#site_22491048 <- all_motif_results[all_motif_results$SNP_id == "chr20:22491048:A:N"]
#site_22491048

#site_22491055 <- all_motif_results[all_motif_results$SNP_id == "chr20:22491055:A:N"]

#FOXA2 s2vss4 h1 
site_22484353 <- all_motif_results[all_motif_results$SNP_id == "chr20:22484353:A:N"]
site_22484353

site_22484354 <- all_motif_results[all_motif_results$SNP_id == "chr20:22484354:A:N"]
site_22484354

site_22484356 <- all_motif_results[all_motif_results$SNP_id == "chr20:22484356:A:N"]
site_22484356

site_22484359<- all_motif_results[all_motif_results$SNP_id == "chr20:22484359:A:N"]
site_22484359

site_22484360<- all_motif_results[all_motif_results$SNP_id == "chr20:22484360:A:N"]
site_22484360

site_22484362<- all_motif_results[all_motif_results$SNP_id == "chr20:22484362:A:N"]
site_22484362

site_22484367<- all_motif_results[all_motif_results$SNP_id == "chr20:22484367:A:N"]
site_22484367

site_22484369<- all_motif_results[all_motif_results$SNP_id == "chr20:22484369:A:N"]
site_22484369


#site_22491055
#site_22402823
#site_22470506
#site_22484337
#site_22484426
#site_22503819
#site_22533981
```



```{r}

#hg19_endoc60_ecog11_ftnucs_hplotagged regions for poster
site_22471691_pvalue <- calculatePvalue(site_22471691, granularity = 1e-6)

#site_22402823_pvalue <- calculatePvalue(site_22402823, granularity = 1e-6)
#site_22484337_pvalue <- calculatePvalue(site_22484337, granularity = 1e-6)
#site_22503819_pvalue <- calculatePvalue(site_22503819, granularity = 1e-6)
#site_22533981_pvalue <- calculatePvalue(site_22533981, granularity = 1e-6)
site_22418853_pvalue <- calculatePvalue(site_22418853, granularity = 1e-6)
site_22418853_pvalue

site_22448990_pvalue <- calculatePvalue(site_22448990, granularity = 1e-6)
site_22448990_pvalue

site_22449256_pvalue <- calculatePvalue(site_22449256, granularity = 1e-6)
site_22449256_pvalue

site_22460280_pvalue <- calculatePvalue(site_22460280, granularity = 1e-6)
site_22460280_pvalue

#all_motif_results_p_value <- calculatePvalue(all_motif_results, granularity = 1e-6)
all_motif_results_p_value

#site_22491048_pvalue <- calculatePvalue(site_22491048, granularity = 1e-6)
site_22491048_pvalue

#site_22491055_pvalue <- calculatePvalue(site_22491055, granularity = 1e-6)
site_22491055_pvalue

#site_22402823_pvalue
#site_22484337_pvalue
#site_22503819_pvalue
#site_22533981_pvalue

```


```{r}
install.packages("Cairo")

library(Cairo)
```
```{r}

FOXA2_binding_1 <- c("chr20:22484353:A:N", "chr20:22484354:A:N", "chr20:22484356:A:N", "chr20:22484359:A:N","chr20:22484360:A:N", "chr20:22484362:A:N", "chr20:22484367:A:N", "chr20:22484369:A:N")

pdf("Region_70_22483998_FOXA2_Motif_Plots.pdf", width = 12, height = 30)

for (dmr in FOXA2_binding_1) {
  print(plotMB(results = all_motif_results, rsid = dmr, effect = "strong"))
} 

dev.off()
```

```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22503819_endoc60_hia5_hg38_.pdf", width = 12, height = 30)

plot_22503819 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22503819:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22503819

dev.off()
```
```{r}

pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22484337_endoc60_hia5_hg38_.pdf", width = 12, height = 30)

plot_22484337 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22484337:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22484337

dev.off()
```

```{r}

pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22470697_endoc60_hia5_hg38_.pdf", width = 12, height = 30)

plot_22470697 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22470697:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22470697

dev.off()

```
```{r}


pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22471018_endoc60_hia5_hg38_.pdf", width = 12, height = 14)

plot_22471018 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22471018:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22471018

dev.off()
```
```{r}


pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22484426_endoc60_hia5_hg38_.pdf", width = 12, height = 12)

plot_22484426 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22484426:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22484426

dev.off()
```


```{r}

pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22471265_endoc60_ecogii_hg19_.pdf", width = 12, height = 30)

plot_22471265 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22471265:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22471265

dev.off()
```

```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22356946_endoc60_ecogii_hg19_FOXA.pdf", width = 12, height = 30)

plot_22356946 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22356946:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22356946

dev.off()
```

```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22433232_endoc60_ecogii_hg19_IRX.pdf", width = 12, height = 30)

plot_22433232 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22433232:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22433232

dev.off()
```


```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22429710_endoc60_ecogii_hg19_hdac2.pdf", width = 12, height = 30)

plot_22429710 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22429710:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22429710

dev.off()
```


```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22429701_endoc60_ecogii_hg19_MAFB.pdf", width = 12, height = 30)

plot_22429701 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22429701:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22429701

dev.off()
```

```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22471258_endoc60_ecogii_hg19.pdf", width = 12, height = 30)

plot_22471258 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22471258:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22471691

dev.off()
```
```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/motifbreakr_plot_22410254_endoc60_ecogii_hg19.pdf", width = 12, height = 30)

plot_22410254 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22410254:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22410254

dev.off()
```


```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/plot_22491055_motif_breakr_hia5_h1_s2_vs_s4_hg38.pdf", width = 12, height = 30)

plot_22491055 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22491055:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22491055

dev.off()
```



```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/plot_22448990_motif_breakr_hia5_h1_hg38.pdf", width = 12, height = 20)

plot_22448990 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22448990:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22448990

dev.off()
```
```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/plot_22418853_motif_breakr_hia5_h1_hg38.pdf", width = 12, height = 20)

plot_22418853 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22418853:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22418853

dev.off()
```

```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/plot_22449256_motif_breakr_hia5_h1_hg38.pdf", width = 12, height = 20)

plot_22449256 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22449256:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22449256

dev.off()
```


```{r}
plot_22402823 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22402823:A:N", 
       effect = "strong", 
       altAllele = "N")
```
```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/plot_22484337_motif_breakr_hia5_endoc60_hg38.pdf", width = 12, height = 20)

plot_22484337 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22484337:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22484337

dev.off()
```
```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/plot_22503819_motif_breakr_hia5_endoc60_hg38.pdf", width = 12, height = 20)

plot_22503819 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22503819:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22503819

dev.off()
```
```{r}
pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/plot_22533981_motif_breakr_hia5_endoc60_hg38.pdf", width = 12, height = 20)

plot_22533981 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22533981:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22533981

dev.off()
```

```{r}

pdf("/Users/j.roberts12/Desktop/Differential_methylation_modkit/plot_22402823_motif_breakr_hia5_endoc60_hg38.pdf", width = 12, height = 20)

plot_22402823 <- plotMB(results = all_motif_results, 
       rsid = "chr20:22402823:A:N", 
       effect = "strong", 
       altAllele = "N")

plot_22402823

dev.off()
```

