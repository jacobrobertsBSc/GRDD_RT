# Flow of information
## Footprint finder - flow of information

Import 6mA pileup 
        |>
Import fire_peaks 
        |>
Footprint finder function 
        |>
For each fire peak  - get all methylation data overlapping that peak 
        |>
Look for consecutive runs (8) of adenines thta have a 
methylation ratio below 0.25 & span more than 20bp
        |>
Treat the low methylation runs as candidate TF footprints
        |>
filter them for those that have flanks with methylation ratios over 0.5
        |>
Export to TSV or BED

## Motif Enrichment 

Obtain JASPAR motifs 
        |>
Get intervals for the footprints from BED import as Granges
        |>
Perform the motifenrichr
        |>
Find which Motifs are most enriched within our large footprints



