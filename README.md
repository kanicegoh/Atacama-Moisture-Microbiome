# Atacama Moisture Microbiome

Bioinformatic and statistical analyses of microbial diversity across soil moisture gradients in the Atacama Desert.

## Workflow

1. Sequence preprocessing
2. Primer removal using Cutadapt
3. ASV inference using DADA2
4. Taxonomic assignment using SILVA v138.2
5. Rarefaction
6. Alpha diversity analyses
7. Statistical modelling
8. Figure generation

## Software

- R
- DADA2
- phyloseq
- ggplot2
- ggeffects

## Notes

The workflow implemented in this repository was adapted from the standard DADA2 pipeline and modified to accommodate multiple published datasets with different primer sets, sequencing configurations, read lengths, and quality profiles.

## References

This workflow was adapted from the official DADA2 pipeline tutorial:

Callahan BJ. DADA2 Pipeline Tutorial.
https://benjjneb.github.io/dada2/tutorial.html

The DADA2 software package is described in:

Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP.
DADA2: High-resolution sample inference from Illumina amplicon data.
Nature Methods. 2016;13:581–583.
https://doi.org/10.1038/nmeth.3869

Taxonomic assignment was performed using the SILVA v138.2 reference database.

## Author

Kanice Goh  
National University of Singapore
