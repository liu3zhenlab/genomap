# package genoseg

This package was dedicated to determine chromosomal regions, each of which harbor multiple markers (SNPs) most of which show the same genotypes. Segments are used to represent chromosomal regions. The algorithm used in DNAcopy (https://bioconductor.org/packages/release/bioc/html/DNAcopy.html) is used to determine segments.

## INSTALL
Download the package, and run the installation.
```{r}
genose.source <- ""~/RProjects/genoseg_0.0.1.tar.gz" # path is subject to change
install.packages(genoseg.source, repos = NULL, type = "source")

library(genoseg)
```

## Input data
Basically, genotyping data are required. Genotype data are organized in a tab-delimited format. For each entry, chromosome and chromsomal coordinate are provided. Genotypes of multiple individuals can be included in the same data file. Here are the first six line of the example data.

```{r, comment = "", echo = F}
head(genodata)
```

## Brief introduction of module usages


