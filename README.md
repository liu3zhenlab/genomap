## genomap

Multiple modules have been developed for genetic and genomic studies in this package.
1. segment genotyping
This part is dedicated to determine chromosomal regions (segments) with approximately the same genotype supported by multipe markers. The algorithm from DNAcopy (https://bioconductor.org/packages/release/bioc/html/DNAcopy.html) is used to determine segments.

2. prediction of genetic positions
The module was designed to predict genetic positions of any physical position based on an input map.

### INSTALLATION
```
# install from GitHub:
# install.packages("devtools")
devtools::install_github("liu3zhenlab/genomap")
```
