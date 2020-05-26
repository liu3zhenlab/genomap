grep -e "^1" -e CHR /data1/home/liu3zhen/A188-GeneticMapping/F2GBS/5-BAF2calli_merge/4o-snps.genotype.txt | cut -f 1-6 | sed 's/XT.II_/S/g'> genodata.txt
### read genodata.txt to R and convert to .rda data format by using save(xxx, file = "genodata.rda")
#genodata <- read.delim("../raw/genodata.txt")
#save(genodata, file = "genodata.rda")
