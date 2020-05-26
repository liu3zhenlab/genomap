### generate segments
segres <- genosegF2(geno = genodata, genocols = c(5,6), chrname = "CHR", posname = "POS",
                    output.common = "test", chromosomes = c(1, 10), min.seg.size = 100000,
                    allele1.name = "A", allele2.name = "B", hetero.name = "H", missing.name = 0)

### genotyping scoring
seg2score(seg.input = "test.seg.txt", segscore.output = "test.segmarker.txt",
          missing.data.code = 0, binsize = 1000000)
