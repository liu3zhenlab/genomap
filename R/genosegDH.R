#### DH segmentation to identify origins of segments
#### Sanzhen Liu
#### 8/25/2019

### requirments:
# 1. genotyping data: chr; pos; lines;
# 2. column NUMBER for lines
# 3. segmentation chromosome list
################################################################################################
# segment merging function (internal function)
################################################################################################
# function for merging the segments if neighboring genotypes are the same:
seg.merge <- function (seg.out) {
  if (nrow(seg.out) < 2) {
    seg.out.merge <- seg.out
  } else {
    # label the equality of the neighbor segments:
    seg.out$Equal2next <- 0 # initiate, 0 means seg.median of this row is not equal to the next
	for (i in 1:(nrow(seg.out) - 1)) {
      if (seg.out$seg.median[i] == seg.out$seg.median[i+1]) {
        seg.out$Equal2next[i] <- 1
      }
    }
    # judge and output:
    rownum <- 0
    seg.out.merge <- NULL
    repeat {
      rownum <- rownum + 1
      if (rownum <= nrow(seg.out)) {
        if (seg.out$Equal2next[rownum] == 1){ ### equal neighboring genotypes
          seg.out$loc.start[rownum+1] <- seg.out$loc.start[rownum]
          nextnum <- seg.out$num.mark[rownum+1]
          seg.out$num.mark[rownum+1] <- seg.out$num.mark[rownum+1] + seg.out$num.mark[rownum] ### total marker number
		  total.seg.values <- (seg.out$seg.mean[rownum+1] * nextnum + seg.out$seg.mean[rownum] * seg.out$num.mark[rownum])
          seg.out$seg.mean[rownum+1] <- total.seg.values / seg.out$num.mark[rownum+1]
        } else {
          seg.out.merge <- rbind(seg.out.merge, seg.out[rownum,])
        }
      } else {
        break
      }
    }

    if (!is.null(seg.out.merge)) {
      seg.out.merge <-  seg.out.merge[, -ncol(seg.out.merge)]
    }
  }
  ### output
  seg.out.merge
}

################################################################################################
# Genotype segmentation
################################################################################################
#' Function for segmentation
#'
#' DH segmentation using genotyping data
#' Input data requirment: genotyping data with columns of chr; pos; >=1 lines
#' The package DNAcopy is required.
#'
#' @param geno genotyping data, containing at least chr; pos; genotyping scores of lines
#' @param genocols column numbers of individuals
#' @param chromosomes chromosomes or contigs of interest for segmentation
#' @param chrname column name/num of chromosomes or contigs
#' @param posname column name/num of marker positions
#' @param data.type binary or logratio, default logratio
#' @param missing.name code for missing data, could be a vector
#' @param alleleX.name allele labels in the output
#' @param alleleX.code allele codes as values for segmentation
#' @param output.common the prefix name of output file
#' @param min.seg.size minimal size of segments, segments with smaller length will be ignored
#' @param seg.mean.cutoffs a numeric vector with two values, c(smaller, larger); default=c(-0.8, 0.8);
#'        default values correspond -1 and 1 as input values of allele1.code and allele2.code,respectively.
#'                  segment mean values <= smaller: allele 1
#'                  segment mean values >= larger: allele 2
#' @param cna.xxx are DNAcopy parameters, please refer DNAcopy mannual
#' @return output file and variable of the genotyping segmentation result
#' @author Sanzhen Liu
#' @seealso \url{https://bioconductor.org/packages/release/bioc/html/DNAcopy.html}
#' @export
#' @examples genosegDH(geno = "")
#'
genosegDH <- function (geno, genocols, chromosomes, chrname = "chr", posname = "pos",
					output.common = "seg", data.type = "logratio",
                    allele1.name = "A", allele2.name = "B", missing.name = c("0"),
					allele1.code = -1, allele2.code = 1,
                    min.seg.size = 100000, cna.alpha = 0.01, cna.nperm = 10000,
                    cna.p.method = "perm", cna.eta = 0.01, cna.min.width = 5,
                    seg.mean.cutoffs = c(-0.8, 0.8)) {

  #if (!"DNAcopy" %in% installed.packages()) {
  #source("http://bioconductor.org/biocLite.R")
  #  biocLite("DNAcopy")
  #}
  #library(DNAcopy)

  ## output
  seg.file <- paste0(output.common, ".seg.txt")
  codes <- c(allele1.name, allele2.name)
  names(codes) <- c(allele1.code, allele2.code)

  allcolnames <- colnames(geno)
  all_lines.seg <- NULL
  all_lines <- allcolnames[genocols]
  # generate segments and breakpoints for each line:
  for (eachline in all_lines) {
    dsort <- geno[geno[, chrname] %in% chromosomes, c(chrname, posname, eachline)]
    dsort[, eachline] <- as.character(dsort[, eachline])
    dsort <- dsort[!dsort[, eachline] %in% missing.name, ] # exclude missing data

    ### convert code to number
    dsort[dsort[, eachline] == allele1.name, eachline] <- allele1.code
    dsort[dsort[, eachline] == allele2.name, eachline] <- allele2.code
    dsort[, eachline] <- as.numeric(dsort[, eachline])

    ### sort SNP via chromosome + position
    dsort <- dsort[order(dsort[, chrname], dsort[, posname]),]
    dsort <- dsort[dsort[, chrname] %in% chromosomes, ]
    ################################################################################################
    # segmentation: SNPs
    #-----------------------------------------------------------------------------------------------
    # segmentation:
    line.CNA <-CNA(dsort[, eachline], dsort[, chrname], dsort[, posname],
                  data.type = data.type, sampleid = eachline)
    # No smoothness for the binary data

    # core step for the segmentation:
    line.CNA.segment <- segment(line.CNA, alpha = cna.alpha, nperm = cna.nperm,
                               p.method = cna.p.method, eta = cna.eta, min.width = cna.min.width);

    # output:
    rec.seg <- segments.summary(line.CNA.segment)

    # merge the intervals that are <= min.seg.size
    # to avoid the assembly error that create pseudo segments, merge some interval:
    rec.seg.merge <- NULL
    for (chrom in chromosomes) {
      chrseg <- rec.seg[rec.seg$chrom == chrom,]

      ### only execute if dataframe has >=1 entries
      if (nrow(chrseg) >= 1) {
          chr.largeseg <- chrseg[(chrseg$loc.end - chrseg$loc.start + 1) >= min.seg.size
                          & (chrseg$seg.mean <= seg.mean.cutoffs[1]
                          | chrseg$seg.mean >= seg.mean.cutoffs[2]), ]
        #### print(chr.largeseg) #####
        if (nrow(chr.largeseg) >= 1) {
          #print(chr.largeseg)
          chr.largeseg.merge <- seg.merge(chr.largeseg)
          ### combine all chromosomes:
          if (!is.null(chr.largeseg.merge)) {
            rec.seg.merge <- rbind(rec.seg.merge, chr.largeseg.merge)
          }
        }
      }
    }

    #nrow(rec.seg.merge)
    #rec.seg.merge <- rec.seg.merge[, -10]
    median.mismatch <- !(rec.seg.merge$seg.median %in% c(allele1.code, allele2.code))
    if (sum(median.mismatch) > 0) {
      warning(paste(median.mismatch, " segments whose genotypes were not clearly determined were removed\n"))
      print(rec.seg.merge[median.mismatch, ])
      rec.seg.merge <- rec.seg.merge[!median.mismatch, ]
    }
    #print(head(rec.seg.merge$seg.median))
    rec.seg.merge$Genotype <- codes[as.character(rec.seg.merge$seg.median)]
    all_lines.seg <- rbind(all_lines.seg, rec.seg.merge)
  }

  # output the result:
  colnames(all_lines.seg) <- c("Individual", "Chr", "segStart", "segEnd", "numMarker",
                               "segMean", "segSD", "segMedian", "segMAD", "Genotype")
  write.table(all_lines.seg, seg.file, row.names=F, quote=F, sep="\t")
  invisible(all_lines.seg)
}

