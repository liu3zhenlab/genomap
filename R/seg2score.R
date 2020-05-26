################################################################################################
# Genotype segmentation
################################################################################################
#' Function for generating common segment markers across all iindividuals
#'
#' read genotyping segmentation result and divide each segment
#' into multiple smaller segments. Basically, breakpoints from
#' all individuals were used to split original segments.
#'
#' @param seg.input input data that are the result from seggenoXX moduless
#' @param segscore.output output file for genotyping data
#' @param missing.data.code missing data code to replace NA. default=NULL.
#' @param binsize bin size for facilitating computation, default = 1000000. NOT recommended to change
#' @return genotyping table for individuals
#' @author Sanzhen Liu
#' @seealso seggenoF2
#' @export
#' @examples seg2score(seg.input = "", segscore.output = "")
#'
##
seg2score <- function(seg.input, segscore.output, missing.data.code = NULL, binsize = 1000000) {
  ### step1: read genotyping segmentation result and divide each segment
  ### into multiple smaller segments. Basically, breakpoints from
  ### all individuals were used to split original segments.
  seg <- read.delim(seg.input)
  chromosomes <- unique(seg$Chr)
  segmarker <- NULL
  for (eachchr in chromosomes) {
    current <- subset(seg, Chr == eachchr)
    allpoint <- c(current[,"segStart"], current[,"segEnd"])
    allpoint <- sort(unique(allpoint))
    allpoint[1] <- max(allpoint[1] - 1, 1)
    count <- 0
    for (i in 2:length(allpoint)) {
      count <- count + 1
      marker <- paste(eachchr,"_seg",count,sep="")
      marker.pos <- round((allpoint[i-1] + 1 + allpoint[i]) / 2)
      currentmarker <- c(marker, eachchr, allpoint[i-1] + 1, allpoint[i], marker.pos)
      names(currentmarker) <- c("Marker", "Chr", "Start", "End", "Pos")
      segmarker <- rbind(segmarker, currentmarker)
    }
  }
  rownames(segmarker) <- 1:nrow(segmarker)
  segmarker <- data.frame(segmarker)

  ### step 2: scoring for segmarkers
  overlap <- function(x) {
    ### determine whether two intervals overlap or not
    x <- as.numeric(as.character(x))
    overlap.val <- (min(x[1:2]) < max(x[3:4])) & (max(x[1:2]) > min(x[3:4]))
  }

  ### employ bin method to reduce the computation time
  segmarker$Bin <- round(as.numeric(as.character(segmarker$Start)) / binsize, 0)
  seg$Bin1 <- round(as.numeric(as.character(seg$segStart)) / binsize, 0)
  seg$Bin2 <- round(as.numeric(as.character(seg$segEnd)) / binsize, 0)
  ### check overlaps
  allgeno <- NULL
  for (eachchr in chromosomes) {
    binrange <- range(c(seg$Bin1[seg$Chr == eachchr], seg$Bin2[seg$Chr == eachchr]))
    binrange
    for (i in binrange[1]:binrange[2]) {
      csegmarker <- segmarker[segmarker$Chr == eachchr & segmarker$Bin == i, ]
      if (nrow(csegmarker) > 0) {
        cseg <- seg[seg$Chr == eachchr & seg$Bin1 <= i & seg$Bin2 >= i, ]
        cmerge <- merge(csegmarker, cseg, by = "Chr", all = T)
        coverlap <- apply(cmerge[, c("Start", "End", "segStart", "segEnd")], 1, overlap)
        cmerge <- cmerge[coverlap, c("Chr", "Marker", "Start", "End", "Pos", "Individual", "Genotype")]
        ### bind all genotyping data
        allgeno <- rbind(allgeno, cmerge)
      }
    }
  }

  ### reorganize data
  cn <- colnames(allgeno)
  common.cn <- cn[1:5]
  allgeno2 <- NULL
  allind <- sort(unique(allgeno$Individual))
  for (eachind in allind) {
    indgeno <- allgeno[allgeno$Individual == eachind, -which(cn == "Individual")]
    colnames(indgeno) <- c(common.cn, eachind)

    ### merge all data:
    if (is.null(allgeno2)) {
      allgeno2 <- indgeno
    } else {
      allgeno2 <- merge(allgeno2, indgeno, by = common.cn, all = T)
      ### change NA to specified missing code
    }
  }

  if (!is.null(missing.data.code)) {
    for (eachind in allind) {
      allgeno2[, eachind] <- as.character(allgeno2[, eachind])
      allgeno2[is.na(allgeno2[, eachind]), eachind] <- missing.data.code
    }
  }
  ### output
  write.table(allgeno2, segscore.output, quote = F, row.names = F, sep = "\t")
  invisible(allgeno2)
}
