################################################################################################
# convert physical positions to genetic positions based on a map
################################################################################################
#' Function for a physcial to genetic conversion
#'
#' @method Here is a brief algorithm description:
#' 1. The genetic and physical curve is locally smoothened using a Generalized Additive Model (gam)
#' 2. Based on the smoothened curve, predict genetic positions for physcial positions
#'
#' @param query  a data frame with at least 2 columns: chromosome and position
#' @param map    a data frame containing a physical and genetic map (4 columns):
#'               1.marker_name; 2.chromosome; 3. physical position; 4. genetic position
#' @param qchr   column (name or number) of chromosomes (1)
#' @param qpos   column (name or number) of chromosome physical positions (2)
#' @param pred   prediction options: 1. "gen" (genetic position) and 2: "rec" (recombination rate) (default=gen)
#' @param cname  colname of the column adding genetic information ("genetic")
#' @param adjust logic value to specify whether predicted genetic values are needed to adjusted (TRUE)
#' @param prefix prefix for output files
#' @param mb     region context in megabase to determine recombination rates if pred="rec" (1)
#' @param k      a smooth factor used in gam function (50)
#' @param mplot  logic value to indicate whether genetic_vs_physical plots per chromosome are produced (TRUE)
#' @return output file and variable of the genotyping segmentation result
#' @author Sanzhen Liu (liu3zhen@gmail.com)
#' @seealso \url{https://cran.r-project.org/web/packages/mpcv/index.html}
#' @details The R script uses a map containing genetic and physical data to identify local models that
#' smoothen the curve and can predict genetic positions for any given physical positions.
#'
#' The script include model fitting via generalized additive models (gam), prediction, adjustment of prediction if
#' "adjust=T" specified. The performance will be affected by the quality and the number of markers of the map.
#'
#' The "gam" function in the package mgcv was used for model fitting. The parameter of "k" is the parameter used in gam.
#' It is "the dimension of the basis used to represent the smooth term.". More details can refer to the function of "s" in mgcv.
#' For most cases, k=50 works well. The result of model fitting can be visualized using observation vs. prediction plots that
#' can be turned on with "mplot=T". Each chromosome has one PDF plot output.
#'
#' Two modes of prediction are designed, including gen and rec, by using the parameter of "pred".
#' The "gen" mode predicts and outputs genetic positions of each physical postions.
#' The "rec" mode predicts and outputs recombination distance per interval around each physical position of interest.
#' The interval length is specified by the parameter of "mb" and the physical position of interest is at the midpoint of
#' of the interval.
#'
#' The reason for adjustment is that the map often starts with a physical position greater than 1 but the genetic positions are
#' usually start with 0. Therefore, prediction at some small physical positions are likely to be negative. Using "adjust=T",
#' all genetic positions are adjusted by substracting the genetic value of physical position of 1 if
#' that genetic value is negative.
#' @export
#' @usage phy2gen(query, map)
#' @examples
#===========================================================================================#
phy2gen <- function (query, map, qchr=1, qpos=2, pred="gen", cname="genetic",
                     adjust = T, prefix="gen_pred", mb=1, k=50, mplot = T) {
  # check parameters:
  stopifnot((is.out & !is.null(outfile)) | !is.out)
  stopifnot(nrow(query) > 0)
  stopifnot(ncol(query) >= 2) # at least 2 columns
  stopifnot(nrow(map) > 0)
  stopifnot(ncol(map) >= 4) # at least 4 columns

  # load required package(s)
  library(mgcv)

  # change map header:
  original_map_header <- colnames(map)
  colnames(map) <- c("marker", "chr", "physical", "genetic")

  # initiate ouput:
  new_query <- NULL
  new_map <- NULL
  if (adjust) {
    adjust_outfile <- paste0(prefix, ".adjust.cM")
    cat("chr", "adjust", sep = "\t", file=adjust_outfile)
    cat("\n", file=adjust_outfile, append=T)
  }

  qry_chrs <- sort(unique(query[, qchr]))
  map_chrs <- unique(map$chr)
  # go through each chromosome in query
  for (chr in qry_chrs) {
    # extract marker data for certain chromosome:
    cur_map <- map[map$chr == chr, ]
    stopifnot(nrow(cur_map) > 0)

    # GAM fitting
    gam_fit <- gam(genetic ~ s(physical, k=k), data = cur_map)

    #####################################################################################
    ### plot both observed map and predicted map
    #####################################################################################
    if (mplot) {
      pdf(paste0(prefix, "_", chr, ".map.pdf"), width = 6, height = 6)
      plot(cur_map$physical, cur_map$genetic, main = paste("obs vs. pred -", chr),
           xlab = "Physical (bp)", ylab = "Genetic (cM)", cex = 0.35, col = "gray30")
      points(cur_map$physical, gam_fit$fitted.values, cex = 0.25, col = "orange")
      legend("bottomright", legend = c("observed", "fitted"), col = c("gray30", "orange"),
             pch = 1, bty = "n")
      dev.off()
    }

    #####################################################################################
    ### adjust
    ####################################################################################
    # adjust
    genetic_adjust <- 0
    if (adjust) {
      # predict genetic position for position 1
      chr_start_point <- data.frame(physical=1)
      genetic_adjust <- as.numeric(predict.gam(gam_fit, chr_start_point))
      genetic_adjust <- min(genetic_adjust, 0) # if >0, no adjust adjusted
      cat(chr, genetic_adjust, sep = "\t", file=adjust_outfile, append=T)
      cat("\n", file=adjust_outfile, append=T)
    }

    ####################################################################################
    # prediction:
    ####################################################################################
    cur_qry <- query[query[, qchr] == chr, ]
    if (pred == "gen") {
      phy <- data.frame(physical = cur_qry[, qpos]) # physical position
      gen <- round(predict.gam(gam_fit, phy) - genetic_adjust, 2) # with adjust
      cur_qry[, cname] <- gen
    } else if (pred == "rec") {
      phy1 <- data.frame(physical = cur_qry[, qpos] - mb * 1000000 / 2) # pos 1 of interval
      phy2 <- data.frame(physical = cur_qry[, qpos] + mb * 1000000 / 2) # pos 1 of interval
      gen1 <- predict.gam(gam_fit, phy1)
      gen2 <- predict.gam(gam_fit, phy2)
      gen_dist <- gen2 - gen1
      gen_dist <- max(gen_dist, 0)
      cur_qry[, cname] <- round(gen_dist, 3)
    }

    ####################################################################################
    # merge prediction
    ####################################################################################
    if (is.null(new_query)) {
      new_query <- cur_qry
    } else {
      new_query <- rbind(new_query, cur_qry)
    }

    # updated map genetic position
    cur_map$genetic <- cur_map$genetic - genetic_adjust

    # merge new map
    if (pred == "gen") {
      if (is.null(new_map)) {
        new_map <- cur_map
      } else {
        new_map <- rbind(new_map, cur_map)
      }
    }
  }

  ####################################################################################
  # output new_map if pred=="gen"
  ####################################################################################
  # if chrs are not in query, add those map to new_map
  if (pred == "gen") {
    chrs_not_in_query <- map_chrs[! map_chrs %in% qry_chrs]
    if (length(chrs_not_in_query) > 0) {
      rest_map <- map[map$chr %in% chrs_not_in_query, ]
      new_map <- rbind(new_map, rest_map)
    }
    new_map_outfile <- paste0(prefix, ".new.map")
    colnames(new_map) <- original_map_header
    write.table(new_map, new_map_outfile, quote=F, row.names=F, sep="\t")
  }

  ####################################################################################
  # output predicted results
  ####################################################################################
  new_qry_outfile <- paste0(prefix, ".pred.gen")
  write.table(new_query, new_qry_outfile, quote=F, row.names=F, sep="\t")

  invisible(new_query)
}
