#' WID clocks
#'
#' THE WID_clocks function uses a methylation beta matrix and returns the WID general, epithelial and immune clocks.
#' It also returns the WID-relative-epithelial-age (WID_rea) and WID-relative-immune-age (WID_ria) (∆WID epithelial or WID immune and WID general clock).
#'
#' Authors: Barrett J., Herzog C. & Kim. Y.-N. et al., 2021
#'
#' @param beta Methylation beta matrix object
#' @param tissue Tissue type; if "breast", the epithelial-fibroblast-fat-immune reference will be used in the EpiDISH algorithm. Defaults to NULL
#' @return Returns a list object with WID clocks.
#' @export

WID_clocks <- function(beta, tissue = NULL){

  #------------------------------#
  # 1. WID general clock
  #------------------------------#

  # load coefficients
  w <- WID_general_clock_coef
  intercept <- w[1]
  w <- w[2:length(w)]

  # compute index
  intersect <- intersect(names(w), rownames(beta))

  if(!identical(intersect, names(w))){
    cat("Warning: not all WID general clock CpGs found in your beta matrix. This could be because you are using the 450k array instead of the EPIC array.\nIndex will still be computed.")
  }

  B <- beta[match(intersect, rownames(beta)),] # subset beta by index CpGs
  w <- w[match(intersect, names(w))]

  B1 <- B*w

  WID_general_clock <- intercept + apply(B1, MARGIN = 2, FUN = 'sum')

  #------------------------------#
  # 2. WID epithelial clock
  #------------------------------#

  # compute IC using epidish
  #     using centEpiFibIC.m as default - tissue is by default NOT breast

  if(is.null(tissue)){
    out.l <- epidish(beta.m = beta,
                     ref.m = centEpiFibIC.m,
                     method = "RPC")
    ic <- out.l$estF[,3]
  } else if(tissue=='breast'){
    out.l <- epidish(beta.m = beta,
                     ref.m = centEpiFibFatIC.m,
                     method = "RPC")
    ic <- out.l$estF[,4]
  }

  # load coefficients + prep matrix
  w <- WID_epithelial_clock_coef
  intercept <- w[1]
  w <- w[2:length(w)]
  intersect <- intersect(names(w[1:168]), rownames(beta))

  B <- beta[match(intersect, rownames(beta)),] # subset beta by index CpGs

  B_ic <- t(t(B)*ic)   # add interaction terms
  rownames(B_ic) <- paste(rownames(B_ic),
                          "_ic", sep="")
  beta_with_interaction <- rbind(B, B_ic)

  # compute index
  ind1 <- match(names(w),rownames(beta_with_interaction))
  w2 <- w[!is.na(ind1)]

  WID_epithelial_clock <- intercept + apply(beta_with_interaction*w2,
                                            MARGIN = 2,
                                            FUN = 'sum')

  #------------------------------#
  # 3. WID immune clock
  #------------------------------#

  out.l <- epidish(beta.m = beta,
                   ref.m = centEpiFibIC.m,
                   method = "RPC")
  ic <- out.l$estF[,3]

  # load coefficients + prep matrix
  w <- WID_immune_clock_coef
  intercept <- w[1]
  w <- w[2:length(w)]
  intersect <- intersect(names(w[1:56]), rownames(beta))

  B <- beta[match(intersect, rownames(beta)),] # subset beta by index CpGs

  B_ic <- t(t(B)*ic)   # add interaction terms
  rownames(B_ic) <- paste(rownames(B_ic),
                          "_ic", sep="")
  beta_with_interaction <- rbind(B, B_ic)

  # compute index
  ind1 <- match(names(w),rownames(beta_with_interaction))
  w2 <- w[!is.na(ind1)]

  WID_immune_clock <- intercept + apply(beta_with_interaction*w2,
                                        MARGIN = 2,
                                        FUN = 'sum')

  #------------------------------#

  return(list(WID_general_clock = WID_general_clock,
              WID_epithelial_clock = WID_epithelial_clock,
              WID_immune_clock = WID_immune_clock,
              WID_rea = WID_epithelial_clock - WID_general_clock,
              WID_ria = WID_immune_clock - WID_general_clock)
  )
}
