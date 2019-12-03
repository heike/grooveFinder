#' Use logistic model to identify groove locations
#'
#' @param x numeric vector of locations (in microns)
#' @param value numeric values of surface measurements in microns
#' @param lasso_method use the 'basic' model or the 'full' model with interaction terms?
#' @param pred_cutoff equal error rate cutoff for classification into GEA or LEA, trained on Hamby set 44
#' @param return_plot return plot of grooves?
#' @importFrom locfit locfit.robust
#' @importFrom locfit locfit
#' @importFrom stats model.matrix
#' @importFrom assertthat assert_that
#' @importFrom assertthat has_name
#' @export
get_grooves_lasso <- function(x, value, lasso_method = "basic", pred_cutoff = ifelse(lasso_method == "basic", .3, 0.34), return_plot = F) {
  land <- data.frame(x = x, value = value)
  original_land <- land

  assert_that(
    has_name(land, "x"), has_name(land, "value"),
    is.numeric(land$x), is.numeric(land$value)
  )

  assert_that(lasso_method %in% c("basic", "full"))
  assert_that(pred_cutoff > 0, pred_cutoff < 1)

  ## generate additional variables
  check_min <- min(land$value[!is.na(land$value)])
  land <- land %>% mutate(value_std = value - check_min)

  ## fit robust loess, calculate residuals
  rlo_fit <- robust_loess_fit(cc = land, iter = 20)
  land$rlo_pred <- predict(rlo_fit, newdata = land)
  land$rlo_resid <- land$value_std - land$rlo_pred
  land$rlo_absresid <- abs(land$rlo_resid)

  mx <- max(land$x, na.rm = T)
  diff_mx <- mx / 2 - land$x
  ## Use this method because some data may have shifted x values
  median <- land$x[which.min(abs(diff_mx))] # changed abs(tst_mx) to abs(diff_mx) - tst_mx not defined
  # median <- median(land$x) # some of the houston data appears to have a shifted x
  land$side <- "right"
  land$side <- ifelse(land$x <= median, "left", land$side)
  land$depth <- abs(land$x - median)

  ## range50 : range of values in a 50-wide band around each data point.
  land$range_50 <- rollapply(land$rlo_resid, width = 50, FUN = function(x) {
    max(x) - min(x)
  }, partial = TRUE)

  ## xint1 and xint2: the predicted locations that the robust LOESS crosses the x-axis.
  xint1 <- min(abs(land$rlo_pred[(land$x < median(land$x))]), na.rm = T)
  xint2 <- min(abs(land$rlo_pred[(land$x > median(land$x))]), na.rm = T)
  ind1 <- which(land$rlo_pred == xint1 | land$rlo_pred == -1 * xint1)
  ind2 <- which(land$rlo_pred == xint2 | land$rlo_pred == -1 * xint2)
  land$xint1 <- land$x[ind1]
  land$xint2 <- land$x[ind2]

  land$ind_edges <- ifelse(land$x < land$xint1 | land$x > land$xint2, 1, 0)

  ## ind_2mad: whether the data point is above the 2*MAR cutoff previously used as an ad-hoc method.
  mad <- mad(land$rlo_resid, na.rm = T)
  land$ind_2mad <- ifelse(land$rlo_resid > 2 * mad, 1, 0)

  ## numpos_50: how many positive residuals there are in a 50-wide band around each data point.
  land$numpos_50 <- rollapply(land$rlo_resid, width = 50, FUN = function(x) {
    sum(x > 0)
  }, partial = TRUE)

  land$numNA_50 <- rollapply(land$rlo_resid, width = 50, FUN = function(x) {
    sum(is.na(x))
  }, partial = TRUE)
  lower <- quantile(land$x, prob = .25)
  upper <- quantile(land$x, prob = .75)
  proxy_dat <- land %>% filter(x < upper & x > lower)
  proxy <- sd(proxy_dat$rlo_resid, na.rm = T)
  land$rlo_resid_std <- land$rlo_resid / proxy
  land$range_50_std <- land$range_50 / proxy

  xrange <- max(land$x) #- min(land$x) # again correcting for shifted houston persistence data
  land$depth_std <- land$depth / xrange
  land$xint1_std <- land$xint1 / xrange
  land$xint2_std <- land$xint2 / xrange

  ## now get logistic predictions
  land <- na.omit(land)

  if (lasso_method == "basic") {
    X <- cbind(1, model.matrix(
      ~ rlo_resid_std + I(rlo_resid_std^2) + side +
        depth_std + side * depth_std + xint1_std +
        xint2_std + range_50 + numNA_50 + ind_2mad +
        numpos_50 + ind_edges - 1,
      land
    ))
    ymean <- X %*% lasso_simple
  } else if (lasso_method == "full") {
    X <- cbind(1, model.matrix(
      ~ (rlo_resid_std + I(rlo_resid_std^2) + side +
        depth_std + xint1_std +
        xint2_std + range_50 + numNA_50 + ind_2mad +
        numpos_50 + ind_edges)^2 - 1,
      land
    ))
    ymean <- X %*% lasso_interactions
  } else {
    stop("invalid lasso_method argument.")
  }

  yhat <- as.vector(exp(ymean) / (1 + exp(ymean)))
  land$pred_val <- yhat
  land$pred_class <- ifelse(land$pred_val < pred_cutoff, "LEA", "GEA")

  groove <- range(land$x[land$pred_class == "LEA"])

  if (return_plot) {
    return(list(
      groove = groove,
      plot = grooves_plot(land = original_land, grooves = groove)
    ))
  } else {
    return(list(groove = groove))
  }
}
