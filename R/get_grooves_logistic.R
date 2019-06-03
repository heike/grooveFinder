#' Use logistic model to identify groove locations
#'
#' @inheritParams get_grooves_quadratic
#' @importFrom locfit locfit.robust
#' @importFrom locfit locfit
#' @importFrom stats model.matrix
#' @export
get_grooves_logisticlegacy <- function(x, value, adjust = 10, # smoothfactor = 15,
                                       # groove_cutoff = 400,
                                       return_plot = F) {
  land <- data.frame(x = x, value = value)
  original_land <- land

  ## generate additional variables

  check_min <- min(land$value[!is.na(land$value)])
  land <- mutate(land, value_std = value - check_min)
  # install.packages("locfit")
  # library(locfit)
  robust_loess_fit <- locfit.robust(value_std ~ x, data = land, alpha = 1, kern = "tcub")
  land$rlo_pred <- predict(robust_loess_fit, newdata = land)

  land$rlo_absresid <- with(land, abs(value_std - rlo_pred))
  land$rlo_resid <- with(land, value_std - rlo_pred)


  median <- median(land$x)
  land$side <- "right"
  land$side <- ifelse(land$x <= median, "left", land$side)
  land$depth <- abs(land$x - median)

  ## range20 : range of values in a 20-wide band around each data point.
  land$range_50 <- rollapply(land$rlo_resid, width = 50, FUN = function(x) {
    max(x) - min(x)
  }, partial = TRUE)

  ## xint1 and xint2: the predicted locations that the robust LOESS crosses the x-axis.
  xint1 <- min(abs(land$rlo_pred[(land$x < median(land$x))]))
  xint2 <- min(abs(land$rlo_pred[(land$x > median(land$x))]))
  ind1 <- which(land$rlo_pred == xint1 | land$rlo_pred == -1 * xint1)
  ind2 <- which(land$rlo_pred == xint2 | land$rlo_pred == -1 * xint2)
  land$xint1 <- land$x[ind1]
  land$xint2 <- land$x[ind2]

  ## ind_2mad: whether the data point is above the 2*MAR cutoff previously used as an ad-hoc method.
  mar <- median(land$rlo_absresid, na.rm = T)
  land$ind_2mad <- ifelse(land$rlo_absresid > 2 * mar, 1, 0)

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

  xrange <- max(land$x) - min(land$x)
  land$depth_std <- land$depth / xrange
  land$xint1_std <- land$xint1 / xrange
  land$xint2_std <- land$xint2 / xrange

  ## now get logistic predictions
  model_all4 <- as.matrix(c(-26.7166509, 0.1727030, 0, -0.1815079, 0, 39.7340095, -1.0473396, 7.0916175, 0.2428548, 0, 1.6039295, 0, 0))


  land <- na.omit(land)
  X <- cbind(1, model.matrix(
    ~rlo_resid_std + I(rlo_resid_std^2) + side +
      depth_std + side * depth_std + xint1_std +
      xint2_std + range_50 + numNA_50 + ind_2mad +
      numpos_50 - 1,
    land
  ))
  ymean <- X %*% model_all4
  yhat <- exp(ymean) / (1 + exp(ymean))
  land$pred_val <- yhat
  land$pred_class <- ifelse(land$pred_val < .25, "LEA", "GEA")

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


#' Fit a robust loess regression
#'
#' Internal function called by get_grooves_lassobasic and get_grooves_lassofull
#' @param cc data frame with columns x and value_std, representing the crosscut
#' @param iter number of iterations
#' @importFrom stats loess
#' @importFrom assertthat assert_that
#' @importFrom assertthat has_name
robust_loess_fit <- function(cc, iter) {
  assert_that(has_name(cc, "x"), has_name(cc, "value_std"))
  n <- nrow(cc)
  weights <- rep(1, n)
  fit <- loess(value_std ~ x, data = cc, span = 1)
  cc$fit <- predict(fit, newdata = cc)
  cc$resid <- cc$value_std - cc$fit
  i <- 1
  while (i < iter) {
    mar <- median(abs(cc$resid), na.rm = T)
    cc$bisq <- pmax(1 - (cc$resid / (6 * mar))^2, 0)^2
    weights <- ifelse(cc$resid > 0, cc$bisq, 1)
    fit <- loess(value_std ~ x, data = cc, span = 1, weights = weights)
    cc$fit <- predict(fit, newdata = cc)
    cc$resid <- cc$value_std - cc$fit
    i <- i + 1
  }
  return(fit)
}
