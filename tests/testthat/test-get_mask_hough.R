context("get_mask_hough")

tmp <- matrix(0, nrow = 40, ncol = 300)

tmp[, c(41, 273)] <- 1
tmp[,c(43, 271)] <- -1

tmp_df <- data.frame(y = rep(1:40, times = 300)*.625, x = rep(1:300, each = 40)*.625, value = tmp[1:length(tmp)])

# Groove 1 is at 26.875 microns
# Groove 2 is at 169.375 microns

tmp_x3p <- x3ptools::df_to_x3p(tmp_df)

test_that("mask works", {
  grooves <- grooveFinder::get_grooves_hough(x3p_to_df(tmp_x3p))
  new.mask <- grooveFinder::get_mask_hough(tmp_x3p, grooves)

  expect_is(grooves, "list")
  expect_is(grooves[[1]], "function")
  expect_is(grooves[[2]], "function")
  expect_lte(grooves$left.groove.fit(6.25), 26.875 + 11)
  expect_gte(grooves$left.groove.fit(6.25), 26.875)
  expect_gte(grooves$right.groove.fit(6.25), 169.375 - 11)
  expect_lte(grooves$right.groove.fit(6.25), 169.375)
  expect_true("mask" %in% names(new.mask))

})
