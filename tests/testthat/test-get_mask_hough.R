context("get_mask_hough")

test.image <- data.frame(x = c(rep(1, 8), rep(2,8), rep(3,8), rep(4,8), rep(5, 8), rep(6, 8),
                               rep(7,8), rep(8,8), rep(9,8), rep(10,8), rep(11, 8), rep(12, 8),
                               rep(13, 8), rep(14,8), rep(15,8), rep(16,8), rep(17,8), rep(18,8),
                               rep(19,8), rep(20,8), rep(21,8), rep(22,8), rep(23,8), rep(24,8)),
                         y = c(rep(seq(1,8), 24)),
                         value = c(rep(1,8), rep(1,7), 0, rep(1,5), rep(0,3), 1, rep(0, 7),
                                   rep(1,8), rep(1,8), rep(1,8), rep(1,8), rep(1,8), rep(1,8),
                                   rep(1,8), rep(1,8), rep(1,8), rep(1,8), rep(1,8), rep(1,8),
                                   rep(1,8), rep(1,8), rep(1,8), rep(1,8), rep(1, 5), rep(0,3),
                                   rep(0,3), rep(0,5), 0, rep(1,7), rep(1,8)))

test.x3p <- df_to_x3p(test.image)

# Tried to create small test image but get_mask_hough kept throwing an error:
# mask works
# only 0's may be mixed with negative subscripts
# but if you run it line by line a """"Valid"""" (looks terrible) mask is created

data(br411)
x3p <- br411

exnames <- c("left.groove.fit", "right.groove.fit")
b4mask <- c("header.info", "surface.matrix", "feature.info", "general.info", "matrix.info")
aftermask <- c("header.info", "surface.matrix", "feature.info", "general.info", "matrix.info", "mask")

test_that("mask works", {
  grooves <- grooveFinder::get_grooves_hough(x3p_to_df(x3p))
  expect_equivalent(names(x3p), b4mask)
  new.mask <- grooveFinder::get_mask_hough(x3p, grooves)
  expect_is(grooves, "list")
  expect_is(grooves[[1]], "function")
  expect_is(grooves[[2]], "function")
  expect_is(grooves$left.groove.fit(100), "numeric")
  expect_equivalent(names(new.mask),aftermask)



})
