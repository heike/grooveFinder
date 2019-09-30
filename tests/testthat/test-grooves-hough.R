context("get_grooves_hough")

data(br411)
x3p <- br411
exnames <- c("left.groove.fit", "right.groove.fit")

test_that("get_grooves_hough works", {
  grooves <- grooveFinder::get_grooves_hough(x3ptools::x3p_to_df(x3p))
  expect_is(grooves, "list")
  expect_is(grooves[[1]], "function")
  expect_is(grooves[[2]], "function")
  expect_is(grooves$left.groove.fit(100), "numeric")
  expect_equivalent(names(grooves), exnames)
})
