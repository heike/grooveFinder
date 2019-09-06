context("test-grooves-hough")

data(br411)
x3p <- br411

test_that("get_grooves_hough works", {
  grooves <- grooveFinder::get_grooves_hough(x3ptools::x3p_to_df(x3p), qu = 0.995)
  expect_is(grooves, "list")
  expect_is(grooves[[1]], "function")
  expect_is(grooves[[2]], "function")
  expect_is(grooves$left.groove.fit(100), "numeric")
})
