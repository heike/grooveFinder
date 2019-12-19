context("get_grooves_hough")




tmp <- matrix(0, nrow = 40, ncol = 300)

tmp[, c(41, 273)] <- 1
tmp[, c(43, 271)] <- -1

tmp_df <- data.frame(y = rep(1:40, times = 300) * .625, x = rep(1:300, each = 40) * .625, value = tmp[1:length(tmp)])
exnames <- c("left.groove.fit", "right.groove.fit")


# Groove 1 is at 26.875 microns
# Groove 2 is at 169.375 microns

tmp_x3p <- x3ptools::df_to_x3p(tmp_df)


tmp_mid <- matrix(0, nrow = 40, ncol = 300)

tmp_mid[, c(81, 220)] <- 1
tmp_mid[, c(83, 218)] <- -1

tmp_mid_df <- data.frame(y = rep(1:40, times = 300) * .625, x = rep(1:300, each = 40) * .625, value = tmp_mid[1:length(tmp_mid)])
exnames <- c("left.groove.fit", "right.groove.fit")
tmp_mid_x3p <- x3ptools::df_to_x3p(tmp_mid_df)

# groove 1 will be placed at 51.875
# groove 2 will be placed at 136.25

tmp2 <- matrix(0, nrow = 40, ncol = 300)


left.groove.fit <- function(yinput) {
  bottom.index <- 15

  left.groove <- (bottom.index + .25 * yinput)

  return(left.groove)
}

right.groove.fit <- function(yinput) {
  bottom.index <- 255 # scale bottom.right to microns

  right.groove <- (bottom.index + .25 * yinput)

  return(right.groove)
}

left <- left.groove.fit(1:(nrow(tmp2)))
left <- floor(left)
right <- right.groove.fit(1:(nrow(tmp2)))
right <- floor(right)

tmp2 <- sapply(length(left):1, FUN = function(i) {
  tmp2[ i, left[i]] <- 1
  tmp2[ i, left[i] + 2] <- -1
  tmp2[i, ]
}) %>% t()

tmp2 <- sapply(length(right):1, FUN = function(i) {
  tmp2[i, right[i]] <- 1
  tmp2[i, right[i] + 2] <- -1
  tmp2[i, ]
}) %>% t()


tmp2_df <- data.frame(y = rep(1:40, times = 300) * .625,
                      x = rep(1:300, each = 40) * .625,
                      value = tmp2[1:length(tmp2)])
tmp2_x3p <- df_to_x3p(tmp2_df)



test_that("get_grooves_hough works", {
  # for straight lines
  grooves <- grooveFinder::get_grooves_hough(x3ptools::x3p_to_df(tmp_x3p))
  expect_is(grooves, "list")
  expect_is(grooves[[1]], "function")
  expect_is(grooves[[2]], "function")
  expect_is(grooves$left.groove.fit(100), "numeric")
  expect_lte(grooves$left.groove.fit(6.25), 26.875 + 11)
  expect_gte(grooves$left.groove.fit(6.25), 26.875)
  expect_gte(grooves$right.groove.fit(6.25), 169.375 - 11)
  expect_lte(grooves$right.groove.fit(6.25), 169.375)
  # slanted lines
  grooves2 <- grooveFinder::get_grooves_hough(x3ptools::x3p_to_df(tmp2_x3p))
  expect_is(grooves2, "list")
  expect_is(grooves2[[1]], "function")
  expect_is(grooves2[[2]], "function")
  expect_is(grooves2$left.groove.fit(100), "numeric")
  # lines within middle 50
  grooves3 <- grooveFinder::get_grooves_hough(x3ptools::x3p_to_df(tmp_mid_x3p))
  expect_is(grooves3, "list")
  expect_is(grooves3[[1]], "function")
  expect_is(grooves3[[2]], "function")
  # groove 1 will be placed at 51.875
  # groove 2 will be placed at 136.25
  expect_is(grooves3$left.groove.fit(100), "numeric")
  expect_gte(grooves3$left.groove.fit(6.25), 51.875)
  expect_lte(grooves3$left.groove.fit(6.25), 51.875 + 11)
  expect_lte(grooves3$right.groove.fit(6.25), 136.25)
  expect_gte(grooves3$right.groove.fit(6.25), 136.25 - 11)


  # test fit around bottom index
  expect_lte(grooves2$left.groove.fit(0), 10.625 + 11)
  expect_gte(grooves2$left.groove.fit(0), 10.625)
  expect_gte(grooves$right.groove.fit(0), 160.625 - 11)
  expect_lte(grooves$right.groove.fit(0), 160.625)
  # text fit around top index
  expect_lte(grooves2$left.groove.fit(40), 16.875 + 11)
  expect_gte(grooves2$left.groove.fit(40), 16.875)
  expect_gte(grooves$right.groove.fit(40), 166.875 - 11)
  expect_lte(grooves$right.groove.fit(40), 166.875)
  expect_equivalent(names(grooves2), exnames)
  # Test that both x3p and dfs work in function
  expect_silent(grooves <- get_grooves_hough(tmp_df))
  expect_silent(grooves <- get_grooves_hough(tmp_x3p))
  expect_error(get_grooves_hough(c(1,2)), "land should be a data frame or an x3p object")
})
