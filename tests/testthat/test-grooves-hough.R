context("get_grooves_hough")




tmp <- matrix(0, nrow = 40, ncol = 300)

tmp[, c(41, 273)] <- 1
tmp[,c(43, 271)] <- -1

tmp_df <- data.frame(y = rep(1:40, times = 300)*.625, x = rep(1:300, each = 40)*.625, value = tmp[1:length(tmp)])
exnames <- c("left.groove.fit", 'right.groove.fit')


# Groove 1 is at 26.875 microns
# Groove 2 is at 169.375 microns

tmp_x3p <- x3ptools::df_to_x3p(tmp_df)

# Test that both x3p and dfs work in function

expect_silent(grooves <- get_grooves_hough(tmp_df))

expect_silent(grooves <- get_grooves_hough(tmp_x3p))

# test_vec <- c(1,2)
# expect_error(get_grooves_hough(test_vec), "Land input is neither an x3p or a dataframe with necessary bullet land data.")

tmp2 <- matrix(0, nrow = 40, ncol = 300)


left.groove.fit <- function(yinput){
  bottom.index <- 15

  left.groove <- (bottom.index + .25*yinput)

  return(left.groove)
}
right.groove.fit <- function(yinput) {

  bottom.index <- 255 # scale bottom.right to microns

  right.groove <- (bottom.index  + .25*yinput)

  return(right.groove)
}

left <- left.groove.fit(1:(nrow(tmp2)))
left <- floor(left)
right <- right.groove.fit(1:(nrow(tmp2) ))
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


tmp2_df <- data.frame(y = rep(1:40, times = 300)*.625, x = rep(1:300, each = 40)*.625, value = tmp2[1:length(tmp2)])
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
})
