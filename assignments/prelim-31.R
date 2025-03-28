vec <- function(...) {
  vals <- list(...) |> unlist()
  matrix(vals, ncol=1)
}
rvec <- function(first, ..., rfun=rnorm) {
  dots <- c(first, list(...)) |> unlist()
  if (length(dots) == 1) {
    dots <- rfun(dots[[1]])
  }

  vec(dots)
}
