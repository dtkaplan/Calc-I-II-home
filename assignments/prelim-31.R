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

show_matrix <- function(M) {
  image(1:ncol(M), 1:nrow(M), t(M),
        col = colorspace::diverging_hcl(
          15, h = c(180, 50), c = 80,
          l = c(20, 95), power = c(0.7, 1.3)),
        axes = FALSE)
}

