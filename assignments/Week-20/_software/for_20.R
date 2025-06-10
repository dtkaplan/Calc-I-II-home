#' Draw simple vector addition problem
#'
#' Run this once to generate a set of candidate target points,
#' then set <keep> to be the indices (A=1, B=2, ...) of the ones you
#' want to keep (e.g. `keep = c(2, 3, 8)`)
#'
#' @param theta angle (in degrees) between the vectors
#' @param len2 relative length of second vector with respect to first
#' @param keep which targets to keep of the 10-or-so generated
#' @param seed integer random seed to use

vec_addition_problem <- function(
    theta = 130, len2=1.5, keep = NULL, vecnames = c("v", "w"),
    labs=LETTERS, linewidth = 2,
    colors =
      sample(c("blue", "black", "tomato", "darkgreen", "darkgray"), 2),
    seed = 2834
    ) {
  set.seed(seed)
  v1 <- rvec(2)
  v1 <- v1 / veclen(v1)
  theta <- theta * pi/180
  rot <- len2*matrix(c(cos(theta), sin(theta),
                  -sin(theta), cos(theta)), nrow=2)
  v2 <- rot %*% v1
  vecframe <- matrix(
    c(0, 0, v1, 0, 0, v2), nrow=2, byrow=TRUE
  )  |> tibble::as_tibble()
  names(vecframe) <- c("x0", "y0", "x1", "y1")
  vecframe$label <- vecnames
  vecframe <- vecframe |> dplyr::mutate(tx = x1*0.5, ty=y1*0.6)
  vecframe$color <- colors

  multipliers <-
    matrix(sample(c(-2,-1,1,2, -2.49, -1.49, -0.49, 1.49, 2.49),
                  2*5, replace = TRUE),
           ncol=2) |> unique()
  targets <- cbind(v1, v2) %*% t(multipliers) |> t()
  keepers <- !duplicated(round(targets))
  targs <- targets[keepers,] |>
    tibble::as_tibble()
  #targs <- targs |> LSTbook::take_sample(n=npts)
  if (!is.null(keep)) targs <- targs[keep, ]
  targs$labels <- labs[1:nrow(targs)]

  P <- vecframe |>
    gf_segment(
      y0 + y1 ~ x0 + x1,
      linewidth = linewidth,
      color = ~ colors,
      arrow = arrow()
    ) |>
    gf_label(ty ~ tx, label = ~ label) |>
    # gf_point(V2 ~ V1, data = targs, inherit=FALSE) |>
    gf_text(V2 ~ V1, label = ~labels, data = targs, inherit=FALSE) |>
    gf_refine(scale_color_identity(), coord_fixed()) |>
    gf_theme(theme_void())

  if (is.null(keep)) report <- multipliers |> tibble::as_tibble()
  else report <- multipliers[keep,] |> tibble::as_tibble()
  names(report) <- vecnames
  report$target <- labs[1:nrow(report)]

  list(P = P,
       multipliers = report, seed = seed)
}
