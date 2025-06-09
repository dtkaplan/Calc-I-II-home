# Defining functions





# Drawing vectors and matrices

draw_matrix <- function(M) {
  image(1:ncol(M), 1:nrow(M), t(M[nrow(M):1, ]),
        col = colorspace::diverging_hcl(
          15, h = c(180, 50), c = 80,
          l = c(20, 95), power = c(0.7, 1.3)),
        axes = FALSE,
        useRaster = TRUE)
}

somevecs <- tibble::tribble(
  ~ rootx, ~ rooty, ~ headx, ~ heady, ~ color, ~ name,
  0, 0, -2, 4, "blue", "u",
  1, -2, 5, -3, "green", "v",
  -1, -1, -4, -1, "orange", "w",
  1, 0, 4, 3, "brown", "x",
  0, -2, 3, -3, "salmon", "y",
  0, 0, 4, 3, "magenta", "b1",
  0, 0, -3, 2, "magenta", "b2",
) |>
  dplyr::mutate(labelx = (rootx + headx)/2,
         labely = (2*rooty + heady)/3)
solve_for <- function(vecnames) {
  somevecs |>
    dplyr::filter(name %in% vecnames) |>
    gf_segment(rooty + heady ~ rootx + headx,
               arrow = grid::arrow(length=unit(0.15, "inches"), type="closed"),
               color = ~ color, linewidth=2) |>
    gf_label(labely ~ labelx, label= ~ name, color = ~ color, size=3) |>
    gf_refine(scale_color_identity(),
              scale_y_continuous(limits=c(-5,5),
                                 breaks=(-5):5),
              scale_x_continuous(limits=c(-5,5),
                                 breaks=(-5):5),
              coord_fixed()) |>
    gf_labs(x="", y = "")
}


solve_graph <- function(seed=101, integers=FALSE) {
  set.seed(seed)
  if (integers) {
    u <- vec(sample(setdiff(-5:5, 0), size = 2))
    v <- vec(sample(setdiff(-5:5, 0), size = 2))
    # if too close in angle, try again
    while ( abs((u %dot% v) / (veclen(u)*veclen(v))) > 0.6 )
      v <- vec(sample(-5:5, size=2))
  } else {
    u <- vec(sample(c(-5:5, runif(3, -3, 3)), size = 2))
    v <- vec(sample(c(-5:5, runif(3, -3, 3)), size = 2))
    # if too close in angle, try again
    while ( abs((u %dot% v) / (veclen(u)*veclen(v))) > 0.6 )
      v <- vec(sample(c(-3:3, runif(3, -3, 3)), size=2))
  }
  ucoef <- sample(setdiff(-2:2, 0), size = 1)
  vcoef <- sample(setdiff(-2:2, 0), size = 1)
  b <- u*ucoef + v*vcoef
  both <- matrix( c(0, 0, u, 0, 0, v, 0, 0, b),
                  ncol=4, byrow = TRUE) |>
    tibble::as_tibble()
  names(both) <- c("rootx", "rooty", "headx", "heady")
  both$label <- c("a", "c", "b")
  both$color <- c(sample(c("blue", "darkorange", "tomato", "brown","magenta", "red"), size=2),
                  "forestgreen")
  both <- both |>
    dplyr::mutate(labelx = (rootx + headx)/2,
                  labely = (2*rooty + heady)/3)

  longest <- round(max(abs(c(u, v, b))) + 1)
  skip <- ifelse(longest > 7, 2, 1)

  both |>
    # dplyr::filter(label != "b") |>
    gf_segment(rooty + heady ~ rootx + headx,
               arrow = grid::arrow(length=unit(0.15, "inches"), type="closed"),
               color = ~ color, linewidth=2) |>
    gf_label(labely ~ labelx, label= ~ label, color = ~ color, size=3) |>
    gf_refine(scale_color_identity(),
              scale_y_continuous(limits=c(-longest, longest),
                                 breaks=seq(-longest, longest, by = skip)),
              scale_x_continuous(limits=c(-longest, longest),
                                 breaks=seq(-longest, longest, by = skip)),
              coord_fixed()) |>
    # gf_label(heady ~ headx,
    #          label = ~ label,
    #          color = ~ color,
    #          data = both |> dplyr::filter(label == "b")) |>
    gf_labs(x="", y = "") |>
    gf_refine(theme_minimal())
}



# Simple vector/matrix operations

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

# Convenience functions on vectors and matrices

# Give a colname to a vector or multiple colnames to a matrix
vec_names <- `colnames<-`
# Vector minus the mean
vec_mm <- function(...)  {
  v <- vec(...)
  v - mean(v, na.rm = TRUE)
}

# arbitrary vec->vec or vec->scalar function applied to each column of a vector
M_map <- function(M, f = "center") {
  choices <- c("center", "unit",
               "unitcenter", "mean", "var", "vlen")

  str <- as.character(substitute(f))
  if (str[1] %in% choices) {
    f <- switch (str[1],
                 center = \(x) (x - mean(x, na.rm = TRUE)),
                 unit = \(x) x / sqrt(sum(x^2)),
                 mean = \(x) mean(x, na.rm = TRUE),
                 var = \(x) var(x, na.rm = TRUE),
                 len = \(x) sqrt(sum(x*x)),
                 len2 = \(x) sum(x*x),
                 unitcenter = function(x) {
                   v <- x - mean(x, na.rm = TRUE)
                   v / sqrt(sum(v^2))
                 }
    )
  }
  if (!is.function(f)) stop("Not a recognized function.")
  first <- f(M[,1])
  res <- matrix(0, length(first), ncol(M))
  for (k in 1:ncol(M)) {
    res[, k] <- f(M[,k])
  }
  res
}



set_col_names <- function(v, nms) {
  colnames(v) <- nms
  v
}
one_sided <- function(tilde) {
 if (length(tilde) == 3)
   stop("The tilde expression should be one-sided, e.g. ~ b, not a ~ b")
}
# create a model matrix pipe style
data_M <- function(.data, tilde) {
  M <- model.matrix(tilde, data = .data |> tibble::remove_rownames())
  M[ , - which(colnames(M) == "(Intercept)"), drop = FALSE]
}

veclen <- function(v) sqrt(v %dot% v)
unitvec <- function(vec) vec/length(vec)
cang <- function(v, w) {(v %dot% w) / (veclen(v) * veclen(w))}
makeQ <- function(...) {
  dots <- list(...) |> bind_cols()
  qr.Q(qr(dots))
}

check_for_svd <- function(S) {
  if (inherits(S, "matrix")) S <- svd(S)
  if (!is.list(S) || !all(c("d","v", "u") %in% names(S)))
    stop("Argument must be matrix or the SVD of a matrix.")
  S
}

rank_mat <- function(S, thresh = 0.01) {
  S <- check_for_svd(S)
  sum(S$d > S$d[1] * thresh)
}

rand_mat <- function(nrow=3, ncol=6, rank = pmin(nrow, ncol)) {
  M <- matrix(runif(nrow*ncol), nrow = nrow, ncol = ncol)
  if (rank >= nrow || rank >= ncol) {
   M
  } else {
    approx_mat(S, n=sample(pmin(rank, nrow, ncol)))
  }
}

# Grab rank 1 matrix from SVD
approx_mat <- function(S, n=1, order = 0) { # input a matrix or the SVD of a matrix
  S <- check_for_svd(S)
  if (order > 0) {
    inds1 <- order(c(S$u[,order]))
    inds2 <- order(c(S$v[,order]))
  } else {
    inds1 <- 1:nrow(S$u)
    inds2 <- 1:nrow(S$v)
  }
  partial <- 0 # initial value
  for (k in n) {
    partial <- partial +  S$d[k] * S$u[inds1,k, drop = FALSE] %*% t(S$v[inds2,k, drop = FALSE])
  }

  partial
}

pretty_mat <- function(M) {
  inds1 <- order(rowSums(M))
  inds2 <- order(colSums(M))
  M[inds1, inds2]
}

#' Generate a matrix whose elements are selected randomly from a set
#' of specified values.
#' @param nrow number of rows for the matrix produced
#' @param ncol number of columns
#' @param values the set from which to draw (randomly) the
#' values in the matrix. Default: integers -9 to 9
values_mat <- function(nrow=4, ncol=3, values = -9:9) {
  matrix(sample(values, size = nrow * ncol, replace = TRUE),
         nrow = nrow, ncol=ncol)
}

# Typeset matrices in LaTeX (for assignments and such)
Lmat <- function(nr, nc) {
  values_mat(nr, nc) |> latex_helper()
}
latex_helper <- function(matr) {
  printmrow <- function(x) {

    cat(cat(x,sep=" & "),"\\\\ \n")
  }

  cat("\\left(\\begin{array}{r}","\n")
  body <- apply(matr,1,printmrow)
  cat("\\end{array}\\right)")
}

# Fourier transform stuff
# sigfft shows only the left-hand side of the FFT
# isigfft restores it by adding back the right-hand side then inverse FFT.

ifft <- function(x) fft(x/length(x), inverse = TRUE)

sigfft <- function(x) {
  tmp <- fft(x)
  nyquist <- round((length(x) + 1.1 - (length(x) %% 2))/2)
  tmp[1:nyquist]
}
isigfft <- function(Lfftx) {
  tmp <- if (length(Lfftx) %% 2) {
    # odd length
    c(Lfftx, Conj(rev(Lfftx[c(-1, -length(Lfftx))])))

  } else {
    c(Lfftx, Conj(rev(Lfftx[-1])))
  }
  ifft(tmp) |> Re()
}

squash_small <- function(x, tol=1e-9) {
  r1 <- Re(x)
  r1 <- ifelse(abs(r1) < tol, 0, r1)
  c1 <- Im(x)
  c1 <- ifelse(abs(c1) < tol, 0, c1)
  # return real part if all of imaginary part is small
  if (all(c1 == 0)) r1
  else r1 + 1i * c1
}

# Plot the amplitude spectrum versus frequency
sig_amp_spec <- function(x, sampfreq=100) {
  Tmp <- tibble::tibble(
    amp = abs(sigfft(x)),
    frequency = seq(0, sampfreq/2, length = length(amp)))
  Tmp |>
    gf_segment(0 + amp ~ frequency + frequency, alpha = 0.2) |>
    gf_point(amp ~ frequency, color = "blue", size = 0.5)
}



# Plotting a function with room for drawing the anti-deriv
drawFpair <- function(f, dom = domain(x = 0:4), bottom = -0.5, alpha=0) {
  # alpha = 0 for problem
  # alpha = 1 for answer
  fname <- as.character(substitute(f))
  vname <- names(dom)
  tilde <- glue::glue("{fname}({vname}) ~ {vname}")
  label <- glue::glue("{fname}({vname})")
  Ftilde <- glue::glue("{toupper(fname)}({vname}) ~ {vname}")
  labelF <- glue::glue("{toupper(fname)}({vname})")
  P1 <- slice_plot(as.formula(tilde), dom, npts = 1000) |>
    gf_labs(y = label, x = "") |>
    gf_lims(y = c(bottom, NA)) |>
    gf_hline(yintercept = 0, color = "blue", linetype = "dashed") |>
    gf_theme(theme_minimal(base_size = 22))
  assign(toupper(fname), antiD(as.formula(tilde)))
  P2 <- slice_plot(as.formula(Ftilde), dom, alpha = alpha) |>
    gf_labs(y = labelF, x = vname) |>
    gf_hline(yintercept = 0, color = "blue", linetype = "dashed") |>
    gf_theme(theme_minimal(base_size = 22))
  list(P1 = P1, P2 = P2)
}

# Draw Riemann boxes on a function:

rplot <- function(tilde, domain, npts = 11, position=c("left", "middle", "right"), ...) {
  position = match.arg(position)
  f <- makeFun(tilde, ...)
  ends <- domain[[1]]
  h <- diff(ends)/(npts - 1)
  points <- seq(ends[1], ends[2], by = h)
  Pts <- tibble::tibble(
    left = points[-length(points)],
    right = points[-1],
    mid = (left + right)/2
  )
  Pts$val <- if (position == "left") f(Pts$left)
  else if (position == "middle") f(Pts$mid)
  else if (position == "right") f(Pts$right)
  else stop("Invalid <position> argument.")

  Pts$color <- ifelse(Pts$val > 0, "blue", "tomato")
  slice_plot(tilde, domain, npts=501) |>
    gf_hline(yintercept = ~ 0, color = "blue") |>
    gf_rect(0 + val ~ left + right, data = Pts,
            inherit = FALSE, fill = NA, alpha = 0,
            color = "black", linewidth = 0.1) |>
    gf_rect(0 + val ~  left + right, data = Pts,
            fill = ~ color,
            inherit = FALSE, alpha = 0.2) |>
    gf_refine(scale_fill_identity()) |>
    gf_theme(theme_void)

}


if (!exists("take_sample")) {
  take_sample <<- function (x, n, replace = FALSE, ...) {
    UseMethod('take_sample')
  }


  #' @export
  take_sample.vector <<- function(x, n=length(x), replace=FALSE, ...) {
    base::sample(x, size = n, replace = replace, ...)
  }

  #' @export
  take_sample.data.frame <<- function(x, n = nrow(x), replace = FALSE, ..., .by = NULL) {

    # slice_sample uses `by` instead of `.by`
    # I can get this to work only by turning `.by` into a character string
    # containing the desired names.
    groups <- substitute(.by) # `groups` will be a call

    if (!is.null(groups)) {
      # handle cases with multiple .by variables
      if (is.call(groups)) {
        if (!is.name(groups[[2]])) groups <- eval(groups)
        else groups <- all.vars(groups)
      }
    }

    dplyr::slice_sample(x, n = n, by = all_of(groups), replace = replace)
  }

  #' @export
  take_sample.datasim <<- function(x, n = 5, replace = FALSE, ..., seed = NULL, report_hidden=FALSE) {
    datasim_run(x, n = n, seed = seed, report_hidden = report_hidden)
  }

  #' @param .by Variables to use to define groups for sampling, as in `{dplyr}`. The sample size
  #' applies to each group.
  #' @param groups Variable indicating blocks to sample within
  #' @param orig.ids Logical. If `TRUE`, append a column named "orig.ids" with the
  #' row from the original `x` that the same came from.
  #' @param prob Probabilities to use for sampling, one for each element of `x`
  #'
  #' @rdname take_sample
  #' @export
  take_sample.default <<- function(x, n = length(x), replace=FALSE, prob=NULL, .by = NULL,
                                  groups = .by, orig.ids=FALSE, ...) {
    size <- n
    missingSize <- missing(size)
    haveGroups <- ! is.null(groups)
    if (length(x) == 1L && is.numeric(x) && x >= 1) {
      n <- x
      x <- 1:n
      if (missingSize)  size <- n
    } else {
      n <- length(x)
      if (missingSize) size <- length(x)
    }
    if (haveGroups && size != n) {
      warning("'size' is ignored when using groups.")
      size <- n
    }
    ids <- 1:n

    if (haveGroups) {
      groups <- rep( groups, length.out = size)  # recycle as needed
      result <- stats::aggregate( ids, by = list(groups), FUN = base::sample,
                                  simplify = FALSE,
                                  replace = replace, prob = prob)
      result <- unlist(result$x)
      if (orig.ids) { nms <- ids[result] }
      result <- x[result]
      if (orig.ids) { names(result) <- nms }
      return(result)
    }
    result <- base::sample(x, size, replace = replace, prob = prob)
    return(result)
  }

}

if (!exists("resample")) {
  resample <<- function(..., replace = TRUE) {
    take_sample(..., replace = replace)
  }
}

  #' A convenience function for shuffling, typically used with
  #' within model_train(), but available elsewhere for, e.g. demonstrations
  #' @rdname take_sample
  #' @export
if (!exists("shuffle")) {
  shuffle <<- function(...) {
    take_sample(...)
  }
}



## Draw a flow field and nullclines
draw_flow <- function(seed=1996, center = c(0,0), width=5,
                      dx = 0.1, dy=0.1, ngrid=31, arrow=1) {
  xrange <- center[1] + width*c(-1,1)
  yrange <- center[2] + width*c(-1,1)
  Grid <- expand.grid(x=seq(xrange[1], xrange[2], length=ngrid),
                      y=seq(yrange[1], yrange[2],  length=ngrid))
  # This one is always the same, so zooming in will work
  Grid2 <- expand.grid(x = seq(-7, 7), y = seq(-7, 7))
  dom <- domain(x=!!xrange, y=!!yrange)
  xraw <- doodle_fun(~ x + y, seed = seed)
  yraw <- doodle_fun(~ x + y, seed = seed + 1)
  Stats <- Grid2 |>
    mutate(xvals = xraw(x,y),
           yvals = yraw(x,y)) |>
    summarize(midx = mean(xvals), midy = mean(yvals),
              sdx = sd(xvals), sdy = sd(yvals))

  xfun <- function(x, y) x + dx*(xraw(x,y)-Stats$midx)/Stats$sdx
  yfun <- function(x, y) y + dy*(yraw(x,y)-Stats$midy)/Stats$sdy
  Arrows <- Grid |>
    mutate(xend = xfun(x,y), yend=yfun(x,y))
  contour_plot(xfun(x,y) - x ~ x + y, domain=dom,
               contour_color="dodgerblue", contours_at = c(0),
               skip=0, labels=FALSE) |>
    contour_plot(yfun(x,y) - y ~ x + y,
                 contour_color="orange3", contours_at = c(0),
                 skip=0, filled=FALSE, labels=FALSE) |>
    gf_segment(y + yend ~ x + xend, data = Arrows, size=arrow*0.5) |>
    gf_point(yend ~ xend, data = Arrows, size = arrow, alpha=0.3) |>
    gf_refine(coord_fixed())
}

## For formatting

make_numbering <- function(word) {
  num <- 1
  function() {
    res <- glue::glue(word)
    num <<- num + 1
    return(res)
  }
}
Q <- make_numbering("[**Question 7.{num}**]{{.underline}}")
vspace <- function(n=2){
  if (knitr::is_latex_output()) {
    paste(rep("\n\n \n\n ", n), collapse = "\n")
  } else ""
}
