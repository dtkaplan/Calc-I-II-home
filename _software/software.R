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
