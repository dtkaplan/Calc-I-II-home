romeo_dynamics <- function(a=1, b=1, c=0, d=0,
                           title="Romeo's Love Dynamics", traj=FALSE) {
  drj <- domain(r = -1:1, j = -1:1)
  love_labels <- c("hate", "dislike", "indifferent", "fond", "in love")
  P <- vectorfield_plot(r ~ a*r + b*j,
                        j ~ c*r + d*j,
                        drj,
                        a = a, b = b, c = c, d = d,
                        transform = {\(x) x^0.5}) |>
    gf_refine(
      coord_fixed(),
      scale_y_continuous(breaks = seq(-1,1,by=0.5),
                         labels = love_labels),
      scale_x_continuous(breaks = seq(-1,1,by=0.5),
                         labels = love_labels)) |>
    gf_labs(x = "Romeo's attitude", y = "Juliet's attitude",
            subtitle = title) |>
    gf_theme(
      axis.text.x=element_text(angle = -45, vjust = 0, hjust=.5),
      axis.text.y=element_text(angle = 45, vjust = 0.5, hjust=0.5)
    )
  if (traj) {

  }
  P
}

eigenvalues_to_matrix <- function(mid, pm, signature = c(1, 1, 0, -1), as.matrix = TRUE) {
  sum <- 2*mid
  if (abs(pm) != abs(Re(pm))) { # It's complex!
    prod <- Re(mid^2 + pm * Conj(pm))
  } else {
    prod <- (mid + pm)*(mid - pm)
  }

  if (prod == 0) {
    warning("Eigenvalues are both zero. Matrix is zero.")
    res <- c(0, 0, 0, 0)
    if (as.matrix) res <- matrix(res, ncol=0)
    return(res)
  }
  a <- signature[1]; b <- signature[2];
  c <- signature[3]; d <- signature[4]
  # find a and d
  if (mid == 0) {
    a <- - d # possibly changes signature
  } else {
    if ((a + d) == 0) {
      a <- a + 0.1 + sign(a)
      d <- d + 0.1 + sign(d)
    }
    d <- sum*sign(d)/2
    a <- sum - d
  }
  # find b and c
  bc <- a*d - prod
  if (bc == 0) {
    if (sign(c) == 0 || sign(b) == 0) next
    else (c <- 0)
  } else {
    if (sign(b)*sign(c) == 0) {
      # must change signature
      c <- c + 0.1 + sign(c)
      b <- b + 0.1 + sign(b)
    }
    c <- sign(c) * sqrt(abs(bc))
    b <- sign(b) * sqrt(abs(bc))
    if (sign(bc) != sign(b*c)) c <- -c
  }

  res <- c(a, b, c, d)
  if (as.matrix) res <- matrix(res, nrow=2)

  res
}

romeo_eigen <- function(mid, pm, title = "", c = 1) {
  abcd <- eigenvalues_to_matrix(mid, pm, signature = c(1,1,c,1))
  romeo_dynamics(a=abcd[1], b=abcd[2], c=abcd[3], d=abcd[4], title = title)
}
