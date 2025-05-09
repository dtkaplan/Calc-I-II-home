draw_matrix <- function(M) {
  image(1:ncol(M), 1:nrow(M), t(M),
        col = colorspace::diverging_hcl(
          15, h = c(180, 50), c = 80,
          l = c(20, 95), power = c(0.7, 1.3)),
        axes = FALSE)
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
  mutate(labelx = (rootx + headx)/2,
         labely = (2*rooty + heady)/3)
solve_for <- function(vecnames) {
  somevecs |>
    filter(name %in% vecnames) |>
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
