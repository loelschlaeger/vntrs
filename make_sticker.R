output_directory <- file.path("man", "figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

angle <- seq(0, 2 * pi, length.out = 240)
ring_radii <- c(0.36, 0.66, 0.96)
neighborhoods <- do.call(
  rbind,
  lapply(seq_along(ring_radii), function(index) {
    radius <- ring_radii[index]
    data.frame(
      x = 0.95 * radius * cos(angle),
      y = 0.62 * radius * sin(angle),
      neighborhood = index
    )
  })
)

search_path <- data.frame(
  x = c(-0.84, -0.48, 0.30, 0.08),
  y = c(-0.43, 0.23, -0.18, 0.02)
)

algorithm_motif <- ggplot2::ggplot() +
  ggplot2::geom_path(
    data = neighborhoods,
    ggplot2::aes(x = x, y = y, group = neighborhood),
    color = "#67E8F9",
    linewidth = 0.9,
    alpha = 0.75
  ) +
  ggplot2::geom_path(
    data = search_path,
    ggplot2::aes(x = x, y = y),
    color = "#F8FAFC",
    linewidth = 1.2,
    lineend = "round",
    arrow = grid::arrow(
      length = grid::unit(0.08, "inches"),
      type = "closed"
    )
  ) +
  ggplot2::geom_point(
    data = search_path[-nrow(search_path), ],
    ggplot2::aes(x = x, y = y),
    shape = 21,
    size = 3.2,
    stroke = 0.8,
    color = "#F8FAFC",
    fill = "#F59E0B"
  ) +
  ggplot2::geom_point(
    data = search_path[nrow(search_path), , drop = FALSE],
    ggplot2::aes(x = x, y = y),
    shape = 23,
    size = 4.2,
    stroke = 0.9,
    color = "#F8FAFC",
    fill = "#FB7185"
  ) +
  ggplot2::coord_fixed(
    xlim = c(-1.05, 1.05),
    ylim = c(-0.72, 0.72),
    clip = "off"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
    panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

grDevices::pdf(NULL)
null_device <- grDevices::dev.cur()

hexSticker::sticker(
  subplot = algorithm_motif,
  package = "vntrs",
  filename = file.path(output_directory, "logo.png"),
  s_x = 1,
  s_y = 0.76,
  s_width = 1.35,
  s_height = 0.94,
  p_x = 1,
  p_y = 1.48,
  p_color = "#F8FAFC",
  p_family = "sans",
  p_fontface = "bold",
  p_size = 44,
  h_fill = "#0B132B",
  h_color = "#5EEAD4",
  h_size = 1.5,
  url = "",
  white_around_sticker = FALSE,
  dpi = 600
)

invisible(grDevices::dev.off(null_device))
