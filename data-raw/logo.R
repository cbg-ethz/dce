library(ggplot2)
library(ggpixel)


set.seed(42)
ggpixel(
  "dce",
  noise.difference = 10,
  color.palette = "viridis",

  overlay = 6,
  overlay.radius = 13,

  margin.left = 5,
  margin.right = 5,
  margin.top = 10,
  margin.bottom = 10,

  letter.spacing = 1,
  aspect.ratio = 1
)
ggsave("man/figures/logo.png")

# needs a couple runs to actually crop. Why? No idea...
for (i in seq_len(5)) {
  knitr::plot_crop("man/figures/logo.png")
}
