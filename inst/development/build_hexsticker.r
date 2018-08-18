
library(hexSticker)
library(ggplot2)
library(ggimage)
p <- ggplot(aes(x = mpg, y = wt), data = mtcars) + geom_point()
p <- p + theme_void() + theme_transparent()
sticker(p, package="rbms", p_size=25, s_x=1, s_y=.75, s_width=1, s_height=1,
        h_fill="#0d51ce", h_color="#740ece", filename="docs/hexrbmslogo.png")
