suppressPackageStartupMessages(library("optparse"))
library(ggplot2)
library(readr)

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--inum"),
  action = "store",
  help = "number",
  metavar = "<NUMEBR>"
)
args1 <- parse_args(parser)

x <- read.delim(paste("maintable1-", args1$inum, ".tsv", sep = ""))

ggplot(x, aes(x = Rate, y = Alpha)) +
  geom_point(color = "blue", alpha = 0.6) + # Scatter plot points
  geom_smooth(method = "lm", color = "red", fill = "pink", se = TRUE) + # Regression line with confidence interval
  labs(
    title = "Scatter Plot with Linear Regression",
    x = "Independent Variable (x)",
    y = "Dependent Variable (y)"
  ) +
  theme_minimal()

m1 <- lm(Alpha ~ Rate, data = x)
summary(m1)
