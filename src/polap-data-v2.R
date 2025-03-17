library(ggplot2)
library(readr)

x <- read.delim("table1.tsv")

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
