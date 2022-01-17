library(lme4)
# install.packages("lmerTest")

library(lmerTest)
library(MASS)
install.packages("ggplot2")
library(ggplot)

# data <- read.csv('Z:/home/ryanh/projects/ripple_heterogeneity/assembly_unit_corrs/results/df.csv')
# summary(data)


# data$basepath = factor(data$basepath)
# data$name = factor(data$name)
# data$membership = factor(data$membership)


# p<-ggplot(data, aes(x=name, y=rho, color=membership)) +
#     geom_boxplot() +
#     geom_jitter(position=position_dodge(1),alpha=.5)


data("mtcars")
cars <- mtcars
print(head(cars,10))

x <- c(0:10)
plot(x,sin(x))
