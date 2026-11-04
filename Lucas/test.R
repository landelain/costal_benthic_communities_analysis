
library(ggplot2)
library(dplyr)
library(tidyr)
library(datasets)
library(ggpubr)

data("iris")

iris2 <- mutate(iris, 
       sepalratio = Sepal.Length / Sepal.Width,
       petalratio = Petal.Length / Petal.Width)

iris
test$sepalratio

plot1 <- ggplot(data = iris) +
  geom_point(mapping = aes(x = Sepal.Width, y = Sepal.Length, color=Species)) + 
  geom_smooth(mapping = aes(x = Sepal.Width, y = Sepal.Length, color=Species))

plot2 = ggplot(data = iris) +
  geom_point(mapping = aes(x = sepalratio, y = petalratio, color=Species)) + 
  geom_smooth(mapping = aes(x = sepalratio, y = petalratio, color=Species))


ggarrange(plot1, plot2, ncol=2, nrow=1)



