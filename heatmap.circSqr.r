#  http://rgraphgallery.blogspot.com/2013/04/rg-heatmap-with-overlayed-circle-size.html

test<-function(){
set.seed (78888)
rectheat = sample(c(rnorm (10, 5,1), NA, NA), 150, replace = T)
circlefill =  rectheat*10 + rnorm (length (rectheat), 0, 3)
circlesize = rectheat*1.5 + rnorm (length (rectheat), 0, 3)
myd <- data.frame (rowv = letters[rep (1:10, 15)], columnv = rep(1:15, each = 10),
          rectheat, circlesize, circlefill)
}
  
require(ggplot2)
 pl1 <-  ggplot(myd, aes(y = factor(rowv),  x = factor(columnv))) +  geom_tile(aes(fill = rectheat)) +  scale_fill_continuous(low = "blue", high = "green")
 pl1  +      geom_point(aes(colour = circlefill,  size =circlesize))  +    scale_color_gradient(low = "yellow",   high = "red")+     scale_size(range = c(1, 20))+   theme_bw()