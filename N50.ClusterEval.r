######
# N50 eval heatmap
#####
library(gplot)
optimalLeafOrder......

N50.eval.heatmap.2<-function(x,labels,optimalLeafOrder=(),...){
	h <- heatmap.2(x,dendrogram=optimalLeafOrder,...)
	return(N50.eval(h))
}

N50.eval <- function(hm){
	
}