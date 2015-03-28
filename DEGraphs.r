###################################################
### code chunk number 2: morelib
###################################################
library("R.utils")
##library(graph)
##library(rrcov) ## for 'T2.test'
library(corpcor)
library(KEGGgraph)
library(Rgraphviz)
##library(RBGL)
library(fields) # For image.plot called in plotValuedGraph
library(lattice)
library(marray)
verbose <- TRUE

###################################################
### code chunk number 1: lib
###################################################
library(DEGraph)



###################################################
### code chunk number 3: data
###################################################
data("Loi2008_DEGraphVignette", package="DEGraph")
classData <- classLoi2008
exprData <- exprLoi2008
annData <- annLoi2008
grList <- grListKEGG




DEGraph.wrapper<-function(classData,exprData,annData,grList ){

### R code from vignette source 'vignettes/DEGraph/inst/doc/DEGraph.Rnw'

###################################################
### code chunk number 4: ttests
###################################################
## Individual t-test p-values
X1 <- t(exprData[, classData==0])
X2 <- t(exprData[, classData==1])
ttpv <- c()
tts <- c()
for(i in 1:ncol(X1)) {
  tt <- t.test(X1[,i],X2[,i])
  ttpv[i]=unlist(tt$p.value)
  tts[i]=unlist(tt$statistic)
}
names(ttpv) <- names(tts) <- rownames(exprData)


###################################################
### code chunk number 5: mvtests
###################################################
prop <- 0.2
## Multivariate tests
resList <- NULL
for (ii in seq(along=grList)) {
  gr <- grList[[ii]]
  res <- testOneGraph(gr, exprData, classData, verbose=verbose, prop=prop)
  resList <- c(resList, list(res))
}
resNames <- names(grList)
pLabels <- sapply(grList, attr, "label")

## get rid of NULL results (no connected component of size > 1)
isNULL <- sapply(resList, is.null)
if (sum(isNULL)) {
  grList[isNULL]
  resList <- resList[!isNULL]
  resNames <- names(grList)[!isNULL]
  pLabels <- pLabels[!isNULL]
}

resL <- sapply(resList, length)
graphNames <- rep(resNames, times=resL)
pathwayNames <- rep(pLabels, times=resL)

graphList <- NULL
for (res in resList) {
  grl <- lapply(res, FUN=function(x) {
    x$graph
  })
  graphList <- c(graphList, as.list(grl))
}

ndims <- NULL
for (res in resList) {
  ndim <- sapply(res, FUN=function(x) {
    x$k
  })
  ndims <- c(ndims, ndim)
}

pKEGG <- NULL
for (res in resList) {
  pp <- sapply(res, FUN=function(x) {
    x$p.value
  })
  pKEGG <- cbind(pKEGG, pp)
}
colnames(pKEGG) <- graphNames
rn <- rownames(pKEGG)
rownames(pKEGG)[grep("Fourier", rn)] <- paste("T2 (", round(100*prop), "% Fourier components)", sep="")

if (exists("maPalette", mode="function")) {
  pal <- maPalette(low="red", high="green", mid="black", k=100)
} else {
  pal <- heat.colors(100)
}
shift <- tts # Plot t-statistics
names(shift) <- translateGeneID2KEGGID(names(tts))
fSignif <- which(pKEGG[2,] < 0.05)
fSignif <- fSignif[order(pKEGG[2,fSignif])]


###################################################
### code chunk number 6: path1
###################################################
gIdx <- fSignif[1]
gr <- graphList[[gIdx]]
mm <- match(translateKEGGID2GeneID(nodes(gr)), rownames(annData))
dn <- annData[mm, "NCBI.gene.symbol"]
res <- plotValuedGraph(gr, values=shift, nodeLabels=dn, qMax=0.95, colorPalette=pal, height=40, lwd=1, cex=0.3, verbose=verbose)
stext(side=3, pos=0, pathwayNames[gIdx])
ps <- signif(pKEGG[, gIdx],2)
txt1 <- paste("p(T2)=", ps[1], sep="")
txt2 <- paste("p(T2F[", ndims[gIdx], "])=", ps[2], sep="")
txt <- paste(txt1, txt2, sep="\n")
stext(side=3, pos=1, txt)
image.plot(legend.only=TRUE, zlim=range(res$breaks), col=pal, legend.shrink=0.3, legend.width=0.8, legend.lab="t-scores", legend.mar=3.3) 


###################################################
### code chunk number 7: path2
###################################################
gIdx <- fSignif[5]
gr <- graphList[[gIdx]]
mm <- match(translateKEGGID2GeneID(nodes(gr)), rownames(annData))
dn <- annData[mm, "NCBI.gene.symbol"]
res <- plotValuedGraph(gr, values=shift, nodeLabels=dn, qMax=0.95, colorPalette=pal, height=40, lwd=1, cex=0.3, verbose=verbose)
stext(side=3, pos=0, pathwayNames[gIdx])
ps <- signif(pKEGG[, gIdx],2)
txt1 <- paste("p(T2)=", ps[1], sep="")
txt2 <- paste("p(T2F[", ndims[gIdx], "])=", ps[2], sep="")
txt <- paste(txt1, txt2, sep="\n")
stext(side=3, pos=1, txt)
image.plot(legend.only=TRUE, zlim=range(res$breaks), col=pal, legend.shrink=0.3, legend.width=0.8, legend.lab="t-scores", legend.mar=3.3) 


###################################################
### code chunk number 8: showpath1
###################################################
gIdx <- fSignif[1]
gr <- graphList[[gIdx]]
mm <- match(translateKEGGID2GeneID(nodes(gr)), rownames(annData))
dn <- annData[mm, "NCBI.gene.symbol"]
res <- plotValuedGraph(gr, values=shift, nodeLabels=dn, qMax=0.95, colorPalette=pal, height=40, lwd=1, cex=0.3, verbose=verbose)
stext(side=3, pos=0, pathwayNames[gIdx])
ps <- signif(pKEGG[, gIdx],2)
txt1 <- paste("p(T2)=", ps[1], sep="")
txt2 <- paste("p(T2F[", ndims[gIdx], "])=", ps[2], sep="")
txt <- paste(txt1, txt2, sep="\n")
stext(side=3, pos=1, txt)
image.plot(legend.only=TRUE, zlim=range(res$breaks), col=pal, legend.shrink=0.3, legend.width=0.8, legend.lab="t-scores", legend.mar=3.3) 


###################################################
### code chunk number 9: showpath2
###################################################
gIdx <- fSignif[5]
gr <- graphList[[gIdx]]
mm <- match(translateKEGGID2GeneID(nodes(gr)), rownames(annData))
dn <- annData[mm, "NCBI.gene.symbol"]
res <- plotValuedGraph(gr, values=shift, nodeLabels=dn, qMax=0.95, colorPalette=pal, height=40, lwd=1, cex=0.3, verbose=verbose)
stext(side=3, pos=0, pathwayNames[gIdx])
ps <- signif(pKEGG[, gIdx],2)
txt1 <- paste("p(T2)=", ps[1], sep="")
txt2 <- paste("p(T2F[", ndims[gIdx], "])=", ps[2], sep="")
txt <- paste(txt1, txt2, sep="\n")
stext(side=3, pos=1, txt)
image.plot(legend.only=TRUE, zlim=range(res$breaks), col=pal, legend.shrink=0.3, legend.width=0.8, legend.lab="t-scores", legend.mar=3.3) 



}
