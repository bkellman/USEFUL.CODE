# data = numeric_matrix_{variables (or genes) x samples}
# annotation = data_frame_{ samples x annotation_variables}

# multiple lables per recursion, not yet implimented
# batch.par list of which annotation to remove e.g. list( 'location', c('location','weight') ) 
#			if batch.par=NULL then the top batch effect is automatically chosen

# batch.par list of which annotation to remove e.g. c( 'location','weight' ) 
#			if batch.par=NULL then the top batch effect is automatically chosen

# QC.par = list( outcex = size of outlier dots , las = orientation of labels: verticle (2) or horizontal (1), scaleV = boolean to scale variables)


specialSave<-function(x,folder,annotation=NULL,fileName=NULL,
	heatmap=F,cor.heatmap=F,cor.method=c("pearson", "kendall", "spearman","all"),
	batch.anyl=F,batch.correct=F,batch.par=NULL,batch.rec=1,dontCorrect=NULL,
	QC=F,QC.par=NULL ){

	#collect outfiles
	outfiles = c()

	print(batch.rec)
	#recursive termination
	if(batch.rec<=0){
		batch.correct=F
	}

	#checks
	annotation_data_row_match(annotation,x)
	cor_method(cor.method)
	batch_annotation(annotation,batch.anyl,batch.correct)

	# establish save directory
	dir.create( folder, showWarnings = FALSE)

	# batch analysis and correction
	if(!is.null(annotation)){
		if (!identical(rownames(annotation), colnames(x))) {
        	warning("Colnames of x are not the same as rownames of annotation")
    	}
		library(swamp)
		# establish save directory
		dir.create( file.path(folder,'batch') , showWarnings = FALSE)
		if(batch.anyl){
			prince <- run.batch.anyl(x,annotation,folder)
		}
		if(batch.correct){
			outfiles<-c(outfiles, 
				run.batch.correct(x,folder,annotation,fileName, 
				heatmap,cor.heatmap,cor.method, batch.anyl,batch.correct,batch.par,batch.rec,dontCorrect, QC,QC.par,prince) )
		}
	}

	#x = t(x) #### TRANSPOSE

	# Quality control
	if(QC){
		run.QC(t(x),folder,QC.par)
	}

	#heatmaps
	if(heatmap | cor.heatmap){
		library(gplots)
	}
	if(heatmap){
		run.heatmap(x,folder)
	}
	if(cor.heatmap){
		run.cor.heatmap(x,folder,cor.method)
	}

	#save
	if(is.null(fileName)){
		save(x,annotation,file=file.path(folder,'data.annotation.rda'))
		write.csv(x,file=file.path(folder,'data.csv'))
		write.csv(annotation,file=file.path(folder,'annotation.csv'))
		return(c(file.path(folder,'data.annotation.rda'),outfiles))
	}else{
		save(x,annotation,file=file.path(folder,paste(fileName,'.data.annotation.rda',sep='')) ) 
		write.csv(x,file=file.path(folder,paste(fileName,'.data.rda',sep='')) )		
		write.csv(annotation,file=file.path(folder,paste(fileName,'.annotation.rda',sep='')) )	
		return(c(file.path(folder,paste(fileName,'.data.annotation.rda',sep='')),outfiles))	
	}
}

########
# batch::swamp
########

run.batch.anyl <- function(g,o,folder){
	if(ncol(g)<50 | nrow(g)){
		if(ncol(g)<10){
			topper=ncol(g)
		}else{
			topper=10
		}
		show.topper=min(ncol(g),nrow(g))
	}else{
		topper=10
		show.topper=50
	}
	res1<-prince(g,o,top=topper,permute=TRUE)

	pdf(file.path(folder,'batch/batch_PCA.pdf'))
	prince.plot(prince=res1)
	dev.off()
	# to plot p values of linear models: lm(principal components ~ sapmle annotations).
	# to see if the variation in the data is associated with sample annotations.

	pdf(file.path(folder,'batch/batch_PCA_var.pdf'))
	res2<-prince.var.plot(g,show.top=show.topper,npermute=10)
	dev.off()
	# to see how many principal components carry informative variation
	## hierarchical clustering analysis

	pdf(file.path(folder,'batch/batch_hca.pdf'))
	hca.plot(g,o)
	dev.off()
	# to show a dendrogram with sample annotations below

	res3<-hca.test(g,o,dendcut=2,test="chi")
	print(paste('hca.test:',res3,sep='\n'))

	## associations between sample annotations
	pdf(file.path(folder,'batch/batch_confounding.pdf'))
	res4<-confounding(o,method="chi")
	dev.off()

	return(res1)
}

run.batch.correct <- function(x,folder,annotation,fileName, heatmap,cor.heatmap,cor.method, batch.anyl,batch.correct,batch.par,batch.rec, dontCorrect, QC,QC.par,prince){

#	batch.correct=F
	batch.par.prev=batch.par
	batch.par=NULL # only the1st batch correction is regulated
	batch.rec = batch.rec-1

	outfiles=c()
	g = x
	o = annotation

	# make for loop for killPC(c(n))
	if(!('pc1') %in% dontCorrect){
		gadj3<-kill.pc(g,pc=1)
		dir.create( file.path(folder,'batch','kill.pc_1') , showWarnings = FALSE)
		print('batch_killpc1')
		dontCorrect1 = c(dontCorrect,'pc12','pc1')
		outfiles<-c(outfiles,
			specialSave(gadj3,file.path(folder,'batch','kill.pc_1'),annotation,fileName,
			heatmap,cor.heatmap,cor.method, batch.anyl,batch.correct,batch.par,batch.rec,dontCorrect1, QC,QC.par))
	}

	if(!('pc12' %in% dontCorrect)){
		gadj4<-kill.pc(g,pc=1:2)
		dir.create( file.path(folder,'batch','kill.pc_12') , showWarnings = FALSE)
		print('batch_killpc12')
		dontCorrect12 = c(dontCorrect,'pc12','pc1')
		outfiles<-c(outfiles,
			specialSave(gadj4,file.path(folder,'batch','kill.pc_12'),annotation,fileName,
			heatmap,cor.heatmap,cor.method, batch.anyl,batch.correct,batch.par,batch.rec,dontCorrect12, QC,QC.par))
	}

	if(is.null(batch.par.prev)){
		lnrsq = prince$rsquared[,1]  # get r^2 and p(F(linear regression)) from 1st PC
		lnp   = prince$linp[,1]
		dont_i= which(names(lnrsq) %in% dontCorrect)
		dont_j= which(names(lnp) %in% dontCorrect)
		if(sum(dont_i==dont_j)!=length(dont_i) | sum(dont_i==dont_j)!=length(dont_j)){stop('dont_i!=dont_j')}
		lnrsq = lnrsq[-dont_i]
		lnp   = lnp[-dont_j]
		batch.indx = c( which(colnames(o) == names(which.max(lnrsq))),
						which(colnames(o) == names(which.min(lnp    ))) )
	}else{
	#	print('! null batchpar')
		batch.indx = which(colnames(o) %in% batch.par.prev)
	}
	print(batch.indx)
	batch.indx = unique(batch.indx)

	for(bp in batch.indx){
		lin1<-adjust.linearmodel(g,o[,bp])
		fn.tmp = paste('lin_',colnames(o)[bp],sep='')
		print(paste('batch_lin',fn.tmp))
		dir.create( file.path(folder,'batch',fn.tmp) , showWarnings = FALSE)
		dontCorrect_batch = c(dontCorrect,colnames(o)[bp])
		outfiles<-c(outfiles,
			specialSave(lin1,file.path(folder,'batch',fn.tmp),annotation,fileName,
			heatmap,cor.heatmap,cor.method, batch.anyl,batch.correct,batch.par,batch.rec, dontCorrect_batch,QC,QC.par))
	}

	return(outfiles)
}

########
# QC
########
# QC.par = list( outcex = size of outlier dots , las = orientation of labels: verticle (2) or horizontal (1), scaleV = boolean to scale variables)
# x = samples x genes/variables
run.QC <- function(x,folder,QC.par){  # checked
	if('outcex' %in% names(QC.par) ){
		outcex = QC.par$outcex
	}else{
		outcex = 1
	}
	if('las' %in% names(QC.par) ){
		las = QC.par$las
	}else{
		las = 2
	}
	xVar = x
	if('scaleV' %in% names(QC.par) ){
		if(QC.par$scaleV){
			xVar = scale(x) #scale columns
		}
	} # else don't scale
	if('varcex' %in% names(QC.par) ){
		varcex = QC.par$varcex
	}else{
		varcex = 1
	}
	if('samcex' %in% names(QC.par) ){
		samcex = QC.par$samcex
	}else{
		samcex = 1
	}
	
	l_col =  split(xVar,rep(1:ncol(x), each = nrow(x))) # use xVar, which is only scaled if QC.par$scaleV
	l_row =  split(t(x),rep(1:ncol(t(x)), each = nrow(t(x))))
	if(!is.null(colnames(x))) { names(l_col) = colnames(x) }
	if(!is.null(rownames(x))) { names(l_row) = rownames(x) }
	# box and whiskers
	pdf(file.path(folder,'boxplot.variables.pdf'))
	boxplot(l_col,main='Boxplot of Variables',las=las,outcex=outcex,cex.axis=varcex); dev.off()
	pdf(file.path(folder,'boxplot.sample.pdf'))
	boxplot(l_row,main='Boxplot of Samples',las=las,outcex=outcex,cex.axis=samcex)  ; dev.off()
	# density
	pdf(file.path(folder,'density.variables.gaussian.pdf'))
	multidense(l_col,main='Density of Variables: Gaussian Kernel',xlab=NULL,legend=NULL,x.leg=NULL,y.leg=NULL,kernel='gaussian')    ; dev.off()
	pdf(file.path(folder,'density.variables.rectangular.pdf'))
	multidense(l_col,main='Density of Variables: Rectangular Kernel',xlab=NULL,legend=NULL,x.leg=NULL,y.leg=NULL,kernel='rectangular'); dev.off()
	pdf(file.path(folder,'density.samples.gaussian.pdf'))
	multidense(l_row,main='Density of Samples: Gaussian Kernel',xlab=NULL,legend=NULL,x.leg=NULL,y.leg=NULL,kernel='gaussian')      ; dev.off()
	pdf(file.path(folder,'density.samples.rectangular.pdf'))
	multidense(l_row,main='Density of Samples: Rectangular Kernel',xlab=NULL,legend=NULL,x.leg=NULL,y.leg=NULL,kernel='rectangular')  ; dev.off()
}

########
# heatmaps
########
run.heatmap <- function(x,folder){  # checked
	pdf(file.path(folder,'data.heatmap.pdf'),height=10,width=10)
	heatmap.2(t(x),trace='none',cexCol=.5,cexRow=.5)
	dev.off()
}

run.cor.heatmap <- function(x,folder,cor.method='all'){  #checked
	if(cor.method=='all'){
		cor.method = c('kendall','spearman','pearson')
	}
	if( 'kendall' %in% cor.method ){
		pdf(file.path(folder,'cor.heatmap.kendall.pdf'),height=10,width=10)
		cr = cor(t(x),method='kendall')
		heatmap.2(cr,trace='none',cexCol=.5,cexRow=.5)
		dev.off()
	}
	if( 'pearson' %in% cor.method ){
		pdf(file.path(folder,'cor.heatmap.pearson.pdf'),height=10,width=10)
		cr = cor(t(x),method='pearson')
		heatmap.2(cr,trace='none',cexCol=.5,cexRow=.5)
		dev.off()
	}
	if( 'spearman' %in% cor.method ){
		pdf(file.path(folder,'cor.heatmap.spearman.pdf'),height=10,width=10)
		cr = cor(t(x),method='pearson')
		heatmap.2(cr,trace='none',cexCol=.5,cexRow=.5)
		dev.off()
	}
}
########
# checks
########
batch_annotation <- function(annotation,batch.anyl,batch.correct){
	if(is.null(annotation)){
		if(batch.anyl | batch.correct){
			stop('batch correction requires annotation')
		}
	}
}
annotation_data_row_match <- function(a,d){
	if(!is.null(a)){ 
		if( nrow(a) != ncol(d) ){
			stop('Data (x) columns and annotation rows match.  Both correspond to the samples')
		}
	}
}
cor_method <- function(par){
	return( par %in% c('kendall','spearman','pearson') | par == 'all')
}

### test
test<-function(){
	x = matrix(rpois(110,lambda=5),10,11,byrow=F)
	x[,1:5] = x[,1:5]*10*runif(5)
	colnames(x) = 1:11

	ann = data.frame(ab=c(1,1,1,1,1,1,0,0,0,0,0),bc=runif(11))

	specialSave(x=x,folder='tmp',annotation=ann,fileName=NULL,
		heatmap=T,cor.heatmap=T,cor.method="all",
		batch.anyl=T,batch.correct=T,batch.par='bc',batch.rec=1,
		QC=T,QC.par=NULL)
	specialSave(x=x,folder='tmp',annotation=ann,fileName=NULL,
		heatmap=T,cor.heatmap=T,cor.method="all",
		batch.anyl=T,batch.correct=T,batch.par=NULL,batch.rec=2,
		QC=T,QC.par=NULL)
}
####