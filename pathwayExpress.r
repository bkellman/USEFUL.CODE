# http://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf

library(biomaRt)
library(KEGGREST)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

###########
# input
# load('Data/rna.rda')
# load('../Grad/Labs/Couchesne/ASDanyl/Data/rna.rda')
exp.mat = rna
collect.genes = function(x) sum(x)/sqrt(length(x))
collect.transcripts = median
###########

#kegg pathway -> genes <- transcripts
# get all genes in each pathway
p = keggLink("pathway","hsa")
map = split(names(p), unname(p))

hsa2ens = list()
ensg2enst = list()
pathExp = list()
ensgExp = list()
#iterate over pathways
for(m in 1:length(map)){ 
	path = names(map)[m]
	print(paste(m,'/',length(map),': ',path,sep=''))
	pathExp[[path]] = 0 			# build pathway expression list
	pathGeneExp = list()		# collect expression of genes in pathway
	#iterate over genes in pathways
	for(h in 1:length(map[[m]])){
		# convert hsa to gene id
		gene.hsa = map[[m]][h]
		if(!(gene.hsa %in% names(hsa2ens))){
			hsa2ens[[gene.hsa]] <- query.kegg.h2g(gene.hsa)
		}
		ensg <- hsa2ens[[gene.hsa]]
		# if gene expression has already been collected from the transcript level
		if(!(ensg %in% names(ensgExp))) {
			# convert gene id to transcript ids]
			ensg2enst[[ensg]] <- query.biomart.g2t(ensg,ensembl)
			enst = ensg2enst[[ensg]]
			# collect expression across transcripts
			# (transcript x sample expression matrix) -> (1 x sample expression matrix)
			exp.mat.t <- exp.mat[which(rownames(exp.mat) %in% enst),]
			if(!is.null(dim(exp.mat.t))) {
				exp.g <- apply(exp.mat.t,2,collect.transcripts) #apply collect function to every sample/column
				ensgExp[[ensg]] = exp.g
			}else if(is.vector(exp.mat.t)) {
				ensgExp[[ensg]] = exp.mat.t
			}else{
				ensgExp[[ensg]] = NULL
			}
		}
		if(!is.null( ensgExp[[ensg]] ) & !is.nan( mean( ensgExp[[ensg]] ) ) & !is.na( mean( ensgExp[[ensg]] ) ) ) {
			pathGeneExp[[ensg]] <- ensgExp[[ensg]]
		}
	}
	pathGeneExp.mat = do.call(rbind,pathGeneExp)				# assemble single pathway expression matrix: genes x samples
	pathExp[[path]] = apply(pathGeneExp.mat,2,collect.genes)	# collect single pathway expression matrix: pathway x samples
}
pathExp.mat = do.call(rbind,pathExp)
save(pathExp.mat,file='~/Documents//Grad/Labs/Couchesne/ASDanyl/Data/rna.pathways.t_median.g_geomMean.rda')
#return(pathExp)

pathwayModel<-function(){
	library(glmnet)
	library(ROCR)

	load('~/Documents/Grad/Labs/Couchesne/ASDanyl/Data/rna.rda')
	load('~/Documents/Grad/Labs/Couchesne/ASDanyl/Data/annot.rna.rda')
	load('~/Documents//Grad/Labs/Couchesne/ASDanyl/Data/rna.pathways.t_median.g_geomMean.rda')
	male.i = which(annot.rna[,2] =='M') # select males
	x = t(rbind(pathExp.mat,rna))[male.i,] # remove males
#	x = t(rna)[male.i,] # remove males
	x = t(pathExp.mat)[male.i,] # remove males
	#plot(mn<-apply(x,2,mean),sd<-apply(x,2,sd),cex=.1)
#	x = x[,which(apply(x,2,mean)>-10 & apply(x,2,sd)>.2)] # filter
	x = x[,which(apply(x,2,mean)>-2 & apply(x,2,sd)>2)] # filter
	y = unlist(annot.rna[male.i,3]) # remove males

	hc<-hclust( d<-dist( x )) # cluster
	cut<-cutree( hc,3)
	indx = as.numeric(names(table(cut))[which.max(table(cut))])
	x = x[cut==indx,] # remove side clusters
	y = y[cut==indx]

	mo.binomial = cv.glmnet(x=x,y=y,family='binomial',type.measure='auc')
	#y.num = ifelse(y=='ASD',1+runif(1)*.01,runif(1)*.01)
	y.num = ifelse(y=='ASD',1,0)
	mo.gaussian = cv.glmnet(x=x,y=y.num,family='gaussian',alpha=.5)
	mo = mo.gaussian
	p= predict(mo,x)
	pred=prediction(p,round(y.num) )
	perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
	plot(perf, col=rainbow(10))
	performance(pred, measure = "auc") 
	# auc = .65-.9 (filter,male,just max cluster,gaussian,rna+paths,alpha=.5)
	# auc = .63 (filter,male,just max cluster,gaussian,rna,alpha=.5)
	# auc = .65 (filter,male,just max cluster,gaussian,paths,alpha=.5)
	
}

############
# query functions
############
query.biomart.g2t<-function(g,mart){
  return( unlist( getBM(attributes = c("ensembl_transcript_id"),
	filters = "ensembl_gene_id", values = g,
	mart = mart) ))
}

query.kegg.h2g<-function(g){
  query <- keggGet(gene.hsa)
  return( strsplit(query[[1]]$DBLINKS[6],' ')[[1]][2] )
}
