load('~/Documents/USEFUL.CODE/hg19/hg19.genes.gtf.rda')
load('~/Documents/USEFUL.CODE/hg19/hg19.transcript.gtf.rda')

gene2transcript<-function(g,transcript){
	tsxL <- lapply(g,function(x) transcript$transcript_id[transcript$gene_id==x])
	names(tsxL) = g
	return(tsxL)
}

transcript2gene<-function(t,transcript){
	gnL <- lapply(t,function(x) transcript$gene_id[transcript$transcript_id==x])
	names(gnL) = t
	return(gnL)
}

gene2name<-function(g,genes){
	gnL <- lapply(g,function(x) genes$gene_name[genes$gene_id==x])
	names(gnL) = g
	return(gnL)
}

transcript2name<-function(t,transcript){
	tsxL <- lapply(t,function(x) transcript$transcript_name[transcript$transcript_id==x])
	names(tsxL) = t
	return(tsxL)
}

transcript2genename<-function(t,transcript){
	tsxL <- lapply(t,function(x) transcript$gene_name[transcript$transcript_id==x])
	names(tsxL) = t
	return(tsxL)
}

transcript_id2bed<-function(t,transcript){
	tsxL <- lapply(t,
		function(x) cbind(	transcript$start[transcript$transcript_id==x],
							transcript$end[transcript$transcript_id==x] )  )
	names(tsxL) = t
	return(tsxL)
}

gene_id2bed<-function(g,genes){
	gnL <- lapply(g,
		function(x) cbind(	genes$start[genes$gene_id==x],
							genes$end[genes$gene_id==x] ) )
	names(gnL) = g
	return(gnL)
}

transcript2gene.exprss<-function(t,m,transcript,collect.genes=max){
	g <- unique(unlist( transcript2gene(t,transcript) ))
	t.g <- gene2transcript(g,transcript)
	g.exprss <- list()
	for(gene in t.g){
		g.exprss[[names(gene)]] <- apply(m[which(rownames(m) %in% gene),] , 2 , collect.genes)
	}
}


pre.proc<-function(){
	if(f){
#	system('grep \tgene\t Homo_sapiens.GRCh38.78.gtf > hg19.genes.gtf
#			grep \ttranscript\t Homo_sapiens.GRCh38.78.gtf > hg19.transcript.gtf
#
#			sed s/gene_id\s// \
#				s/gene_name\s// \
#				s/gene_source\s// \
#				 hg19.genes.gtf > hg19.genes.sed.gtf')
#	system("grep -v '\tensembl_havana\t' hg19/hg19.transcript.gtf > tmp.txt")
	}

	gtf = read.csv('hg19/Homo_sapiens.GRCh38.78.rmheader.gtf',sep='\t',header=F)
	save(gtf,file='hg19/Homo_sapiens.GRCh38.78.gtf.rda')

	genes = read.csv('hg19/hg19.genes.gtf',sep='\t')
	transcript = read.csv('hg19/hg19.transcript.gtf',sep='\t')

	save(genes,file='hg19/hg19.genes.gtf.rda')
	save(transcript,file='hg19/hg19.transcript.gtf.rda')
}
