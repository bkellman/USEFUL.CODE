multihist<-function(l,main,xlab,legend=NULL,x.leg=NULL,y.leg=NULL){
  if(length(l)<=1){stop('There is only 1 distribution')}
  if(length(l)>3){warning('This function is not optimized for more than 3 distributions')}
  if(class(l)!='list'){stop('l must be a list')}
  if(length(l)!=length(legend) & !is.null(legend)){warning('List and legend length should match')}

  color = rainbow(length(l),alpha=1/4)

  pi<-list()
  param<-list()
  for(i in 1:length(l)){
	pi[[i]] <- hist(l[[i]],plot=F)
	bks = pi[[i]]$breaks
	dns = pi[[i]]$counts
	bkN = length(bks)
	param[[i]] <- c(xmin=bks[1],xmax=bks[bkN],
			ymin=0     ,ymax=max(dns))
  }
  param <- do.call(rbind,param)

  xlim = c(min(param[,1]) - min(param[,1])*.1 , max(param[,2]) + max(param[,2])*.1 )
  ylim = c(min(param[,3]) - min(param[,3])*.1 , max(param[,4]) + max(param[,4])*.1 )
  yorder = order(-param[,4])

  add=F
#  for(i in yorder){
  for(i in 1:length(l)){
	plot( pi[[i]], col=color[i], add=add,main=main,xlab=xlab,xlim=xlim,ylim=ylim)  # second
	add=T
  }
  if(is.null(x.leg)){ x.leg=max(param[,2])*.75 }
  if(is.null(y.leg)){ y.leg=max(param[,4])*.85 }
  if(is.null(legend)){ legend=names(l) }
  legend(x=x.leg,y=y.leg,legend=legend,fill=color)
}


################


multidense<-function(l,main,xlab,legend=NULL,x.leg=NULL,y.leg=NULL,kernel='gaussian'){
  if(length(l)<=1){stop('There is only 1 distribution')}
#  if(length(l)>3){warning('This function is not optimized for more than 3 distributions')}
  if(class(l)!='list'){stop('l must be a list')}
  if(is.null(legend)){ legend=names(l) }
  if(length(l)!=length(legend) & !is.null(legend)){warning('List and legend length should match')}

  color = rainbow(length(l))

  pi<-list()
  param<-list()
  for(i in 1:length(l)){
	pi[[i]] <- density(l[[i]],kernel=kernel)
	x = pi[[i]]$x
	y = pi[[i]]$y
	param[[i]] <- c(xmin=min(x),xmax=max(x),
	   ymin=min(y),ymax=max(y))
  }
  param <- do.call(rbind,param)
  #print(param)
  xlim = c(min(param[,1]) - min(param[,1])*.1 , max(param[,2]) + max(param[,2])*.1 )
  ylim = c(min(param[,3]) - min(param[,3])*.1 , max(param[,4]) + max(param[,4])*.1 )
  yorder = order(-param[,4])

  add=F
#  for(i in yorder){
  for(i in 1:length(l)){
	if(!add){
	  plot( pi[[i]], col=color[i],main=main,xlab=xlab,xlim=xlim,ylim=ylim)  # second
	}else{
	  lines( pi[[i]], col=color[i]) #,ylim=ylim)  # second
	}
	add=T
  }
  if(is.null(x.leg)){ x.leg=max(param[,2])*.75 }
  if(is.null(y.leg)){ y.leg=max(param[,4])*.85 }
  legend(x=x.leg,y=y.leg,legend=legend,fill=color)
}






