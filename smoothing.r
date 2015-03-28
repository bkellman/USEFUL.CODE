########
## boxcar smooting
#######
# walks along a rows axis
# applies specificed function (default = mean) to values within specified window along the column
# vars is a vector corresponding to the labels that define the borders of each sub-matrix boxcar smoothing

# vars = [ a a a b b b ]
# a 1 1 2			# a 1 1 2
# a 1 3 2	  a		# a 1 3 2
# a 1 1 2	-->		# a 1 1 2
# b 2 2 1			
# b 1 1 1		b	# b 2 2 1
# b 3 2 1	  -->	# b 1 1 1
					# b 3 2 1

# x = data matrix, variables
boxcar <- function(x,window=2,func=general.mean,vars=NULL,annotation=NULL){
	if(is.null(vars)){
		x_out = boxcar_private(x,window,func)
	}else{
		if(length(vars)!=nrow(x)){stop("vars must correspond to the columns of x")}
		if(sum(table(vars)< window+1)){stop("window size must be greater than the minimum number of repeats in the vars vector")}
		uvars <- unique(vars)
		xL<-list()
		for(v in uvars){
			xL[[v]] <- boxcar_private(x[which(vars==v),],window,func)
		}
		x_out = do.call(rbind,xL)
	}
	outL = list()
	for(i in 1:ncol(x_out)){
		outL[[colnames(x_out)[i]]] = type.convert(x_out[,i])
	}
	x_out = data.frame(outL)
	return(x_out)
}

general.mean <- function(x){
#	if(is.factor(x) | is.character(x)){
	if(sum( is.na(as.numeric(x)) )>0){
		x = as.character(x)
		xu = unique(x)
		return(paste(xu,collapse='-'))
	}else{
		return(mean(as.numeric(x)))
	}
}

boxcar_private <- function(x,window,func){
	xL<-list()
	save(x,file='tmp.rda')
	if(window%%2==0){
		half  = window/2
		min_i  = half+1
		max_i  = nrow(x) - (half  - 1)
		#max_i  = nrow(x) - half  - 1
		#max_i  = nrow(x) - half  
		car_min = -half
		car_max = half-1
 	}else{
 		half  = floor(window/2)
		min_i  = half+1
		max_i  = nrow(x) - half
		car_min = -half
		car_max = half
 	}

	for(i in min_i:max_i){
		x_i = x[(i+car_min):(i+car_max),]
#		print(x_i)
  		xL[[i]] = apply(x_i,2,func)
  	}
  	x_out = do.call(rbind,xL)

  	return(x_out)
}
