library(ROCR)
#library(survcomp)
### prediction evaluation
pred.eval<-function(predi,true,id){
        pdf(paste("ROC.",id,".pdf",sep=""))
        pred<-prediction(predi, true)
        #ci<-round(concordance.index(predi,true,array(TRUE,length(predi)))$c.index,digits=3)
        perf <- performance(pred,"tpr","fpr")
        # changing params for the ROC plot - width, etc
        par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
        # plotting the ROC curve
        plot(perf,col="black",lty=3, lwd=3)
        # calculating AUC & MCC
        auc <- performance(pred,"auc")
        mcc <- performance(pred,"mat")
        # now converting S4 class to vector
        auc <- unlist(slot(auc, "y.values"))
        mcc <- unlist(slot(mcc, "y.values"))[2]
        # adding min and max ROC AUC to the center of the plot
        auc<-round(auc, digits = 3)
        mcc<-round(mcc, digits = 3)

        auct <- paste(c("AUC  = "),auc,sep="")
        mcct <- paste(c("MCC  = "),mcc,sep="")
        #cict <- paste(c("CI  = "),ci,sep="")
        legend(0.3,0.6,c(auct,mcct,"\n"),border="white",cex=1.7,box.col = "white")
        dev.off()
        return(c(auc,mcc))
}
