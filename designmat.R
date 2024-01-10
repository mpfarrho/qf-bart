design.matrix <- function(data_raw,sl.cN,P,fhorz,
         mod,sl.norm,
         length.hold.out, h){
  require(zoo)
  
  rownames(data_raw) <- as.character(zoo::as.Date.ts(data_raw))
  train.start <- start(data_raw)+c(0,P) # timing for training set wrt. lag structure
  
  # extract data
  Yraw <- data_raw[1:(nrow(data_raw)-length.hold.out+h-1),grepl(paste0(sl.cN,collapse="|"),colnames(data_raw))]
  if(sl.norm){
    Ymu <- apply(Yraw,2,mean)[paste0(sl.cN,"_GDP")]
    Ysd <- apply(Yraw,2,sd)[paste0(sl.cN,"_GDP")]
    Yraw <- apply(Yraw,2,function(x) (x-mean(x))/sd(x))
  }else{
    Ymu <- Ysd <- rep(1,length(sl.cN))
  }
  
  Ylag <- embed(Yraw,fhorz+1)
  colnames(Ylag) <- paste0(rep(colnames(Yraw),fhorz+1),"_lag",rep(0:fhorz,each=ncol(Yraw)))
  rownames(Ylag) <- rownames(Yraw)[-c(1:fhorz)]
  
  N <- length(sl.cN)
  Y <- Ylag[,paste0(sl.cN,"_GDP_lag",0)]#apply(Ylag[,paste0(sl.cN,"_GDP_lag",0)],2,avg_h,n=fhorz)
  X <- cbind(Ylag[,c(paste0(c(sl.cN),"_GDP_lag",fhorz),
                     paste0(c(sl.cN),"_FSI_lag",fhorz))],1)
  
  # exclude initial observations
  Y <- Y[fhorz:nrow(Y),]
  X <- X[fhorz:nrow(X),]
  
  colnames(X)[ncol(X)] <- "cons"
  X.out <- c(Ylag[nrow(Ylag),c(paste0(c(sl.cN),"_GDP_lag",0),
                               paste0(c(sl.cN),"_FSI_lag",0))],1)
  
  # construct selection indexes
  if(mod == "AR-abg" || mod == "AR-abg-fq"){
    sl.X <- NULL
    for(nn in 1:N) sl.X <- rbind(sl.X,which(grepl(sl.cN[nn],colnames(X)))) # only select domestic variables
    sl.X <- cbind(sl.X,rep(ncol(X),N))
  }else if(mod %in% c("AR-abg-m")){
    sl.X <- matrix(rep(1:ncol(X),each=N),N,ncol(X)) # select all variables and lags
    sl.X <- sl.X[,!grepl(paste0(c("PC_GDP","PC_FSI"),collapse="|"),colnames(X))] # exclude principal components
  }
  T <- nrow(Y)
  
  design <- list("Y"=Y,"X"=X,
                 "X.out"=X.out,"train.start"=train.start,
                 "sl.X"=sl.X,"Ymu"=Ymu,"Ysd"=Ysd,
                 "T"=T,"N"=N)
  return(design)
}
