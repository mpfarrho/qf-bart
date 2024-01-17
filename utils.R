# remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 2 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

# useful for plotting boxplot with mean
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data = fun, geom = geom, width = 0.7, ...)
}

stat_sum_single <- function(fun, geom="crossbar", ...) {
  stat_summary(fun=fun, geom=geom, width = 0.7, linewidth = 0.3, ...)
}

# define parameters of the beta distribution by mean and variance
beta_ab <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}

# function to save median and mean alongside quantiles
quantmean.fun <- function(x){
  return(c(quantile(x, c(0.16,0.5,0.84),na.rm=TRUE),mean(x,na.rm=TRUE)))
}

# diebold mariano test wrapper
dm.test_dplyr <- function(x,y,alternative="less",h,power,ret="stars"){
  require(forecast)
  e1 <- ts(x[(is.na(x) + is.na(y)) == 0])
  e2 <- ts(y[(is.na(x) + is.na(y)) == 0])
  
  if(sum(is.na(e1))==length(e1) || sum(is.na(e2))==length(e2)){
    return(NA)
  }else if(sum(e1 - e2) == 0){
    return(NA)
  }else{
    pval <- dm.test(e1=e1,e2=e2,alternative="two.sided",h=h,power=power)$p.value
    if(ret=="stars"){
      if(pval>0.1){
        return("")
      }else if(pval > 0.05 & pval <= 0.10){
        return("'")
      }else if(pval > 0.01 & pval <= 0.05){
        return("°")
      }else if(pval <= 0.01){
        return("*")
      }
    }else{
      return(pval)
    }
  }
}

# compute random numbers from implied distribution based on quantiles
rnum_from_p <- function(yp,n=3000){
  dens <- density(yp)
  return(sample(x = dens$x, n, prob = dens$y, replace=TRUE) + rnorm(n, 0, dens$bw))
}

dens_from_p <- function(y=NULL,fcst_p,mode="p",n=3000){
  dens <- density(fcst_p,n=n)
  dens.xy.norm <- cbind(dens$x,dens$y/sum(dens$y))
  dens.xy <- cbind("x"=dens$x,"y"=dens$y)
  dens.p <- cbind("x"=dens$x,"y"=cumsum(dens.xy.norm[,2]))
  
  if(is.null(y)){
    return(dens.xy)
  }else{
    if(mode=="d"){
      lps <- as.numeric(log(approx(dens.xy[,"x"],dens.xy[,"y"],y)$y))
      return(lps)
    }else if(mode=="p"){
      pit <- approx(dens.p[,"x"],dens.p[,"y"],y)
      return(pit)
    }
  }
}

# compute average growth rate for horizon
avg_h <- function(x, n){stats::filter(x, rep(1 / n, n), method="convolution", sides = 1)}

# function to draw the factor loadings (basic linear regression)
get.facload <- function(yy,xx,l_sd){
  V_prinv <- diag(NCOL(xx))/l_sd
  V_lambda <- solve(crossprod(xx) + V_prinv)
  lambda_mean <- V_lambda %*% (crossprod(xx,yy))
  
  lambda_draw <- lambda_mean + t(chol(V_lambda)) %*% rnorm(NCOL(xx))
  return(lambda_draw)
}

# factor loadings draw
get.Lambda <- function(eps,fac,S,pr,m,q,id.fac){
  L <- matrix(0,m,q)
  if(id.fac){
    for(jj in 1:m){
      if (jj<=q){
        normalizer <- sqrt(S[,jj])
        yy0 <- (eps[,jj]-fac[,jj])/normalizer
        xx0 <- fac[,1:(jj-1),drop=FALSE]/normalizer
        if (jj>1){
          l_sd <- pr[jj,1:(jj-1)]
          lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
        }else{
          lambda0 <- 1
        }
        
        if (jj>1){
          L[jj,1:(jj-1)] <- lambda0
          L[jj,jj] <- 1
        }else if (jj==1){
          L[jj,jj] <- 1
        }
      }else{
        normalizer <- sqrt(S[,jj])
        yy0 <- (eps[,jj])/normalizer
        xx0 <- fac[,,drop=FALSE]/normalizer
        l_sd <- pr[jj,]
        lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
        L[jj,] <- lambda0
      }
    }
  }else{
    for(jj in 1:m){
      normalizer <- sqrt(S[,jj])
      yy0 <- (eps[,jj])/normalizer
      xx0 <- fac[,,drop=FALSE]/normalizer
      l_sd <- pr[jj,]
      lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
      L[jj,] <- lambda0
    }
  }
  return(L)
}

# forecast eval functions
QS <- function(true,Qtau,tau,mode="crps"){
  tau_len <- length(tau)
  true_rep <- rep(true,tau_len)
  if(mode=="crps"){
    QS.vec <- 2*(((true_rep<=Qtau)*1)-tau)*(Qtau-true_rep)
  }else{
    QS.vec <- (((true_rep<=Qtau)*1)-tau)*(Qtau-true_rep)
  }
  return(QS.vec)
}

qwCRPS_point <- function(true,Qtau,tau,weighting){
  require(pracma)
  
  tau_len <- length(tau)
  true_rep <- rep(true,tau_len)
  QS.vec <- (true_rep-Qtau)*(tau-((true_rep<=Qtau)*1))
  
  weights <- switch(tolower(weighting),
                    "none" = 1,
                    "tails" = (2*tau-1)^2,
                    "right" = tau^2,
                    "left" = (1-tau)^2,
                    "center" = tau*(1-tau))
  wghs <- QS.vec*weights
  return(pracma::trapz(tau,wghs))
}

# Horseshoe posterior
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}
