qfbart <- function(Y,X,sl.X,X.out,train.start,
                   Ymu,Ysd,prior_obj){
  list2env(prior_obj,envir=globalenv())
  ntot <- nsave+nburn
  
  T <- NROW(Y)
  K <- ncol(sl.X)
  mcmc_output <- TRUE
  
  # frequentist variant of the model
  if(mod == "AR-abg-fq"){
    pe.quants <- array(NA,dim=c(4,length(set.p),ncol(Y)))
    quant.mean <- array(NA,dim=c(4,length(set.p),T,ncol(Y)))
    dimnames(pe.quants) <- list(c("low","med","high","mean"),set.p,colnames(Y))
    dimnames(quant.mean) <- list(c("low","med","high","mean"),set.p,NULL,colnames(Y))
    for(i in 1:N){
      Y.i <- Y[,i]
      X.i <- X[,sl.X[i,],drop=F]
      X.out.i <- X.out[sl.X[i,]]
      sim <- rq(Y.i~X.i-1,tau=set.p)
      
      for(mm in 1:4){
        pe.quants[mm,,i] <- (apply(as.matrix(sim$coefficients)*X.out.i,2,sum) * Ysd[i]) + Ymu[i]
        quant.mean[mm,,,i] <- (apply(sim$fitted.values,1,sort) * Ysd[i]) + Ymu[i]
      }
    }
    
    ret.list <- list("quants"=quant.mean,            # insample fitted quantiles
                     "predq"=pe.quants               # quantile estimates of quantile forecast
    )
  }
  
  # bayesian implementation and nested models
  if(mod != "AR-abg-fq"){
    # prior on the weights
    M <- length(set.p)
    b <- rep(1, M)
    a <- rep(1, M)
    
    if (as.character(wghs) == "non-parm"){
      est.weights <- FALSE
      omega.start <- 1
    }else if (as.character(wghs) == "parm"){
      est.weights <- FALSE
      omega.start <- 0
    }else if(as.character(wghs) == "mix"){
      est.weights <- FALSE
      omega.start <- 0
    }
    omega <- matrix(omega.start, M, N)
    
    # factor model specification
    if (R>0 & fm){
      f <- matrix(rnorm(R*T, 0, 0.01), T, R)
      Lambda <- matrix(0, N*M, R)
    }else{
      f <- matrix(0, T, R)
      Lambda <- matrix(0, N*M, R)
    }
    
    list.tree.eq <- list()
    list.sampler.run  <- list()
    
    # pooling
    sl.dom <- matrix(FALSE,K,N)
    for(i in 1:N){
      sl.dom[,i] <- grepl(paste(c(substr(colnames(Y)[i],1,3),"cons"),collapse="|"),colnames(X)[sl.X[i,]])
    }
    K.dom <- unique(apply(sl.dom,2,sum))
    
    # initialize BART
    control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                             keepTrees = FALSE, n.samples = 1,
                             n.cuts = 100L, n.burn = 0, n.trees = num.trees, n.chains = 1,
                             n.threads = 1, n.thin = 1L, printEvery = 1,
                             printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                             updateState = FALSE)
    
    for (ii in seq_len(N)){
      sampler.list <- list()
      sampler.run <- list()
      for (jj in seq_len(M)){
        prior.sig = c(10000^50, 0.5)
        sigma.init <- 1
        
        sampler.list[[jj]] <- dbarts(Y[,ii]~X[,sl.X[ii,]], control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=sigma.init, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))
      }
      list.tree.eq[[ii]] <- sampler.list
      list.sampler.run[[ii]] <- sampler.run
    }
    
    # initialize stochvol and set priors
    if (sv){
      sv.draw.idio <- list()
      sv.latent.idio <- list()
      
      sv.draw.fac <- list()
      sv.latent.fac <- list()
      for (mm in seq_len(R)){
        sv.draw.fac[[mm]] <- list(mu = 0, phi = 0.95, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
        sv.latent.fac[[mm]] <- rep(0,T)
      }
      
      sv_priors <- list()
      for(mm in 1:R){
        sv_priors[[mm]] <- specify_priors(
          mu = sv_constant(0), # -4, 1e-3
          phi = sv_beta(shape1 = 5, shape2 = 1.5),
          sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*B_h)),
          nu = sv_infinity(),
          rho = sv_constant(0)
        ) 
      }
    }
    H <- matrix(0, T, R)
    
    beta.mat <- array(NA, c(K, M, N))
    Vprior.array <- array(1, c(K, M, N))
    for (ii in seq_len(N)){
      for (jj in seq_len(M)){
        beta.mat[,jj,ii] <- lm(Y[,ii]~X[,sl.X[ii,]]-1)$coefficients
      }
    }
    beta.prior <- beta.mat*0 # stays zero if not updated in mode pool=TRUE
    
    # compute theta and tau^2  (this is the same across equations)
    theta <- (1-2*set.p)/(set.p * (1 - set.p))
    tau2 <- 2/(set.p * (1-set.p))
    sigma <- matrix(1, M, N)
    
    # Starting values for auxiliary variable in Gaussian representation of ALD
    v <- array(1, c(T, M, N))
    
    # storage
    omega.store <- array(NA, c(nsave, M,N))
    beta.store <- array(NA, c(nsave, K, M, N))
    quant.store <- array(NA, c(nsave, T, M, N))
    f.store <- array(NA, c(nsave, T, R))
    Lambda.store <- array(NA, c(nsave, M, R*N))
    pred.store <- array(NA, c(nsave,M,N))
    count.store <- array(NA, c(nsave, K-1, N, M))
    
    y.quantiles <- array(NA, c(T, M, N))
    y.quantiles_parts <- array(NA, c(T, M, N, 2))
    eta  <- array(NA, c(T, M, N))
    
    #Starting values for the factors and loadings
    sigma2.draw <- rep(1, M*N)
    
    # starting values for the HS prior
    bigpi <- 1
    nu <- 1
    xi <- 1
    zeta <- 1
    
    if(beta.setup=="flat"){
      V.prior.mat <- array(10^10, c(K, M,N))
    }else if(beta.setup=="tight"){
      V.prior.mat <- array(1, c(K, M,N))
    }else{
      V.prior.mat <- array(1, c(K, M,N))
    }
    ind.restr <- matrix(seq(N*M), N, M)
    
    # ------------------------------------------------------------------------------
    # starting sampling algorithm
    pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
    for (irep in seq_len(ntot)){
      # This block samples the parameters of the equations associated with the different quantiles
      var.mat <- array(NA, c(T, M, N))
      count.mat <- array(NA, c(K-1,M,N))
      
      for (i in seq_len(N)){
        X.i <- X[, sl.X[i,],drop=F]
        Lambda.i <- Lambda[seq(i, N*M, by=N),]
        for (q in seq_len(M)){
          # These are quantities necessary to render the regression and BART part conditionally Gaussian and homoscedastic
          v.i <- v[,q,i]
          var.i <- tau2[q]*sigma[q,i]*v.i
          var.mat[,q,i] <- var.i
          
          if(as.character(wghs) == "mix"){
            ytildeq <- (Y[,i] - (X.i%*%beta.mat[, q,i]) - theta[q]*v.i - (f%*%t(Lambda.i))[,q])
            var.qt <- var.i
            
            # We need to re-define the response and the scaling parameter for BART to work
            list.tree.eq[[i]][[q]]$setResponse(ytildeq)
            list.tree.eq[[i]][[q]]$setWeights(1/var.qt)
            
            rep_mm <- list.tree.eq[[i]][[q]]$run(0L, 1L)
            list.sampler.run[[i]][[q]] <- rep_mm
            
            count.mat[,q,i] <- t(rep_mm$varcount)/num.trees
            
            # Now, conditional on the tree sample the regression coefficients 
            yhatq <- (Y[,i] - rep_mm$train - theta[q]*v.i  - (f%*%t(Lambda.i))[,q])/(sqrt(var.i))
            Xhatq <- X.i/(sqrt(var.i))
            
            if (K > 1) V.prior.inv <- diag(1/V.prior.mat[,q,i]) else V.prior.inv <- 1/V.prior.mat[,q,i]
            V.q <- solve(crossprod(Xhatq) + V.prior.inv)
            m.q <- V.q %*% (V.prior.inv %*% beta.prior[,q,i] + crossprod(Xhatq,yhatq))
            m.draw <- m.q + t(chol(V.q))%*%rnorm(K)
            beta.mat[,q,i] <- m.draw
            
            # Draw the parameters used for the Gauss approximation part
            d.q.2 <- as.numeric(sqrt(theta[q]^2 + 2*tau2[q]) / abs(Y[,i] - rep_mm$train - (X.i%*%m.draw) - (f%*%t(Lambda.i))[,q]))
            g.q.2 <- (theta[q]^2 + 2*tau2[q]) / (sigma[q,i] * tau2[q])
            v.i <- 1/rinvgauss(T,mean=d.q.2,dispersion=1/g.q.2)
            v[,q,i] <- v.i
            
            # Draw sigma
            if(ald.scale){
              n.tilda <- (n0 + 3*T)/2
              s.tilda <- (s0 + 2*sum(v.i) + sum((Y[,i] - rep_mm$train - (X.i%*%m.draw) - theta[q]*v.i - (f%*%t(Lambda.i))[,q])^2/(tau2[q] * v.i)))/2
              sigma[q,i] <- 1/rgamma(1, n.tilda, s.tilda)
            }else{
              sigma[q,i] <- 1
            }
            
            y.quantiles[ ,q,i] <- rep_mm$train + (X.i%*%m.draw) + (f%*%t(Lambda.i))[,q]
            y.quantiles_parts[,q,i,1] <- rep_mm$train
            y.quantiles_parts[,q,i,2] <- (X.i%*%m.draw)
            
            if (R > 0){
              # Construct shocks for estimating the latent factor
              eta[, q,i] <- Y[,i] - rep_mm$train - (X.i%*%m.draw)  - v.i*theta[q]
            }
          }else{
            if (omega[q,i]==1){
              ytildeq <- (Y[,i] - (X.i%*%beta.mat[, q,i]) * (1-omega[q,i]) - theta[q]*v.i - (f%*%t(Lambda.i))[,q])/omega[q,i]
              var.qt <- var.i/omega[q,i]^2
              
              # We need to re-define the response and the scaling parameter for BART to work
              list.tree.eq[[i]][[q]]$setResponse(ytildeq)
              list.tree.eq[[i]][[q]]$setWeights(1/var.qt)
              
              rep_mm <- list.tree.eq[[i]][[q]]$run(0L, 1L)
              list.sampler.run[[i]][[q]] <- rep_mm
              
              count.mat[,q,i] <- t(rep_mm$varcount)/num.trees
              m.draw <- rep(0,K)
            }else{
              rep_mm <- list()
              rep_mm$train <- rep(0, T)
            }
            
            if(omega[q,i]==0){
              # Now, conditional on the tree sample the regression coefficients 
              yhatq <- (Y[,i] - rep_mm$train * omega[q,i] - theta[q]*v.i  - (f%*%t(Lambda.i))[,q])/(sqrt(var.i))
              Xhatq <- X.i * (1-omega[q,i])/(sqrt(var.i))
              
              if (K > 1) V.prior.inv <- diag(1/V.prior.mat[,q,i]) else V.prior.inv <- 1/V.prior.mat[,q,i]
              V.q <- solve(crossprod(Xhatq) + V.prior.inv)
              m.q <- V.q %*% (V.prior.inv %*% beta.prior[,q,i] + crossprod(Xhatq,yhatq))
              m.draw <- m.q + t(chol(V.q))%*%rnorm(K)
              beta.mat[,q,i] <- m.draw
            }
            
            # Draw the parameters used for the Gauss approximation part
            d.q.2 <- as.numeric(sqrt(theta[q]^2 + 2*tau2[q]) / abs(Y[,i] - rep_mm$train*omega[q,i] - (X.i%*%m.draw)*(1-omega[q,i]) - (f%*%t(Lambda.i))[,q]))
            g.q.2 <- (theta[q]^2 + 2*tau2[q]) / (sigma[q,i] * tau2[q])
            v.i <- 1/rinvgauss(T,mean=d.q.2,dispersion=1/g.q.2)
            v[,q,i] <- v.i
            
            # Draw sigma
            if(ald.scale){
              n.tilda <- (n0 + 3*T)/2
              s.tilda <- (s0 + 2*sum(v.i) + sum((Y[,i] - rep_mm$train*omega[q,i] - (X.i%*%m.draw)*(1-omega[q,i]) - theta[q]*v.i - (f%*%t(Lambda.i))[,q])^2/(tau2[q] * v.i)))/2
              sigma[q,i] <- 1/rgamma(1, n.tilda, s.tilda)
            }else{
              sigma[q,i] <- 1
            }
            
            y.quantiles[ ,q,i] <- rep_mm$train*omega[q,i] + (X.i%*%m.draw)*(1-omega[q,i]) + (f%*%t(Lambda.i))[,q]
            y.quantiles_parts[,q,i,1] <- rep_mm$train
            y.quantiles_parts[,q,i,2] <- (X.i%*%m.draw)
            
            if (R > 0){
              # Construct shocks for estimating the latent factor
              eta[, q,i] <- Y[,i] - rep_mm$train*omega[q,i] - (X.i%*%m.draw)*(1-omega[q,i])  - v.i*theta[q] 
            }
          }
        }
      }
      
      # Step II: Sample shrinkage parameters for the HS
      if(!(beta.setup %in% c("tight","flat"))){
        if(pool){
          for(q in 1:M){
            V.q <- matrix(V.prior.mat[,q,][sl.dom],K.dom,N)
            beta.q <- matrix(beta.mat[,q,][sl.dom],K.dom,N)
            
            V.q.inv <- diag(1/(apply(1/V.q,1,sum) + 1/beta.pool.tau))
            beta.q.mu <- V.q.inv %*% apply(beta.q/V.q,1,sum)
            beta.pool <- beta.q.mu + t(chol(V.q.inv)) %*% rnorm(K.dom,0,1)
            beta.prior[,q,][sl.dom] <- beta.pool
          }
        }
        
        hs_draw <- get.hs(bdraw=(as.vector(beta.mat)-as.vector(beta.prior)),lambda.hs=bigpi,nu.hs=nu,tau.hs=xi,zeta.hs=zeta)
        bigpi <- hs_draw$lambda
        nu <- hs_draw$nu
        xi <- hs_draw$tau
        zeta <- hs_draw$zeta
        V.prior.mat <- array(hs_draw$psi, c(K, M, N))
        V.prior.mat[V.prior.mat<1e-8] <- 1e-8
        V.prior.mat[V.prior.mat>10] <- 10
      }
      
      # Step III: Sample the common factor in included in X
      # Before we do this, create matrices eta.mat and Var.mat, we need to do this such that eta.mat = (eta_1p1, eta_2p1, .., eta_Np1, eta_1p2, eta_2p2, ..)
      eta.mat <- NULL
      Var.mat <- NULL
      for (jj in seq_len(M)){
        eta.mat <- cbind(eta.mat, eta[,jj,])
        Var.mat <- cbind(Var.mat, var.mat[,jj,])
      }
      
      if (R>0 & fm){
        for (t in seq_len(T)){
          normalizer <- 1/sqrt(Var.mat[t, ])
          Lt <- Lambda * normalizer
          yt <- eta.mat[t, ] * normalizer
          
          if (R > 1) Q <- diag(1/exp(H[t,])) else Q <- 1/exp(H[t])
          fac.Sigma <-  solve(crossprod(Lt)+Q)
          fac.mean <- fac.Sigma%*%crossprod(Lt,yt)
          
          if (R > 1 ){
            f[t,] <- fac.mean+t(chol(fac.Sigma))%*%rnorm(R) 
          } else {
            f[t] <-  fac.mean+t(chol(fac.Sigma))%*%rnorm(R)
          } 
        }
        f <- apply(f,2,function(x){(x-mean(x))/sd(x)}) # normalize factor draw
        
        if (sv){
          # Step IV: Sample the factor volatilities
          sv.para <- matrix(NA, R, 3) # first dim = mu, second = rho, third = sigma
          for (jj in 1:R){
            svdraw_mm <- svsample_general_cpp(f[,jj], startpara = sv.draw.fac[[jj]], startlatent = sv.latent.fac[[jj]], priorspec = sv_priors[[jj]])
            sv.draw.fac[[jj]][c("mu", "phi", "sigma")] <- as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
            sv.latent.fac[[jj]] <- svdraw_mm$latent
            H[,jj] <- svdraw_mm$latent
            sv.para[jj,] <- svdraw_mm$para[, c("mu", "phi", "sigma")]
          }
        }else{
          H[] <- 0
        }
        
        for (j in seq_len(N*M)){
          sl.fac <- which(ind.restr==j,arr.ind = TRUE)[[2]]
          normalizer <- 1/sqrt(Var.mat[, j])
          yj <- eta.mat[,j]*normalizer
          fj <- f[,sl.fac]*normalizer
          prior.v <- fac.var
          
          V.f <- solve(crossprod(fj)+prior.v)
          m.lambda.f <- V.f %*% crossprod(fj, yj)
          lambda.draw <- m.lambda.f + t(chol(V.f))%*%rnorm(1)
          Lambda[j,sl.fac] <- lambda.draw
        }
        
        # identify sign of the factor
        for(rr in 1:R){
          L.tmp <- Lambda[,rr]
          L.sign <- sign(mean(L.tmp[L.tmp!=0]))
          Lambda[,rr] <- Lambda[,rr]*L.sign
          f[,rr] <- f[,rr]*L.sign
        }
      }
      
      pred.mat <- array(0, c(M,N)); dimnames(pred.mat) <- list(set.p,colnames(Y))
      if (irep > nburn){
        pred.mat <- array(NA, c(M,N)); dimnames(pred.mat) <- list(set.p,colnames(Y))
        
        # Shock DFM and predict log-volas
        if (fm){
          if (sv){
            HT1 <- sv.para[,1] + sv.para[,2]*(H[T,] - sv.para[,1]) + sv.para[,3]*rnorm(R)
          }else{
            HT1 <- H[T,]
          }
          HT1[exp(HT1/2)>20] <- 2*log(20) # offsetting for pandemic
          f.shock <- Lambda%*%rnorm(R, 0, exp(HT1/2))
        }else{
          f.shock <- Lambda%*%rep(0,R)
        }
        
        # use samples from holdout
        X.p <- X.out
        for (i in seq_len(N)){
          f.i <- f.shock[seq(i, N*M, by=N)]
          
          for (j in seq_len(M)){
            pred.tree <- list.tree.eq[[i]][[j]]$predict(X.p[sl.X[i,-ncol(sl.X)]])
            pred.reg <- X.p[sl.X[i,]]%*%beta.mat[,j,i]
            
            if(as.character(wghs) == "mix"){
              pred.t <- pred.tree + pred.reg + f.i[j]
            }else{
              pred.t <- pred.tree*omega[j,i] + pred.reg*(1-omega[j,i]) + f.i[j]
            }
            
            pred.mat[j,i] <- pred.t
          }
        }
        
        # rescaling
        for(i in 1:N){
          pred.mat[,i] <- (pred.mat[,i] * Ysd[i]) + Ymu[i]
          y.quantiles[,,i] <- (y.quantiles[,,i] * Ysd[i]) + Ymu[i]
        }
        
        # storage
        pred.store[irep-nburn,,] <- pred.mat
        
        if(as.character(wghs) == "mix"){
          var_quantiles <- apply(y.quantiles_parts,c(2,3,4),sd)^2
          omega.store[irep-nburn,,] <- var_quantiles[,,1] / apply(var_quantiles,c(1,2),sum)
        }else{
          omega.store[irep-nburn,,] <- omega
        }
        
        beta.store[irep-nburn,,,] <- beta.mat
        quant.store[irep-nburn,,,] <- y.quantiles
        f.store[irep-nburn,,] <- f
        Lambda.store[irep-nburn,,] <- Lambda
        count.store[irep-nburn,,,] <- count.mat
      }
      
      # progress bar
      if(!silent) setTxtProgressBar(pb, irep)
    }
    
    pe.quants <- apply(pred.store, c(2,3), quantmean.fun)
    dimnames(pe.quants) <- list(c("low", "med", "high", "mean"),
                                  set.p,
                                  colnames(Y))
    
    quant.mean <- apply(quant.store, c(2,3,4), quantmean.fun) # insample estimate of the quantiles
    dimnames(quant.mean) <- list(c("low","med","high","mean"),
                                  as.character(rownames(Y)),
                                  set.p,
                                  colnames(Y))
    
    weights.quant <- apply(omega.store,c(2,3),quantmean.fun)
    beta.quant <- apply(beta.store, c(2,3,4), quantmean.fun)
    f.quant <- apply(f.store, c(2,3), quantmean.fun)
    L.quant <- apply(Lambda.store, c(2,3), quantmean.fun)
    count.quant <- apply(count.store, c(2,3,4), quantmean.fun)
    
    dimnames(weights.quant) <- list(c("low","med","high","mean"),set.p,colnames(Y))
    dimnames(beta.quant) <- list(c("low","med","high","mean"),
                                 NULL,
                                 set.p,
                                 colnames(Y))
    dimnames(f.quant) <- list(c("low","med","high","mean"),
                              as.character(rownames(Y)),
                              set.p)
    
    ret.list <- list("weights"=weights.quant,        # weight on the non-parametric part
                     "beta"=beta.quant,              # "linear" part of the coefficients
                     "quants"=quant.mean,            # insample fitted quantiles
                     "predq"=pe.quants,              # quantile estimates of quantile forecast
                     "factors"=f.quant,              # latent factors wrt. cross-section
                     "loadings"=L.quant,             # loadings on the factors
                     "count"=count.quant            # splitting rule count
    )
  }
  return(ret.list)
}

