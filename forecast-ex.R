run <- 1

library(dbarts)
library(GIGrvg)
library(stochvol)
library(quantreg)
library(statmod)

source("utils.R")
source("designmat.R")
source("qf-bart.R")
load("data_raw.rda")

silent <- FALSE
sl.norm <- TRUE
length.hold.out <- 130
sl.cN <- c("DE","FR","IT","UK","US","AT","DK","ES","FI","NL","SE")

fm.grid <- c(FALSE,TRUE)
horz.grid <- 1:4
weights.sample <- c("parm","non-parm","mix")

h.out <- seq(1, length.hold.out)
mod.grid <- c("AR-abg","AR-abg-m")
pool.grid <- c(FALSE,TRUE)

# set up grids
grid0 <- expand.grid("h"=h.out,
                     "fhorz"=horz.grid,
                     "fm"=FALSE,
                     "wghs"="parm",
                     "mod"="AR-abg-fq",
                     "pool"=FALSE,
                     "ald"=c("fixed"),
                     "bart"=c("v1"),
                     stringsAsFactors = FALSE)
grid1 <- expand.grid("h"=h.out,
                     "fhorz"=horz.grid,
                     "fm"=fm.grid,
                     "wghs"=weights.sample,
                     "mod"=mod.grid,
                     "pool"=pool.grid,
                     "ald"="regular",
                     "bart"=c("v1"),
                     stringsAsFactors = FALSE)

grid <- rbind(grid0,grid1)
rownames(grid) <- 1:nrow(grid)

# select corresponding specification
grid.slct <- grid[run,]
h <- as.numeric(grid.slct[["h"]])
fhorz <- as.numeric(grid.slct[["fhorz"]])
fm <- as.logical(grid.slct[["fm"]])
wghs <- grid.slct[["wghs"]]
mod <- grid.slct[["mod"]]  
pool <- grid.slct[["pool"]]

# prior setup for BART and ALD
ald <- grid.slct[["ald"]]
bart.v <- grid.slct[["bart"]]

# storage of results
foldername <- paste0("results")
dir.create(foldername, showWarnings = FALSE)
filen <- paste0(formatC(h,flag="0",width=3),"_h",fhorz,
                "_fm",fm,"_w",wghs,"_mod",mod,"_pool",pool,
                "_ald",ald,"_bart",bart.v)
filename1 <- paste0(foldername,"/",filen,".rda")

# ---------------------------------------------------------------------------------------------------
# model selection
P <- 1 # lags
sv <- TRUE # stochastic volatility
set.p <- c(0.05,0.1,0.16,0.25,0.4,0.5,0.6,0.75,0.84,0.9,0.95) # quantile function argument # seq(0.05,0.95,by=0.05)
R <- length(set.p) # Number of latent factors = one factor per quantile (cross-section information)

# prior settings
fac.var <- 1
beta.pool.tau <- 10
B_h <- 1 # prior on state innovations in H 

# settings for ALD scale
if(ald=="regular"){
  ald.scale <- TRUE
  n0 <- rep(1,R)
  s0 <- rep(1,R)
}else if(ald=="fixed"){
  ald.scale <- TRUE
  n0 <- rep(0,R)
  s0 <- rep(0,R)
}

# set hyperparameters for BART
cgm.level <- 0.95 # default value
cgm.exp <- 2 # default value

if(bart.v=="v1"){
  sd.mu <- 2 # default
  num.trees <- 250 # default
  
  beta.setup <- "hs"
}

# ---------------------------------------------------------------------------------------------------
# compute design matrix
design_ls <- design.matrix(data_raw=data_raw,sl.cN=sl.cN,P=P,fhorz=fhorz,mod=mod,
                           sl.norm=sl.norm,
                           length.hold.out=length.hold.out,h=h)
list2env(design_ls,envir=globalenv())

# construct prior settings object in list
prior_obj <- list("nsave"=5000,"nburn"=5000,"silent"=silent,
                  "mod"=mod,"fm"=fm,"wghs"=wghs,"pool"=pool,
                  "set.p"=set.p,"R"=length(set.p),"sv"=sv,
                  fac.var=fac.var,"n0"=n0,"s0"=s0,"ald.scale"=ald.scale,"B_h"=B_h,
                  "cgm.level"=cgm.level,"cgm.exp"=cgm.exp,"sd.mu"=sd.mu,"num.trees"=num.trees,
                  "beta.setup"=beta.setup,"beta.pool.tau"=beta.pool.tau
)

if(!file.exists(filename1)){
  ret.list <- qfbart(Y=Y,X=X,sl.X=sl.X,X.out=X.out,train.start=train.start,
                     Ymu=Ymu,Ysd=Ysd,prior_obj=prior_obj)
  save(file=filename1, ret.list)
}


