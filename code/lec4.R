library(tidyverse)
library(ggmap)
library(classInt)
library(MBA)
library(fields)
library(geoR)
library(RColorBrewer)
library(mgcv)
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
library(bayesplot)
library(loo)
library(MCMCpack)

#### running Bayesian estimation for a subset of dataset 3 in Stan ####
library(rstan)
#options(mc.cores = parallel::detectCores())
data3=read.csv("../data/dataset3.csv")
data3=data3[1:200,]

## recalling the MLE
mle <- likfit(coords=data3[,1:2], data=data3[,4], trend = trend.spatial(~x,data3),
  ini.cov.pars=c(0.12,0.2),nugget = 0.02,cov.model="exponential",nospatial=TRUE)

mle

## stan data file
data3list=list(n=nrow(data3),y=data3$y,
  x=data3$x,dmat=as.matrix(dist(data3[,c("sx","sy")])),id=diag(nrow(data3)))

## running stan
set.seed(1)
fit1 <- stan(
  file = "bayesgp1.stan",  # Stan program
  data = data3list,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  #open_progress = T,
  cores = 4             # number of cores (could use one per chain)
  #refresh = 10             # no shown every 100 iterations
  )

## save(file="../data/data3subset200_bayesgp1.Rdata",list=c("fit1"))
## load("../data/data3subset200_bayesgp1.Rdata")

## run time
get_elapsed_time(fit1)

### posterior estimates, traceplots, densities
params=c("alpha","beta","sigmasq","tausq","phi")
plot(fit1, plotfun = "trace",pars=params, inc_warmup = TRUE)
mcmc_dens_overlay(fit1, pars = params)
#plot(fit1, plotfun = "dens",pars=params, inc_warmup = TRUE)
plot(fit1,pars=params,ci_level=0.95)


### Rhat
round(summary(fit1)[[1]][params,c("Rhat","n_eff","se_mean")],4)

### comparing with the MLE
#round(summary(fit1)[[1]][params,c("mean","2.5%","97.5%")],4)
#mle

tab=round(cbind(summary(fit1)[[1]][params,c("mean","2.5%","97.5%")],
      c(mle$beta,mle$sigmasq,mle$tausq,1/mle$phi)),4)
colnames(tab)[4]="mle"


## running stan latent model
set.seed(1)
fit1_latent <- stan(
  file = "bayesgp1_latent.stan",  # Stan program
  data = data3list,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  #open_progress = T,
  cores = 4             # number of cores (could use one per chain)
  #refresh = 10             # no shown every 100 iterations
  )

## save(file="../data/data3subset200_bayesgp1_latent.Rdata",list=c("fit1_latent"))
## load("../data/data3subset200_bayesgp1_latent.Rdata")

get_elapsed_time(fit1_latent)

### posterior estimates, traceplots, densities
plot(fit1_latent, plotfun = "trace",pars=params, inc_warmup = TRUE)
#plot(fit1_latent, plotfun = "dens",pars=params, inc_warmup = TRUE)
mcmc_dens_overlay(fit1_latent, pars = params)
plot(fit1_latent,pars=params,ci_level=0.95)

tab_latent=round(cbind(summary(fit1_latent)[[1]][params,c("mean","2.5%","97.5%")],
      c(mle$beta,mle$sigmasq,mle$tausq,1/mle$phi)),4)
colnames(tab_latent)[4]="mle"


### combining results from marginalized and latent model MCMC samples and MLE estimates
all_results=as.data.frame(rbind(summary(fit1)[[1]][params,c("mean","2.5%","97.5%")],
summary(fit1_latent)[[1]][params,c("mean","2.5%","97.5%")],
cbind(c(mle$beta,mle$sigmasq,mle$tausq,1/mle$phi),rep(NA,4),rep(NA,4))))
row.names(all_results)=NULL
all_results$param=rep(params,3)
all_results$method=c(rep("marginalized MCMC",5),
  rep("latent MCMC",5),rep("marginalized MLE",5))

ggplot(all_results,aes(x=param,y=mean,col=method,group=method)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`,color=method)) +
  ylab("Estimates") + xlab("Parameters")

### Rhat
round(summary(fit1_latent)[[1]][params,c("Rhat","n_eff","se_mean")],4)

#### predictions ####
data3=read.csv("../data/dataset3.csv")
data3out=data3[101:120,]
data3=data3[1:100,]
cross_dist=rdist(data3out[,c("sx","sy")],data3[,c("sx","sy")])
## stan data file
data3list_new=list(n=nrow(data3),y=data3$y,
  x=data3$x,dmat=as.matrix(dist(data3[,c("sx","sy")])),
  n_new=nrow(data3out),
  x_new=data3out$x,dmat_new=cross_dist)

## running latent model with predictions
set.seed(1)
fit1_latent_with_pred <- stan(
  file = "bayesgp1_latent_with_pred.stan",  # Stan program
  data = data3list_new,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  #open_progress = T,
  cores = 4             # number of cores (could use one per chain)
  #refresh = 10             # no shown every 100 iterations
  )

## save(file="../data/data3subset100_bayesgp1_latent_with_pred.Rdata",list=c("fit1_latent_with_pred"))
## load("../data/data3subset100_bayesgp1_latent_with_pred.Rdata")

pred_params=paste0('y_new[',1:20,"]")
pred_mat=as.data.frame(summary(fit1_latent_with_pred)[[1]][pred_params,c("mean","2.5%","97.5%")])
pred_mat$truth=data3out$y

ggplot(pred_mat,aes(x=truth,y=mean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`)) +
  ylab("True y") + xlab("Predicted y") + 
  geom_abline(intercept=0,slope=1,col="red")

log_lik_latent <- extract_log_lik(fit1_latent_with_pred)
waic_latent <- waic(log_lik_latent)
waic_latent$waic

### non-spatial model
data3list_ns=list(n=nrow(data3),y=data3$y,
  x=data3$x,n_new=nrow(data3out),x_new=data3out$x)

## running non-spatial linear model with predictions
set.seed(1)
fit1_non_spatial <- stan(
  file = "bayeslinearmodel.stan",  # Stan program
  data = data3list_ns,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  #open_progress = T,
  cores = 4             # number of cores (could use one per chain)
  #refresh = 10             # no shown every 100 iterations
  )

## save(file="../data/data3subset100_bayeslinearmodel.Rdata",list=c("fit1_non_spatial"))
## load("../data/data3subset100_bayeslinearmodel.Rdata")

pred_mat_non_spatial=as.data.frame(summary(fit1_non_spatial)[[1]][pred_params,c("mean","2.5%","97.5%")])
pred_mat_non_spatial$truth=data3out$y

ggplot(pred_mat_non_spatial,aes(x=truth,y=mean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`)) +
  ylab("True y") + xlab("Predicted y") + 
  geom_abline(intercept=0,slope=1,col="red")

log_lik_non_spatial <- extract_log_lik(fit1_non_spatial)
waic_non_spatial <- waic(log_lik_non_spatial)
waic_non_spatial$waic

##### spBayes package ######
## BEF data ###
## Data preliminaries
data(BEF.dat)
BEF.dat <- BEF.dat[BEF.dat$ALLBIO02_KGH>0,]
bio <- BEF.dat$ALLBIO02_KGH*0.001;
log.bio <- log(bio)
## Extract the coordinates
coords <- as.matrix(BEF.dat[,c("XUTM","YUTM")])

## Make a surface plot
x.res <- 100; y.res <- 100

surf <- mba.surf(cbind(coords, log.bio), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)")
points(coords)

### variogram on the residuals ###
BEF.dat$log.bio=log.bio
ols=lm(log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,data=BEF.dat)
BEF.dat$res=ols$residuals

max.dist <- 0.5*max(iDist(BEF.dat[,2:3]))
bins <- 20

vario <- variog(coords=BEF.dat[,2:3], data=BEF.dat$res, 
                 uvec=(seq(0, max.dist, length=bins)))

dev.new()
plot(vario,pch=16)

vfit <-variofit(vario, ini.cov.pars=c(0.04,3/6000), ##sigma^2 and 1/phi 
                cov.model="exponential", minimisation.function="optim",
                nugget=0.08, weights="equal")
vfit

p <- 6 ## This is the number of columns in the design matrix
## Set the prior mean and precision for the regression
beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

## For use with bayesGeostatExact, do the following
phi <- 1/vfit$cov.pars[2] ## Set the spatial range (from the variogram)
alpha <- vfit$nugget/vfit$cov.pars[1] ## Set the nugget/partial-sill ratio
sigma.sq.prior.shape <- 2.0 ## Set IG shape for sigma.sq (partial sill)
sigma.sq.prior.rate <- 0.08 ## Set IG scale for sigma.sq (partial sill)

## Run bayesGeostatExact to deliver exact posterior samples
set.seed(1)
sp.exact <- bayesGeostatExact(
    log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,
    data=BEF.dat, coords=coords, n.samples=1000,
    beta.prior.mean=beta.prior.mean,
    beta.prior.precision=beta.prior.precision,
    cov.model="exponential",
    phi=phi, alpha=alpha,
    sigma.sq.prior.shape=sigma.sq.prior.shape,
    sigma.sq.prior.rate=sigma.sq.prior.rate,
    sp.effects=FALSE)

##Produce the posterior summaries
round(summary(sp.exact$p.samples)$quantiles,3)

## Run spLM to deliver MCMC samples from marginal posterior distributions
n.samples <- 1000
set.seed(1)
bef.sp <- spLM(log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,
               data=BEF.dat, coords=coords, starting=list("phi"=3/200,"sigma.sq"=0.08,
                                                          "tau.sq"=0.02), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
               priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2, 0.08),
                           "tau.sq.IG"=c(2, 0.02)), cov.model="exponential",n.samples=n.samples)

round(summary(mcmc(bef.sp$p.theta.samples))$quantiles,3)

## Recover spatial residuals using spRecover
burn.in <- floor(0.75*n.samples)
bef.sp <- spRecover(bef.sp, start=burn.in, thin=2)

## The posterior samples of the regression coefficients and the spatial effects can then be obtained as
beta.samples = bef.sp$p.beta.recover.samples
w.samples = bef.sp$p.w.recover.samples

## Obtain trace plots for regression coefficients
dev.new()
par(mfrow=c(3,2))
plot(beta.samples, auto.layout=TRUE, density=FALSE)

round(summary(mcmc(bef.sp$p.beta.recover.samples))$quantiles,3)
round(summary(mcmc(bef.sp$p.theta.recover.samples))$quantiles,3)

## Obtain posterior means and sd's of spatial residuals for each location
w.hat.mu <- apply(w.samples,1,mean)
w.hat.sd <- apply(w.samples,1,sd)

## Plot the spatial residual mean surface and a map of sd's
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords, ols$residuals), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="LM residuals")
surf <- mba.surf(cbind(coords, w.hat.mu), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="Mean spatial effects (w(s))")

## loading shp files for predictions
BEF.shp <- readOGR("../data/BEF-data/BEF_bound.shp")
shp2poly <- BEF.shp@polygons[[1]]@Polygons[[1]]@coords
BEF.poly <- as.matrix(shp2poly)
BEF.grids <- readGDAL("../data/BEF-data/dem_slope_lolosptc_clip_60.img")

## Construct the prediction design matrix for the entire grid extent.
pred.covars <- cbind(BEF.grids[["band1"]], BEF.grids[["band2"]], BEF.grids[["band3"]], BEF.grids[["band4"]], BEF.grids[["band5"]])
pred.covars <- cbind(rep(1, nrow(pred.covars)), pred.covars)


## Extract the coordinates of the BEF bounding polygon vertices and use the pointsInPoly (spBayes) function to obtain the desired subset of the prediction design matrix and associated prediction coordinates (i.e., pixel centroids).
pred.coords <- SpatialPoints(BEF.grids)@coords
pointsInPolyOut <- pointsInPoly(BEF.poly, pred.coords)
pred.covars <- pred.covars[pointsInPolyOut,]
pred.coords <- pred.coords[pointsInPolyOut,]

bef.bio.pred <- spPredict(bef.sp, start=burn.in, thin=2, pred.coords=pred.coords, pred.covars=pred.covars)

## Mapping the predicted values
bef.bio.pred.mu = apply(bef.bio.pred$p.y.predictive.samples,1,mean)
bef.bio.pred.sd = apply(bef.bio.pred$p.y.predictive.samples,1,sd)

surf <- mba.surf(cbind(coords, log.bio), no.X=x.res, no.Y=x.res, extend=TRUE, sp=TRUE)$xyz.est
#surf <- surf [!is.na(over(surf, BEF.shp)),]
surf <- surf [!is.na((over(surf, BEF.shp)))[,1],]
surf <- as.image.SpatialGridDataFrame(surf)
z.lim <- range(surf[["z"]], na.rm=TRUE)

pred.grid <- as.data.frame(list(pred.coords, pred.mu=bef.bio.pred.mu, pred.sd=bef.bio.pred.sd))
coordinates(pred.grid) = c("x", "y")
gridded(pred.grid) <- TRUE
pred.mu.image <- as.image.SpatialGridDataFrame(pred.grid["pred.mu"])

par(mfrow=c(1,2))
image.plot(surf, axes=TRUE, zlim=z.lim, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Log metric tons of biomass")
plot(BEF.shp, add=TRUE)
image.plot(pred.mu.image, zlim=z.lim, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Mean predicted log metric tons of biomass")
plot(BEF.shp, add=TRUE)

