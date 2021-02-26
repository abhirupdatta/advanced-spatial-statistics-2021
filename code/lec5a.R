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
library(spNNGP)
library(rstan)
library(wesanderson)

#### spatial model with a non-Gaussian outcome ####
#### Eastern hemlock (Tsuga canadensis L.) analyzed in Lany et al. (2019)
data("MI_TSCA")

### subset of the data for illustration
sub_data_large =MI_TSCA %>% 
  dplyr::filter(long > quantile(MI_TSCA$long,0.5),lat < quantile(MI_TSCA$lat,0.5))
presence=which(sub_data_large$TSCA==1)
absence=which(sub_data_large$TSCA==0)
set.seed(1)
index=c(sample(presence,50),sample(absence,150))
sub_data=sub_data_large[index,]

### plot of the locations
presence_sub=which(sub_data$TSCA==1)
plot(MI_TSCA$long,MI_TSCA$lat,col="grey",xlab="long",ylab="lat")
points(sub_data$long[presence_sub],sub_data$lat[presence_sub],col="red")
points(sub_data$long[-presence_sub],sub_data$lat[-presence_sub],col="blue")

## stan data file
data_binom_small=list(n=nrow(sub_data),p=ncol(sub_data[,-(1:3)]),
  y=sub_data$TSCA,x=as.matrix(sub_data[,-(1:3)]),
  dmat=as.matrix(dist(as.matrix(sub_data[,c(2,1)]))))



## running stan for non-spatial model
set.seed(1)
fit_binom_small_nonspatial <- stan(
  file = "bayeslinearmodel_binomial.stan",  # Stan program
  data = data_binom_small,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  #open_progress = T,
  cores = 4             # number of cores (could use one per chain)
  #refresh = 10             # no shown every 100 iterations
  )

## save(file="../data/data_binom_small_nonspatial_bayesgp1.Rdata",list=c("fit_binom_small_nonspatial"))
## load("../data/data_binom_small_nonspatial_bayesgp1.Rdata")

## run time
get_elapsed_time(fit_binom_small_nonspatial)

### posterior estimates, traceplots, densities
params=c("alpha",paste0("beta[",1:6,"]"))
plot(fit_binom_small_nonspatial, plotfun = "trace",pars=params, inc_warmup = TRUE)
plot(fit_binom_small_nonspatial,pars=params,ci_level=0.95)
tab_non_spatial=round(summary(fit_binom_small_nonspatial)[[1]][params,c("mean","2.5%","97.5%")],4)

### Rhat
round(summary(fit_binom_small_nonspatial)[[1]][params,c("Rhat","n_eff","se_mean")],4)



### waic
log_lik_latent_nonspatial <- extract_log_lik(fit_binom_small_nonspatial)
waic_latent_nonspatial <- waic(log_lik_latent_nonspatial)
waic_latent_nonspatial$waic

## running stan for Binomial GP model
beta_init=round(summary(fit_binom_small_nonspatial)[[1]][params,c("mean")],4)
names(beta_init)=NULL
init_beta <- function(chain_id) {
  list(alpha= beta_init[1], beta = beta_init[-1])
}


set.seed(1)
fit_binom_small <- stan(
  file = "bayesgp1_binomial.stan",  # Stan program
  data = data_binom_small,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  #open_progress = T,
  cores = 4,             # number of cores (could use one per chain)
  init = init_beta
  #refresh = 10             # no shown every 100 iterations
  )

## save(file="../data/data_binom_small_bayesgp1.Rdata",list=c("fit_binom_small"))
## load("../data/data_binom_small_bayesgp1.Rdata")

## run time
get_elapsed_time(fit_binom_small)

### posterior estimates, traceplots, densities
params=c("alpha",paste0("beta[",1:6,"]"),"sigmasq","phi")
plot(fit_binom_small, plotfun = "trace",pars=params, inc_warmup = TRUE)
plot(fit_binom_small,pars=params,ci_level=0.95)
tab=round(summary(fit_binom_small)[[1]][params,c("mean","2.5%","97.5%")],4)

### Rhat
round(summary(fit_binom_small)[[1]][params,c("Rhat","n_eff","se_mean")],4)

### waic
log_lik_latent <- extract_log_lik(fit_binom_small)
waic_latent <- waic(log_lik_latent)
waic_latent$waic

### comparing the parameter estimates ###
### combining results from spatial and non-spatial model 
all_results=as.data.frame(rbind(tab_non_spatial,tab[1:7,]))
row.names(all_results)=NULL
all_results$param=rep(head(params,-2),2)
all_results$method=c(rep("non-spatial",7),rep("spatial GP",7))

ggplot(all_results,aes(x=param,y=mean,col=method,group=method)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`,color=method)) +
  ylab("Estimates") + xlab("Parameters")



#### spatially varying coefficient model ####
### code taken from https://www.sciencedirect.com/science/article/abs/pii/S1364815219310412
library(spBayes)
data(PM10.dat)

### data with observed PM10 and CTM output
PM10.mod <- PM10.dat[!is.na(PM10.dat$pm10.obs),]

### data only with CTM output where we will predict PM10
PM10.pred <- PM10.dat[is.na(PM10.dat$pm10.obs),]

#### plot of locations with regulatory PM10 data and CTM grid
PM10.dat %>% 
  ggplot(aes(x=x.coord,y=y.coord)) + 
  geom_point(shape=1,col="grey") +
  geom_point(data=PM10.mod,aes(x=x.coord,y=y.coord),col="red")

### plot of regulatory data
pal <- wes_palette("Zissou1", 100, type = "continuous")
PM10.dat %>% 
  ggplot(aes(x=x.coord,y=y.coord)) + 
  geom_point(shape=1,col="grey") +
  geom_point(data=PM10.mod %>% mutate(pm10.obs=pm10.obs^2),
             aes(x=x.coord,y=y.coord,col=pm10.obs)) +
  scale_color_gradientn(colours = pal)

### plot of CTM output
PM10.pred %>% mutate(pm10.ctm=pm10.ctm^2) %>%
  ggplot(aes(x=x.coord,y=y.coord,col=pm10.ctm)) + 
  geom_point() +
  scale_color_gradientn(colours = pal)

### scatter plots of true and CTM predicted PM10 (square-root scale)
PM10.mod %>% 
  ggplot(aes(x=pm10.obs,y=pm10.ctm,col=x.coord)) + 
  geom_point() +
  scale_color_viridis_c() +
  geom_abline(intercept=0,slope=1,col="red") 
  
PM10.mod %>% 
  ggplot(aes(x=pm10.obs,y=pm10.ctm,col=y.coord)) + 
  geom_point() +
  scale_color_viridis_c() +
  geom_abline(intercept=0,slope=1,col="red") 

### setting priors and initial values
d.max <- max(iDist(PM10.mod[,c("x.coord","y.coord")]))
d.max #km

r <- 2
priors <- list("phi.Unif"=list(rep(3/(0.75*d.max), r), rep(3/(0.001*d.max), r)),
"sigma.sq.IG"=list(rep(2, r), rep(1, r)),
"tau.sq.IG"=c(2, 1))

starting <- list("phi"=rep(3/(0.1*d.max), r), "sigma.sq"=rep(1, r), "tau.sq"=1)

tuning <- list("phi"=rep(0.1, r), "sigma.sq"=rep(0.05, r), "tau.sq"=0.1)
n.samples <- 10000

### running MCMC for a spatially varying coefficient model
m.3 <- spSVC(pm10.obs ~ pm10.ctm, coords=c("x.coord","y.coord"),
data=PM10.mod, starting=starting, svc.cols=c(1,2),
tuning=tuning, priors=priors, cov.model="exponential",
n.samples=n.samples, n.report=5000, n.omp.threads=4)

### recoveering the parameters using composition sampling
m.3 <- spRecover(m.3, start=floor(0.75*n.samples), thin=2,
n.omp.threads=4, verbose=FALSE)

### parameter summary
round(summary(m.3$p.beta.recover.samples)$quantiles[,c(3,1,5)],3)

round(summary(m.3$p.theta.recover.samples)$quantiles[,c(3,1,5)],3)


tilde.beta.0 <- apply(
m.3$p.tilde.beta.recover.samples[["tilde.beta.(Intercept)"]],
1, median)

tilde.beta.ctm <- apply(
m.3$p.tilde.beta.recover.samples[["tilde.beta.pm10.ctm"]],
1, median)

PM10.mod$intercept=tilde.beta.0
PM10.mod$slope=tilde.beta.ctm

## plotting the spatially varying intercept
PM10.dat %>% 
  ggplot(aes(x=x.coord,y=y.coord)) + 
  geom_point(shape=1,col="grey") +
  geom_point(data=PM10.mod,aes(x=x.coord,y=y.coord,col=intercept),size=2) +
  scale_color_gradientn(colours = pal)

## plotting the spatially varying slope
PM10.dat %>% 
  ggplot(aes(x=x.coord,y=y.coord)) + 
  geom_point(shape=1,col="grey") +
  geom_point(data=PM10.mod,aes(x=x.coord,y=y.coord,col=slope),size=2) +
  scale_color_gradientn(colours = pal)

### predictions at all the CTM grid locations
m.3.pred <- spPredict(m.3, pred.covars=cbind(1, PM10.pred$pm10.ctm),
pred.coords=PM10.pred[,1:2], thin=25,
joint=F, n.omp.threads=4, verbose=FALSE)

## save(file="../data/pm10_svc.Rdata",list=c("m.3","tilde.beta.0","tilde.beta.ctm","m.3.pred"))
## load("../data/pm10_svc.Rdata")

PM10.pred$pm10.pred=rowMeans(m.3.pred$p.y.predictive.samples^2)
PM10.pred$pm10.pred.exceed.prob=rowMeans(m.3.pred$p.y.predictive.samples^2>50)

## plotting the predicted PM10
PM10.pred %>% 
  ggplot(aes(x=x.coord,y=y.coord,col=pm10.pred)) + 
  geom_point(shape=16) +
  #geom_point(data=PM10.mod,aes(x=x.coord,y=y.coord,col=slope),size=2) +
  scale_color_gradientn(colours = pal)

## plotting the predicted probability of PM10 exceeding 50 mu g / m^3
PM10.pred %>% 
  ggplot(aes(x=x.coord,y=y.coord,col=pm10.pred.exceed.prob)) + 
  geom_point(shape=16) +
  #geom_point(data=PM10.mod,aes(x=x.coord,y=y.coord,col=slope),size=2) +
  scale_color_gradientn(colours = pal)
