##Predictive process example
rm(list=ls())
library(MBA)
library(fields)
library(spBayes)
library(spNNGP)

## first part same as last class (full GP and predictive process)
##Simulated some data
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)

n <- 500
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))
p <- length(B)

sigma.sq <- 2
tau.sq <- 0.1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X%*%B + w, sqrt(tau.sq))

##Set up spLM call
n.samples <- 2000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)

tuning <- list("phi"=0.05, "sigma.sq"=0.05, "tau.sq"=0.05)

priors <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
               "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
               "tau.sq.IG"=c(2, 0.1))

cov.model <- "exponential"

n.report <- 100
verbose <- TRUE

##Call full GP and predictive process GP mosels
burn.in <- floor(0.75*n.samples)

## loading results of full GP and predictive process runs from last lecture
## load("../data/pred_proc.Rdata")

## NNGP (m=15) 
m.s.5 <- spNNGP(y~X-1, coords=coords, starting=starting, method="latent", n.neighbors=5,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, return.neighbor.info = TRUE, n.omp.threads=2)


m.s.15 <- spNNGP(y~X-1, coords=coords, starting=starting, method="latent", n.neighbors=15,
                 tuning=tuning, priors=priors, cov.model=cov.model,
                 n.samples=n.samples, return.neighbor.info = TRUE, n.omp.threads=2)


## Timing
m.gp$run.time
m.pp.gp.25$run.time
m.pp.gp.64$run.time
m.s.5$run.time
m.s.15$run.time

## WAIC and other model comparison metrics for the NNGP models
spDiag(m.s.5)$WAIC
spDiag(m.s.15)$WAIC

## Summary cov parameters
round(summary(m.gp$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.pp.gp.25$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.pp.gp.64$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.s.5$p.theta.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.s.15$p.theta.samples)$quantiles[,c(3,1,5)],2)


## Summary random effects
m.gp.w.hat <- apply(m.gp$p.w.recover.samples, 1, median)

m.pp.25.w.hat <- apply(m.pp.gp.25$p.w.recover.samples, 1, median)

m.pp.64.w.hat <- apply(m.pp.gp.64$p.w.recover.samples, 1, median)

m.s.15.w.hat <- apply(m.s.15$p.w.samples, 1, median)

m.s.5.w.hat <- apply(m.s.5$p.w.samples, 1, median)

## Interpolate
surf.w <- mba.surf(cbind(coords, w), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.gp <- mba.surf(cbind(coords, m.gp.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.pp.25 <- mba.surf(cbind(coords, m.pp.25.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.pp.64 <- mba.surf(cbind(coords, m.pp.64.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.nngp.15 <- mba.surf(cbind(coords, m.s.15.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.nngp.5 <- mba.surf(cbind(coords, m.s.5.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est

dev.new()
par(mfrow=c(2,3))
image.plot(surf.w, main="True w")
image.plot(surf.gp, main="GP estimated w")
image.plot(surf.pp.25, main="PPGP knots 25 w"); points(m.pp.gp.25$knot.coords, pch=19)
image.plot(surf.pp.64, main="PPGP knots 64 w"); points(m.pp.gp.64$knot.coords, pch=19)
image.plot(surf.nngp.5, main="NNGP (m=5) estimated w")
image.plot(surf.nngp.15, main="NNGP (m=15) estimated w")

dev.new()
par(mfrow=c(2,3))
plot(w, m.gp.w.hat,  xlab="True w", ylab="Full GP Posterior median w")
plot(w, m.pp.25.w.hat,  xlab="True w", ylab="PPGP 36 Posterior median w")
plot(w, m.pp.64.w.hat,  xlab="True w", ylab="PPGP 64 Posterior median w")
plot(w, m.s.15.w.hat,  xlab="True w", ylab="NNGP 5 Posterior median w")
plot(w, m.s.15.w.hat,  xlab="True w", ylab="NNGP 15 Posterior median w")
plot(m.gp.w.hat, m.s.15.w.hat,  xlab="Full GP Posterior median w", ylab="NNGP 15 Posterior median w")


### BCEF dataset: A forest canopy height dataset at n=188,717 locations
## code from https://arxiv.org/pdf/2001.09111.pdf
## we will analyze a subset of the data
data(BCEF)

###
set.seed(1)
in.sample=sample(1:nrow(BCEF),10000)
BCEF.mod <- BCEF[in.sample,]
BCEF.pred <- BCEF[-in.sample,]

## code for latent model 
n.samples <- 5000
starting <- list("phi"=3/2, "sigma.sq"=40, "tau.sq"=1)
priors <- list("phi.Unif"=c(3/10, 3/0.1), "sigma.sq.IG"=c(2, 40),
"tau.sq.IG"=c(2, 10))
cov.model <- "exponential"
tuning <- list("phi"=0.02)
bcef.s <- spNNGP(FCH~PTC, coords=c("x","y"), data=BCEF.mod, starting=starting,
method="latent", n.neighbors=10,
tuning=tuning, priors=priors, cov.model=cov.model,
n.samples=n.samples, n.omp.threads=18, n.report=2500,
fit.rep=TRUE, sub.sample=list(start=4000, thin=10))

### response model
tuning <- list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.005)
bcef.r <- spNNGP(FCH~PTC, coords=c("x","y"), data=BCEF.mod, starting=starting,
method="response", n.neighbors=10,
tuning=tuning, priors=priors, cov.model=cov.model,
n.samples=n.samples, n.omp.threads=18, n.report=2500,
fit.rep=TRUE, sub.sample=list(start=4000, thin=10),
verbose=FALSE)

### conjugate model
theta.alpha <- as.matrix(expand.grid(seq(0.1,1,length.out=15),
seq(3/10,3/0.1,length.out=15)))
colnames(theta.alpha) <- c("alpha","phi")
bcef.c <- spConjNNGP(FCH~PTC, coords=c("x","y"), data=BCEF.mod,
cov.model="exponential", sigma.sq.IG=c(2, 40),
n.neighbors=10,theta.alpha=theta.alpha,
k.fold = 2, score.rule = "crps",
fit.rep=TRUE, n.samples=200,
n.omp.threads=18,
verbose=FALSE)

## load("../data/nngp.Rdata")

## run times
bcef.s$run.time
bcef.r$run.time
bcef.c$run.time

### predictions 
BCEF.pred.mod=BCEF[sample(1:nrow(BCEF.pred),2000),]
BCEF.pred.s <- predict(bcef.s, X.0 = as.matrix(cbind(1,BCEF.pred.mod[,c("PTC")])),
  coords.0=as.matrix(BCEF.pred.mod[,c("x","y")]))
BCEF.pred.r <- predict(bcef.r, X.0 = as.matrix(cbind(1,BCEF.pred.mod[,c("PTC")])),
  coords.0=as.matrix(BCEF.pred.mod[,c("x","y")]))
BCEF.pred.c <- predict(bcef.c, X.0 = as.matrix(cbind(1,BCEF.pred.mod[,c("PTC")])),
  coords.0=as.matrix(BCEF.pred.mod[,c("x","y")]))

###### comparing the predictions
dev.new()
predmat=data.frame(true=BCEF.pred.mod$FCH,latent=rowMeans(BCEF.pred.s$p.y.0),
response=rowMeans(BCEF.pred.r$p.y.0),conjugate=as.vector(BCEF.pred.c$y.0.hat))

library(GGally)
ggpairs(predmat)

### analyzing the same data using BRISC 
library(BRISC)
coords=as.matrix(BCEF.mod[,c("x","y")])
x=cbind(1,as.matrix(BCEF.mod[,c("PTC")]))
y=as.vector(BCEF.mod[,c("FCH")])
estimation_result <- BRISC_estimation(coords, y, x)
bootstrap_result <- BRISC_bootstrap(estimation_result, n_boot = 100)

## save(file="../data/nngp.Rdata",list=c("bcef.s","bcef.r","bcef.c","BCEF.pred.s","BCEF.pred.r","BCEF.pred.c","estimation_result","bootstrap_result"))

### comparing parameter estimates ###
library(tidyverse)
estimates=as.data.frame(rbind(colMeans(bcef.s$p.beta.samples),
colMeans(bcef.r$p.beta.samples),
bcef.c$beta.hat,
estimation_result$Beta)) %>% 
  pivot_longer(everything(),names_to = "params",values_to = "estimate") %>%
  as.data.frame()

low=as.data.frame(rbind(apply(bcef.s$p.beta.samples,2,quantile,0.025),
apply(bcef.r$p.beta.samples,2,quantile,0.025),
bcef.c$beta.hat - 1.96*sqrt(diag(bcef.c$beta.var)),
bootstrap_result$confidence.interval[1,1:2]))  %>% 
  pivot_longer(everything(),names_to = "params",values_to = "low") %>%
  as.data.frame()

high=as.data.frame(rbind(apply(bcef.s$p.beta.samples,2,quantile,0.975),
apply(bcef.r$p.beta.samples,2,quantile,0.975),
bcef.c$beta.hat + 1.96*sqrt(diag(bcef.c$beta.var)),
bootstrap_result$confidence.interval[2,1:2]))  %>% 
  pivot_longer(everything(),names_to = "params",values_to = "high") %>%
  as.data.frame()

param_df=cbind(estimates,low=low$low,high=high$high)
param_df$method=unlist(lapply(c("spNNGP.latent","spNNGP.response","spNNGP.conjugate","BRISC"),rep,2))

param_df %>% ggplot(aes(x=1,y=estimate,fill=method)) +
 facet_wrap(. ~ params,scales = "free_y") +
 geom_col(position = position_dodge(width=0.9)) +
 geom_errorbar(aes(ymin = low, ymax = high), position = position_dodge(width=0.9), width = 0.5)


### large scale simulation of Gaussian Process realizations using BRISC package

mycorrplot=function(M,t="",o="viridis",d=1){
  #dev.new()
  melted_cormat <- melt(M)
  ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_fill_viridis_c(option=o,direction=d) + 
  ylab("") + xlab("") +
  ggtitle(t)
  }


### simulation using Matern NNGP for heatmaps
set.seed(1)
n <- 1000
coords <- cbind(runif(n,0,1), runif(n,0,1))
coords=coords[order(rowSums(coords)),]

sigma.sq = 1
phi = 5
nu = 1/2

set.seed(1)
nsim = 10000
nei=100
simulation_result <- BRISC_simulation(coords, sigma.sq = sigma.sq, phi = phi, nu=nu, sim_number = nsim, n.neighbors = nei)
M=cov(t(simulation_result$output.data))

dmat=as.matrix(dist(coords))
Mtrue=Matern(dmat,alpha=phi,smoothness=nu,phi=sigma.sq)

w_GP=rmvnorm(nsim,rep(0,n),Mtrue)
Mtrue_est=cov(w_GP)


#par(mfrow=c(1,3))
p1=mycorrplot(M,t=paste0("NNGP: Matern, nu=",nu,", phi=",phi))
p5=mycorrplot(Mtrue_est,t=paste0("GP: Matern, nu=",nu,", phi=",phi))
dev.new()
plot(p1)
dev.new()
plot(p5)

diff_exp=M-Mtrue_est

### same sim with Matern 3/2
phi = 5
nu = 3/2 

set.seed(1)
nsim = 10000
nei=100
simulation_result <- BRISC_simulation(coords, sigma.sq = sigma.sq, phi = phi, nu=nu, sim_number = nsim, n.neighbors = nei, cov.model = "matern")
M=cov(t(simulation_result$output.data))

dmat=as.matrix(dist(coords))
Mtrue=Matern(dmat,alpha=phi,smoothness=nu,phi=sigma.sq)

w_GP=rmvnorm(nsim,rep(0,n),Mtrue)
Mtrue_est=cov(w_GP)

max(abs(M - Mtrue))
mean(abs(M - Mtrue))
mean(abs(M - Mtrue_est))

dev.new()
#par(mfrow=c(1,3))
p1=mycorrplot(M,t=paste0("NNGP: Matern, nu=",nu,", phi=",phi))
p5=mycorrplot(Mtrue_est,t=paste0("GP: Matern, nu=",nu,", phi=",phi))
plot(p1)
plot(p5)

diff_mat=M-Mtrue_est

plot(density(diff_mat),xlab="Difference",main="",col="red")
lines(density(diff_exp),col="blue")
abline(v=0)
legend("topleft",legend=c("exponential","Matern 3/2"),col=c("blue","red"),lty=1)

