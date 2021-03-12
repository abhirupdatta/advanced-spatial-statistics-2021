library(rstan)
library(mvtnorm)
library(tidyverse)
library(loo)
library(sf)
library(reshape2)
library(ggrepel)
library(spdep)

# spatial analyses - for later
library(raster) # raster
library(USAboundaries) 

# color scales
library(viridis)

#### chloropleth maps
#### Covid cases in Maryland counties
mdcovid=read.csv("../data/mdCovid.csv")

md.state<-us_states(resolution = 'high', states='maryland')
glimpse(md.state) # look at spatial object - basically a dataframe
md.cos<-us_counties(resolution = 'high', states='maryland')
### separating Baltimoree city and county which has the same name "Baltimore"
md.cos$name[which(md.cos$geoid=="24005")]="Baltimore County"
md.cos$name[which(md.cos$geoid=="24510")]="Baltimore City"
md.map<-ggplot(md.cos)+
  geom_sf(fill='blue', color='white')
md.cos$lon<-st_coordinates(st_centroid(md.cos))[,1] # add longitude to sf
md.cos$lat<-st_coordinates(st_centroid(md.cos))[,2] # add latitude to sf
md.map + geom_label_repel(data=md.cos, aes(x=lon, y=lat, label=name))

### merging the sf object with data
md.cos.covid=md.cos %>% 
  mutate(County=name) %>% 
  left_join(mdcovid) %>%
  mutate(Cases=as.numeric(gsub(",", "", Cases)),
         Deaths=as.numeric(gsub(",", "", Deaths)))

### cases maps
ggplot(md.cos.covid,aes(fill=Cases))+
  geom_sf(color='white') +
  scale_fill_viridis(alpha=.5) +
  geom_label_repel(data=md.cos.covid, aes(x=lon, y=lat, label=name),col="black",fill=NA,alpha=1)

### deaths maps
ggplot(md.cos.covid,aes(fill=Deaths))+
  geom_sf(color='white') +
  scale_fill_viridis(alpha=.5) +
  geom_label_repel(data=md.cos.covid, aes(x=lon, y=lat, label=name),col="black",fill=NA,alpha=1)


### Moran's I 
md.nei <- poly2nb(md.cos.covid)
md.nei.list <- nb2listw(md.nei)
moran.test(md.cos.covid$Cases,md.nei.list)
moran.test(md.cos.covid$Deaths,md.nei.list)

#### Modeling areal datasets using DAGAR and ICAR ####
### simulations 

#################################
#####   Datasets (2d grid)  #####
#################################

### produces adjacency graph for a 2d grid
inclattice=function(m){
  n=m^2
  Minc=matrix(0,n,n)
  for(i in 1:(m-1))	for(j in 1:(m-1)) Minc[(i-1)*m+j,(i-1)*m+j+1]=Minc[(i-1)*m+j,i*m+j]=1
  for(i in 1:(m-1)) Minc[(i-1)*m+m,i*m+m]=1
  for(j in 1:(m-1)) Minc[(m-1)*m+j,(m-1)*m+j+1]=1
  Minc+t(Minc)
}

### produces vector with list of neighbors
neighbors=function(Minc){
  n=nrow(Minc)
  unlist(lapply(2:n,function(i) which(Minc[i,1:(i-1)]==1)))
}

### produces number of directed neighbors and a vector where each the adjacency of each node ends in the nei vector
adj_index=function(Minc){
  Minc.low=Minc
  Minc.low[upper.tri(Minc)]=0
  list(N_nei=rowSums(Minc.low),adj_index=cumsum(rowSums(Minc.low)))
  }

m=20
Minc=inclattice(m)
N=nrow(Minc) ## number of nodes
N_edges=sum(Minc)/2 ### number of edges
nei=neighbors(Minc) 
adj.ends=adj_index(Minc)$adj_index 
N_nei=adj_index(Minc)$N_nei

s=cbind(rep(1:m,m),kronecker(1:m,rep(1,m)))
dmat=as.matrix(dist(s))  ## the scaling by 6 is done to keep the average intersite distance similar to the 1-d path graph with m^2 vertices

phi=0.5 ### smooth settng 
#phi=2   ### rough setting
set.seed(1)
w=as.vector(rmvnorm(1,rep(0,N),exp(-phi*dmat)))
y=rpois(rep(1,N),exp(w))

plot_tab=data.frame(sx=s[,1],sy=s[,2],true=exp(w))
plot_tab %>% ggplot(aes(x=sx,y=sy,fill=true)) + 
  geom_tile() +
  scale_fill_viridis_c(direction = -1)

### fitting dagar model
dagar_dat <- list(N = N, N_edges=N_edges, nei=nei, 
                  adjacency_ends=adj.ends, y = y, N_nei=N_nei)

## load("../data/grid_smooth.Rdata")
## load("../data/grid_rough.Rdata")
fit_dagar_grid <- stan(file = 'dagar_poisson.stan', 
            data = dagar_dat, 
            iter=2000, 
            cores = 4)
mat_dagar <- as.matrix(fit_dagar_grid)
apply(mat_dagar[,1:2],2,summary)

params=c("rho","std_dev_w")
plot(fit_dagar_grid, plotfun = "trace",pars=params, inc_warmup = TRUE)
plot(fit_dagar_grid,pars=params,ci_level=0.95)

diag.mat_dagar=round(summary(fit_dagar_grid)[[1]][paste0("w[",1:N,"]"),c("Rhat","n_eff","se_mean")],4)
plot(diag.mat_dagar[,"Rhat"],ylim=c(0.95,2),ylab="Rhat",xlab="w_i")
abline(h=1.05,col="red")

what_dagar=colMeans(mat_dagar[,paste0("w[",1:N,"]")])
plot(exp(w),exp(what_dagar))
abline(a=0,b=1,col="red")

plot_tab$dagar=exp(what_dagar)

### model metrics
### rmse
dagar_rmse=sqrt(mean((w-what_dagar)^2))
### waic
log_lik_latent_dagar <- extract_log_lik(fit_dagar_grid)
waic_latent_dagar <- waic(log_lik_latent_dagar)

### fitting icar model
node1=unlist(lapply(2:N,function(i) rep(i,N_nei[i])))
node2=nei

icar_dat <- list(N = N, node1=node1, node2=node2,y = y, N_edges=N_edges)

fit_icar_grid <- stan(file = 'icar_poisson.stan', 
            data = icar_dat, 
            iter=2000, 
            cores = 4)
mat_icar <- as.matrix(fit_icar_grid)
apply(mat_icar[,1:2],2,summary)

params=c("std_dev_w")
plot(fit_icar_grid, plotfun = "trace",pars=params, inc_warmup = TRUE)
plot(fit_icar_grid,pars=params,ci_level=0.95)

diag.mat_icar=round(summary(fit_icar_grid)[[1]][paste0("w[",1:N,"]"),c("Rhat","n_eff","se_mean")],4)
plot(diag.mat_icar[,"Rhat"],ylim=c(0.95,2),ylab="Rhat",xlab="w_i")
abline(h=1.05,col="red")

what_icar=colMeans(mat_icar[,paste0("w[",1:N,"]")])
plot(exp(w),exp(what_icar))
abline(a=0,b=1,col="red")


### model 
icar_rmse=sqrt(mean((w-what_icar)^2))
### waic
log_lik_latent_icar <- extract_log_lik(fit_icar_grid)
waic_latent_icar <- waic(log_lik_latent_icar)

plot_tab$icar=exp(what_icar)
plot_tab %>% 
  pivot_longer(cols=!c(sx,sy),names_to="model",values_to="fit") %>%
  ggplot(aes(x=sx,y=sy,fill=fit)) + 
  geom_tile() +
  scale_fill_viridis_c(direction = -1) +
  facet_wrap(. ~ model)

data.frame(model=c("dagar","icar"),rmse=round(c(dagar_rmse,icar_rmse),2),
waic=round(c(waic_latent_dagar$waic,waic_latent_icar$waic),2))

## save(file="../data/grid_smooth.Rdata",list=c("fit_dagar_grid","fit_icar_grid"))
## save(file="../data/grid_rough.Rdata",list=c("fit_dagar_grid","fit_icar_grid"))




##### Slovenia stomach cancer analysis #####
# install.packages("remotes")
# remotes::install_github("DouglasMesquita/RASCO")
library(RASCO)
data(slovenia)

# plot(slovenia["E"])
# plot(slovenia["SE"])
# plot(slovenia["O"])

slovenia$ObyE <- cut(slovenia$O/slovenia$E, 
                 breaks = c(-0.1, 0.75, 1, 1.25, 1.5, 4.1), 
                 labels = c('0-0.75','0.75-1','1-1.25','1.25-1.5','1.5-4'))

dev.new()
ggplot(data = slovenia) +
    geom_sf() +
    geom_sf(data = slovenia, aes(fill = ObyE)) +
    scale_fill_viridis_d(direction=-1) +
    labs(fill="O/E")

dev.new()
ggplot(data = slovenia) +
    geom_sf() +
    geom_sf(data = slovenia, aes(fill = SE)) +
    scale_fill_viridis_c(direction=-1)

### obtaining neighbors 
slov.nei <- poly2nb(slovenia)  ## gets the border-based neighbors
coords <- st_coordinates(st_centroid(st_geometry(slovenia)))
plot(st_geometry(slovenia), border="grey")
plot(slov.nei, coords, pch = 19, cex = 0.4, add=TRUE)
N=length(slov.nei)
Minc=matrix(0,N,N)
for(i in 1:N) Minc[i,slov.nei[[i]]]=1

### DAGAR model

N=nrow(Minc) ## number of nodes
N_edges=sum(Minc)/2 ### number of edges
nei=neighbors(Minc) 
adj.ends=adj_index(Minc)$adj_index 
N_nei=adj_index(Minc)$N_nei 


dagar_dat <- list(N = N, N_edges=N_edges, nei=nei, 
                  adjacency_ends=adj.ends, y = slovenia$O,
                  offset=slovenia$E,X=slovenia$SEc,N_nei=N_nei)
## load(file="../data/slov.Rdata")

fit_dagar_slov <- stan(file = 'dagar_poisson_offset_covariates.stan', 
            data = dagar_dat, 
            iter=2000, 
            cores = 4)
mat_dagar <- as.matrix(fit_dagar_slov)

params=c("alpha","beta","rho","std_dev_w")
plot(fit_dagar_slov, plotfun = "trace",pars=params, inc_warmup = TRUE)
plot(fit_dagar_slov,pars=params,ci_level=0.95)

diag.mat_dagar=round(summary(fit_dagar_slov)[[1]][paste0("w[",1:N,"]"),c("Rhat","n_eff","se_mean")],4)
plot(diag.mat_dagar[,"Rhat"],ylim=c(0.95,2),ylab="Rhat",xlab="w_i")
abline(h=1.05,col="red")

### waic
log_lik_latent_dagar <- extract_log_lik(fit_dagar_slov)
waic_latent_dagar <- waic(log_lik_latent_dagar)


### non-spatial model

N=nrow(Minc) ## number of nodes

ns_dat <- list(N = N, y = slovenia$O,
                  offset=slovenia$E,X=slovenia$SEc)

fit_ns_slov <- stan(file = 'nonspatial_poisson_offset_covariates.stan', 
            data = ns_dat, 
            iter=2000, 
            cores = 4)
mat_ns <- as.matrix(fit_ns_slov)

params=c("alpha","beta")
plot(fit_ns_slov, plotfun = "trace",pars=params, inc_warmup = TRUE)
plot(fit_ns_slov,pars=params,ci_level=0.95)

### waic
log_lik_latent_ns <- extract_log_lik(fit_ns_slov)
waic_latent_ns <- waic(log_lik_latent_ns)


### ICAR model
node1=unlist(lapply(2:N,function(i) rep(i,N_nei[i])))
node2=nei

icar_dat <- list(N = N, N_edges=N_edges, node1=node1, node2=node2, y = slovenia$O,
                  offset=slovenia$E,X=slovenia$SEc)

fit_icar_slov <- stan(file = 'icar_poisson_offset_covariates.stan', 
            data = icar_dat, 
            iter=2000, 
            cores = 4)
mat_icar <- as.matrix(fit_icar_slov)

params=c("alpha","beta","std_dev_w")
plot(fit_icar_slov, plotfun = "trace",pars=params, inc_warmup = TRUE)
plot(fit_icar_slov,pars=params,ci_level=0.95)

diag.mat_icar=round(summary(fit_icar_slov)[[1]][paste0("w[",1:N,"]"),c("Rhat","n_eff","se_mean")],4)
plot(diag.mat_icar[,"Rhat"],ylim=c(0.95,2),ylab="Rhat",xlab="w_i")
abline(h=1.05,col="red")

### waic
log_lik_latent_icar <- extract_log_lik(fit_icar_slov)
waic_latent_icar <- waic(log_lik_latent_icar)

tab1=data.frame(model=c("Non-spatial","DAGAR","ICAR"),
  waic=round(c(waic_latent_ns$waic,waic_latent_dagar$waic,waic_latent_icar$waic),2))

## save(file="../data/slov.Rdata",list=c("fit_ns_slov","fit_dagar_slov","fit_icar_slov"))

#### comparing coeffs
coefs=function(mat){
  apply(mat[,c("alpha","beta")],2,
  function(x) paste0(round(mean(x),2)," (",round(quantile(x,0.025),2),",",
  round(quantile(x,0.975),2),")"))
  }
tab2=cbind(data.frame(model=c("Non-spatial","DAGAR","ICAR")),data.frame(rbind(coefs(mat_ns),coefs(mat_dagar),coefs(mat_icar))))

tab2 %>% left_join(tab1)

## plot w's
slovenia$w.dagar=as.vector(apply(mat_dagar[,paste0("w[",1:N,"]")],2,mean))
slovenia$w.icar=as.vector(apply(mat_icar[,paste0("w[",1:N,"]")],2,mean))

dev.new()
ggplot(data = slovenia) +
    geom_sf() +
    geom_sf(data = slovenia, aes(fill = w.dagar)) +
    scale_fill_viridis_c(direction=-1) +
    labs(fill="w.dagar")

dev.new()
ggplot(data = slovenia) +
    geom_sf() +
    geom_sf(data = slovenia, aes(fill = w.icar)) +
    scale_fill_viridis_c(direction=-1) +
    labs(fill="w.icar")