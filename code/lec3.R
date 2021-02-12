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

### Baltimore boundaries and tools to plot maps
left=-76.76
bottom=39.18
right=-76.42
top=39.43

### Baltimore city map
baltmap = ggmap(get_stamenmap(bbox=c(left=left,bottom=bottom,right=right,top=top),zoom=12))

tab=read.csv("../data/baltpm.csv")
tab


bmap=baltmap+geom_point(data=tab,
    aes(x=Lon,y=Lat,fill=PM25),shape=21,alpha=0.8,size=5,col="black",stroke=2) + 
    scale_fill_viridis_c(direction=-1)
dev.new()
plot(bmap)

### generating UTM projections
### code courtesy: https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
library(sp)
library(rgdal)

#Function
LongLatToUTM<-function(x,y,zone){
 xy <- data.frame(Easting = x, Northing = y)
 coordinates(xy) <- c("Easting", "Northing")
 proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## sets the co-ordinate reference system
 res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
 return(as.data.frame(res))
}

# Example
C=LongLatToUTM(tab$Lon,tab$Lat,18)
tab=cbind(tab,C)

dev.new()
ggplot(data=tab,aes(x=Easting,y=Northing,fill=PM25)) + 
         geom_point(shape=21,alpha=0.8,size=5,col="black",stroke=2) + 
    scale_fill_viridis_c(direction=-1)

### comparing different distances ###

latlonggrid=expand.grid(seq((-179),179,length=60),seq(0,89,length=15))

dmat.geo=rdist.earth(latlonggrid,latlonggrid)

### Earth's radius
R=as.numeric(rdist.earth(t(as.matrix(c(-180,0))),t(as.matrix(c(0,0))))/pi) ## the value used by fields package

### Naive Euclidean distance
dmat.naive.euclid=rdist(latlonggrid,latlonggrid)*pi*R/180

plot(dmat.geo[upper.tri(dmat.geo)],
     dmat.naive.euclid[upper.tri(dmat.geo)],
     xlab="Geodesic",ylab="Naive Euclidean")
abline(a=0,b=1,col="red")

## Mercator projection based distance
mercator.grid=t(apply(latlonggrid,1,function(x) R*c(pi*x[1]/180,log(tan(pi/4+pi*x[2]/360)))))
dmat.mercator=rdist(mercator.grid,mercator.grid)

plot(dmat.geo[upper.tri(dmat.geo)],
     dmat.mercator[upper.tri(dmat.geo)],
     xlab="Geodesic",ylab="Mercator-projection based")
abline(a=0,b=1,col="red")

## chordal distance
euclid.grid=t(apply(latlonggrid,1,function(x)
    R*c(cos(pi*x[1]/180)*cos(pi*x[2]/180),sin(pi*x[1]/180)*cos(pi*x[2]/180),sin(pi*x[2]/180))))
dmat.euclid=rdist(euclid.grid,euclid.grid)

plot(dmat.geo[upper.tri(dmat.geo)],
     dmat.euclid[upper.tri(dmat.geo)],
     ylim=range(dmat.geo),
     xlab="Geodesic",ylab="Chordal")
abline(a=0,b=1,col="red")

#### running Bayesian estimation for a subset of dataset 3 in Stan ####
library(rstan)
#options(mc.cores = parallel::detectCores())
data3=read.csv("../data/dataset3.csv")
data3small=data3[1:50,]

## recalling the MLE
mle <- likfit(coords=data3small[,1:2], data=data3small[,4], trend = trend.spatial(~x,data3small),
  ini.cov.pars=c(0.12,0.2),nugget = 0.02,cov.model="exponential",nospatial=TRUE)

mle

## stan data file
data3list=list(n=nrow(data3small),y=data3small$y,
  x=data3small$x,dmat=as.matrix(dist(data3small[,c("sx","sy")])),id=diag(nrow(data3small)))

## rrunning stan
set.seed(1)
fit1 <- stan(
  file = "bayesgp1.stan",  # Stan program
  data = data3list,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  open_progress = T,
  #cores = 4             # number of cores (could use one per chain)
  #refresh = 10             # no shown every 100 iterations
  )

## save(file="../data/data3small_bayesgp1.Rdata",list=c("fit1"))
## load("../data/data3small_bayesgp1.Rdata")

## run time
get_elapsed_time(fit1)

### posterior estimates, traceplots, densities
plot(fit1, plotfun = "trace",pars=c("alpha","beta","sigmasq","tausq","phi"), inc_warmup = TRUE)
plot(fit1, plotfun = "dens",pars=c("alpha","beta","sigmasq","tausq","phi"), inc_warmup = TRUE)
plot(fit1,pars=c("alpha","beta","sigmasq","tausq","phi"),ci_level=0.95)

### comparing with the MLE
#round(summary(fit1)[[1]][c("alpha","beta","sigmasq","tausq","phi"),c("mean","2.5%","97.5%")],4)
#mle

tab=round(cbind(summary(fit1)[[1]][c("alpha","beta","sigmasq","tausq","phi"),c("mean","2.5%","97.5%")],
      c(mle$beta,mle$sigmasq,mle$tausq,1/mle$phi)),4)
colnames(tab)[4]="mle"
