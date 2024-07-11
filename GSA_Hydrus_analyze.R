rm(list=ls())
##############################################################################################
library(data.table)
library(stringr)
library(tidyr)
library(SAFER)
library(tidyverse)
library(caTools)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(FRACTION)
library(matrixStats)
library(devtools)
library(calibrater)
library(readxl)
library(plotly)
library(matrixStats)
##############################################################################################
#
conc.list<-as.data.frame(list.files(path='G:/SA_barrier/donut/SiCl/',pattern = '*Conc'))
Nodes<-lapply(conc.list, function(x) str_extract(x,'\\d+'))
Nodes<-as.numeric(Nodes[[1]])
#

domain<-fread('G:/SA_barrier/donut/SaLo/DOMAIN.orig',skip = 5)
mesh<-fread('G:/SA_barrier/donut/SaLo/MESHTRIA.TXT',skip=1)
xyNodes<-mesh[Nodes,2:3]

Nodes<-cbind(Nodes,xyNodes)
colnames(Nodes)<-c('Nodes','x','z')
# 
profs<-c(20,40,70,100)
tols<-c(2,2,3,3)
for (i in 1:4) {
  locs<-which(abs(Nodes$x-profs[i])<tols[i])
  Nodes[locs,2]<-profs[i]
}
#
profs<-c(-10,-30,-60,-90,-120,-150)
tols<-c(2,2,3,3,3,3)
for (i in 1:6) {
  locs<-which(abs(Nodes$z-profs[i])<tols[i])
  Nodes[locs,3]<-profs[i]
}
#
# # Joining files from parallel simulations
# dir_paral<-'F:/SA_barrier/Puebla_SaLo/paral/paral'  
# file.remove(list.files(path='F:/SA_barrier/Puebla_SaLo/paral/Nodes/',full.names = TRUE))
# for (i in 1:4) {
#   filenames<-paste(dir_paral,as.character(i),'/ConcNode',Nodes$Nodes,'_',as.character(i),'.out',sep='')
#   for (j in 1:length(filenames)) {
#     datos<-fread(filenames[j])
#     Nnode<-Nodes$Nodes[j]
#     write.table(datos,file=paste('F:/SA_barrier/Puebla_SaLo/paral/Nodes/ConcNode',as.character(Nnode),'.out',sep=''),append = TRUE,quote = FALSE,col.names = FALSE)
#   }
# }
#
# # allocation of nodes
# xdist<-c(0,20,40,70,100) # distance from the tree
# zdist<-c(-10,-30,-60,-90,-120,-150) # soil depth
# #
# Nodes.matrix<-matrix(nrow=6,ncol = 5)
# for (i in 1:5) {
#   Nodes.coli<-Nodes[xyNodes$V2>=0.9*xdist[i]&xyNodes$V2<=1.1*xdist[i]]
#   for (j in 1:6) {
#     Nodes.matrix[j,i]<-as.numeric(Nodes.coli[Nodes.coli$V3<=0.9*zdist[j]&Nodes.coli$V3>=1.1*zdist[j],1])
#   }
# }
# 
# ###################################################
# # selection of the depth of the evaluation window
# depth<--10
# # selection of the width of the evaluation window
# width<-40
# pos1<-c(which(zdist==depth,arr.ind = TRUE),which(xdist==width,arr.ind = TRUE))
# p1<-c(1:3)
# p2<-6
# Nodes.sel<-Nodes.matrix[p2,p1]
# #
# # output variables for analysis
# barr.conc<-data.frame(matrix(nrow=400,ncol=length(Nodes.sel)))
# raw.conc<-data.frame(matrix(nrow=1,ncol=length(Nodes.sel)))
# ratio.conc<-data.frame(matrix(nrow=400,ncol=length(Nodes.sel)))
# for (i in 1:length(Nodes.sel)) {
#   barr.sel<-fread(file=paste('F:/SA_barrier/Puebla_SaLo/paral/Nodes/ConcNode',as.character(Nodes.sel[i]),'.out',sep='')) # data from the textural barrier tests
#   raw.sel<-fread(file=paste('F:/SA_barrier/Puebla_SaLo/raw/ConcNode',as.character(Nodes.sel[i]),'_raw.out',sep='')) # data from regular soil (1 Material)
#   barr.conc[,i]<-rowMeans(barr.sel[,-1])
#   raw.conc[,i]<-rowMeans(raw.sel)
#   ratio.conc[,i]<-barr.conc[,i]/raw.conc[,i]
# }
# 
# Meanratio.conc<-rowMeans(ratio.conc)
######################################################
r<-100
M <- 3 
N <- r*(M+1)
distr_fun <- c("unif","unif","unif")
distr_par <- list(c(0, 50), c(0, 40),c(0,30))
design_type <- "radial" # other option is "trajectory"
samp_strat <- "lhs"
#
file.bar<-'G:/SA_barrier/donut/Barrier.par'
X<-fread(file.bar)
X_labels <- c('Distance', 'Width', 'Depth')
EET.means<-data.frame(matrix(nrow=dim(Nodes)[1],ncol=3))
EET.sd<-data.frame(matrix(nrow=dim(Nodes)[1],ncol=3))
Si<-data.frame(matrix(nrow=dim(Nodes)[1],ncol=3))
barrier.ratios<-data.frame(matrix(nrow=dim(X)[1],ncol=dim(Nodes)[1]))

for (i in 1:dim(Nodes)[1]) {
  barr.sel<-fread(file=paste('G:/SA_barrier/donut/ClLo/ConcNode',as.character(Nodes[i,1]),'.out',sep='')) # data from the textural barrier tests
  raw.sel<-fread(file=paste('G:/SA_barrier/donut/raw/ConcNode',as.character(Nodes[i,1]),'.out',sep='')) # data from regular soil (1 Material)
  barr.conc.mean<-rowMeans(barr.sel[,-1])
  raw.conc.mean<-rowMeans(raw.sel)
  barr.conc.end<-barr.sel[,366]
  raw.conc.end<-as.numeric(raw.sel[,366])
  barr.conc.max<-apply(barr.sel,1,max)
  raw.conc.max<-max(raw.sel)
  Yi<-barr.conc.max/raw.conc.max
#  Yi<-barr.conc.max/raw.conc.mean
#  Yi<-barr.conc.end/raw.conc.end
#  Yi<-Yi$V366
  barrier.ratios[,i]<-Yi
  EETind <- EET_indices(r, distr_par, X, Yi, design_type)
  EET.means[i,]<-colMeans(EETind$EE)
  EET.sd[i,]<-colSds(EETind$EE)
  Si[i,]<-as.matrix(EETind[[1]])/sapply(EETind[1],max)
}


colnames(Si)<-c('Si_dist','Si_width','Si_depth')
Si<-data.frame(Nodes,Si)
EET.means<-data.frame(seq(1:dim(Nodes)[1]),Nodes,EET.means)
names(EET.means)[c(1,3:7)]<-c('n','X','Z','dist','width','depth')

locs<-which(Nodes$x<50&Nodes$z>-100)
barrier.ratios.root<-barrier.ratios[,locs]
mean.ratios.root<-rowMeans(barrier.ratios.root)
hist(mean.ratios.root)

#######
bestConfig<-which.max(mean.ratios.root)
X[which.max(mean.ratios.root),]

barr.EC.sel<-data.frame(matrix(nrow=dim(Nodes.sel)[1],ncol=366))
for (i in 1:dim(Nodes.sel)[1]) {
  barr.sel<-fread(file=paste('G:/SA_barrier/donut/ClLo/ConcNode',as.character(Nodes[i,1]),'.out',sep='')) # data from the textural barrier tests
  raw.sel<-fread(file=paste('G:/SA_barrier/donut/raw/ConcNode',as.character(Nodes[i,1]),'.out',sep='')) # data from regular soil (1 Material)
  barr.EC.sel[i,]<-barr.sel[bestConfig,]
  
}

mean.barr.EC.sel<-cbind(mean.barr.EC.sel,colMeans(barr.EC.sel))

colnames(mean.barr.EC.sel)<-c('t','EC_SiCl','EC_SaLo','EC_SaCl','EC_ClLo')
plot(seq(1,366),mean.barr.EC.sel)
#######


Nodes.sel<-Nodes[locs,]
Si.sel<-Si[Nodes$x<50&Nodes$z>-100,]
EET.sd.sel<-EET.sd[Nodes$x<50&Nodes$z>-100,]
EET.sd.sel<-data.frame(Nodes.sel,EET.sd.sel)
orden<-order(-EET.sd.sel$z,EET.sd.sel$x)
EET.sd.sel2<-EET.sd.sel[orden,]

step1<-EET.sd.sel2[, 1:4]
step2<-EET.sd.sel2[, c(1:3,5)]
step3<-EET.sd.sel2[, c(1:3,6)]
colnames(step1)[4]<-'sd.EET'
colnames(step2)[4]<-'sd.EET'
colnames(step3)[4]<-'sd.EET'
step4<-as.factor(rep(c('D','W','H'),each=12))
EET.sd.sel3<-rbind(step1,step2,step3)
EET.sd.sel3<-data.frame(EET.sd.sel3,step4)

ggplot(EET.sd.sel3,aes(x=Nodes,y=sd.EET,fill=step4),legend=list(x=0.8,y=0.2,font=list(size=18)))+
  geom_bar(stat='identity')
 

barrier.ratios.sel<-barrier.ratios[,Nodes$x<50&Nodes$z>-100]
      sapply(barrier.ratios.sel, min)
summary(barrier.ratios.sel)



ecdf_EET.sd<-data.frame(var1=EET.sd[,1],var2=EET.sd[,2],var3=EET.sd[,3],ecdf_1(EET.sd[,1]),ecdf_2(EET.sd[,2]),ecdf_3(EET.sd[,3]))
colnames(ecdf_EET.sd)<-c('dist','width','depth','P_dist','P_width','P_depth')
ggplot(ecdf_EET.sd, aes(x = dist, y = P_dist, color = "red")) +
  geom_line() +
  geom_line(data = ecdf_EET.sd, aes(x = width, y = P_width, color = "blue")) +
  geom_line(data = ecdf_EET.sd, aes(x = depth, y = P_depth, color = "green")) +
  labs(x = "desv.est EET", y = "Probabilidad acumulada") +
  theme_minimal()

best.barr<-matrix(nrow = 12,ncol=3)  
i<-1  
for (i in 1:12) {
  ss<-X[which.min(barrier.ratios.sel[,i]),]
  best.barr[i,]<-as.matrix(ss)
}

X[which.min(barrier.ratios.sel[,i]),]
class(X)
a<-as.matrix(X[which.min(barrier.ratios.sel[,i]),])

ggplot(EET.means,aes(x=X,y=Z,color=dist))+
  geom_point(size=3)+
  scale_colour_gradientn(colours = terrain.colors(10))

ggplot(EET.means,aes(x=X,y=Z,color=width))+
  geom_point(size=3)+
  scale_colour_gradientn(colours = terrain.colors(10))

ggplot(EET.means,aes(x=X,y=Z,color=depth))+
  geom_point(size=3)+
  scale_colour_gradientn(colours = terrain.colors(10))

resultado<-data.frame(Nodes[,2:3],apply(barrier.ratios,2,min))
names(resultado)<-c('X','Z','Min_ratio')

ggplot(resultado,aes(x=X,y=Z,color=Min_ratio))+
  geom_point(size=3)+
  scale_colour_gradientn(colors=terrain.colors(10))

Nodes2<-data.frame(Nodes,apply(barrier.ratios,2,median))
names(Nodes2)<-c('Node','X','Y','Median_ratio')

ggplot(Nodes2,aes(x=X,y=Y,color=Median_ratio))+
  geom_point(size=3)+
  scale_colour_gradientn(colours = terrain.colors(10))

Nodes<-cbind(seq(1,30),Nodes)
names(Nodes)<-c('n','Node','X','Z')
ggplot(data=Nodes,aes(x=X,y=Z))+
  geom_point(size=3)+
  geom_text(aes(label=n),vjust=-0.5)+
  geom_text(aes(label=Node),vjust=1.5)

################################ EET sensitivity analysis ######################
############ Generate inputs: Latin Hypercube sampling - SAFE #########
# 3 inputs, Uniform (continuous) distributions

Yi<-Meanratio.conc
#  Compute indices:



EE <- EETind$EE

mi <- EETind$mi

sigma <- EETind$sigma 

dev.new()

EET_plot(mi, sigma, labels = X_labels)

Nboot <-100

EETind100 <- EET_indices(r, distr_par, X, Yi, design_type, Nboot)

EE <- EETind100$EE

mi <- EETind100$mi

sigma <- EETind100$sigma

mi_lb <- EETind100$mi_lb

mi_ub <- EETind100$mi_ub

sigma_lb <- EETind100$sigma_lb

sigma_ub <- EETind100$sigma_ub

# Plot bootstrapping results in the plane (mean(EE),std(EE)):

dev.new()


EET_plot(mi, sigma, mi_lb, mi_ub, sigma_lb, sigma_ub, labels = X_labels)

rr <- seq(r / 5, r, by = r / 5)

EETconv<- EET_convergence(EE, rr)

m_r <- EETconv$m_r

# Plot the sensitivity measure (mean of elementary effects) as a function 
# of model execution:

dev.new()

plot_convergence(rr * (M + 1), m_r, xlab = "no of model executions", ylab = "mean of EEs", labels = X_labels)
