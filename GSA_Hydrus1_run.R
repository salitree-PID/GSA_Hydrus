rm(list=ls())
#
################## Installation of packages ##################
# install.packages('caTools')
# install.packages('FRACTION')
# install.packages('matrixStats')
# install_github("SAFEtoolbox/SAFE-R")
# install.packages('rlang')
# install.packages('data.table')
# install.packages('tidyverse')
# install.packages('cowplot')
# install.packages('gridExtra')
# install.packages('devtools')
# library(devtools)
# install_github("SAFEtoolbox/SAFE-R")
# install.packages('plotly')
#
################## Loading packages ###########################
library(tidyverse)
library(caTools)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(FRACTION)
library(matrixStats)
library(devtools)
library(calibrater)
library(SAFER)
library(readxl)
library(data.table)
library(plotly)
#
# ############ Generate inputs: Latin Hypercube sampling - SAFE #########
# # 3 inputs, Uniform (continuous) distributions
# r<-100
# M <- 3 
# N <- r*(M+1)
# distr_fun <- c("unif","unif","unif")
# distr_par <- list(c(0, 50), c(0, 40),c(0,35))
# samp_strat <- "lhs"
# X <- AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
# #
# # redefinition of minimum and maximum values for each parameter
# X[,1]<-X[,1]+30
# X[,2]<-X[,2]+10
# X[,3]<-X[,3]+5
# barrera<-as.data.frame(X)
# names(barrera)<-list('distance','width','depth')
# write.table(x=X, file = 'E:/donut/Barrier.par',
#             quote = FALSE,row.names = FALSE, 
#             col.names = list('distance','width','depth'),append = FALSE,sep = ',')
# # Plot of the sampling space
# # plot_ly(barrera,x=~distance,y=~width,z=~depth)
# #
barrera<-fread(file='E:/donut/Barrier.par')
##### Archivos comunes Mesh / Domain / Output
direct<-'E:/donut/ClLo/'
setwd(direct)
# mesh
filename1<-'MESHTRIA.TXT'
linea1 <- read_lines(file=paste(direct,filename1,sep=''),skip=0,n_max = 1)
idummy<-as.numeric(substr(linea1,1,10))
NumNP<-as.numeric(substr(linea1,11,20))
nEdges<-as.numeric(substr(linea1,21,30))
NumEl<-as.numeric(substr(linea1,31,40))
#
head_mesh=list('n','x','z')
mesh <- read.fwf(file=paste(direct,filename1,sep=''),
                 skip=1,
                 widths=c(7, 14, 14))

locs<-c(1:dim(mesh)[1])
nn<-as.numeric(mesh$V1)
locs2<-locs[is.na(nn)]
# Extraction of the first part of the mesh file containing coordinates of Nodes
mesh2<-mesh[1:(locs2[1]-1),]
Nodes<-as.numeric(mesh2$V1)
xcoord<-as.numeric(mesh2$V2)
zcoord<-as.numeric(mesh2$V3)
#
# Domain
filename2<-'DOMAIN.orig'
# loads the domain.dat template file (domain.orig)
domain_orig<-fread(paste(direct,filename2,sep=''),skip = 5)
#
# fraccionamiento del fichero domain en 
# encabezado del fichero - encabezado_domain
# sector 1: Nodal Information - domain_orig
# sector 2: BLOCK I: ELEMENT INFORMATION - domain_orig2
domain_orig2<-read.csv(paste(direct,filename2,sep=''),skip=5+NumNP)
encabezado_domain <- read_lines(file=paste(direct,filename2,sep=''),skip=0,n_max = 6)
# Nodes
ObsNodesListF<-'BOUNDARY.IN'
ObsNodesList<-read_lines(file=ObsNodesListF,skip=55,n_max = 1)
ObsNodesList<-as.numeric(t(fread(file=ObsNodesListF,skip=54,nrows = 2)))
NobsNodes<-length(ObsNodesList)
#
columnash<-seq(2,NobsNodes*5+1,5)
columnasTheta<-seq(3,NobsNodes*5+1,5)
columnasConc<-seq(5,NobsNodes*5+1,5)
#
xNodes<-xcoord[ObsNodesList]
zNodes<-zcoord[ObsNodesList]
#
# Nodes at different depths
Nodes10<-ObsNodesList[zNodes>-10]
Nodes50<-ObsNodesList[zNodes>(-53)&zNodes<(-48)]
Nodes100<-ObsNodesList[zNodes>(-104)&zNodes<(-98)]
Nodes200<-ObsNodesList[zNodes>(-205)&zNodes<(-195)]
#
# Root zone data at different depths
Nodes10RZ<-ObsNodesList[zNodes>-10&xNodes<50]
Nodes50RZ<-ObsNodesList[zNodes>(-53)&zNodes<(-48)&xNodes<50]
Nodes100RZ<-ObsNodesList[zNodes>(-104)&zNodes<(-98)&xNodes<50]
Nodes200RZ<-ObsNodesList[zNodes>(-205)&zNodes<(-195)&xNodes<50]
#
# storing vars
Conc10av<-data.frame(matrix(data=NA,nrow=366,ncol=N))
Conc50av<-data.frame(matrix(data=NA,nrow=366,ncol=N))
Conc100av<-data.frame(matrix(data=NA,nrow=366,ncol=N))
Conc200av<-data.frame(matrix(data=NA,nrow=366,ncol=N))
Conc10rz<-data.frame(matrix(data=NA,nrow=366,ncol=N))
Conc50rz<-data.frame(matrix(data=NA,nrow=366,ncol=N))
Conc100rz<-data.frame(matrix(data=NA,nrow=366,ncol=N))
Conc200rz<-data.frame(matrix(data=NA,nrow=366,ncol=N))
# 
hnodei<-data.frame()
thetanodei<-data.frame()
Concnodei<-data.frame()
################## Beginning of the loop ####################
# Modification of materials location in the Domain.dat file #
#
ns<-dim(barrera)[1]
#
for (s in 1:ns) {
  #
  # location of Nodes that will be assigned to the barrier
  Mat2_nodes<-Nodes[xcoord>barrera$distance[s]&xcoord<(barrera$distance[s]+barrera$width[s])&zcoord>-barrera$depth[s]]
  # Cleans the materials and sets all nodes to Mat=1 (original soil)
  domain_orig$M<-1
  # edits material in the limits of the barrier
  domain_orig$M[Mat2_nodes]<-2 
  #
  # # Plot of materials in the soil profile
  # plot(xcoord,zcoord)
  # points(xcoord[Mat2_nodes],zcoord[Mat2_nodes],col='red')
  # points(xNodes,zNodes,col='blue',cex=3,pch=16)
  #
  # Writes the modified DOMAIN file to run Hydrus-2D
  fileOUT<-'E:/donut/ClLo/DOMAIN.DAT'
  file.remove(fileOUT)
  write.table(x=encabezado_domain, file = fileOUT,sep='',
              quote = FALSE,row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(x=domain_orig, file = fileOUT,
              quote = FALSE,row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(x=domain_orig2, file = fileOUT,
              quote = FALSE,row.names = FALSE, col.names = FALSE,append = TRUE)
  #
  #################### Hydrus 2D model execution ####################
  #
  system2('E:/donut/ClLo/H2D_Calc.exe')
  #
  #################### Reading and processing of output data ####################
  #
  ObsNodFile<-'ObsNod.out' 
  ObsNod<-fread(ObsNodFile,skip=5) # reads Observation nodes data. Nodes in the file are sorted
  #
  hNodes<-data.frame(ObsNod[,1],subset(ObsNod,select=columnash))
  thetaNodes<-data.frame(ObsNod[,1],subset(ObsNod,select=columnasTheta))
  ConcNodes<-data.frame(ObsNod[,1],subset(ObsNod,select=columnasConc))
  #
  for (i in 1:NobsNodes) {
    colnames(hNodes)[i+1]<-paste('h','No',ObsNodesList[i],sep='')
    colnames(thetaNodes)[i+1]<-paste('theta','No',ObsNodesList[i],sep='')
    colnames(ConcNodes)[i+1]<-paste('Conc','No',ObsNodesList[i],sep='')
    
    nombre_h<-paste("hNode", ObsNodesList[i], '.out', sep = "")
    write.table(x=t(hNodes[,i+1]),file=nombre_h,row.names=FALSE,col.names=FALSE,append=TRUE)
    nombre_theta<-paste("ThetaNode", ObsNodesList[i],'.out', sep = "")
    write.table(x=t(thetaNodes[,i+1]),file=nombre_theta,row.names=FALSE,col.names=FALSE,append=TRUE)
    nombre_Conc<-paste("ConcNode", ObsNodesList[i],'.out', sep = "")
    write.table(x=t(ConcNodes[,i+1]),file=nombre_Conc,row.names=FALSE,col.names=FALSE,append=TRUE)
  }
  concBottom<-fread('solute1.out',skip=5)
  write.table(x=t(concBottom$V29),file='concBottom.out',row.names=FALSE,col.names=FALSE,append=TRUE)
  #
}
# #####
# # Plot of EC time series at different depths
# fig1<-plot_ly(Conc10,x=Conc10[,1],y=rowMeans(Conc10[,-1]), name='z= -10 cm',type='scatter',mode='lines')
# fig1<- fig1 %>% add_trace(Conc50,y=rowMeans(Conc50[,-1]), name='z= -50 cm', mode='lines')
# fig1<- fig1 %>% add_trace(Conc100,y=rowMeans(Conc100[,-1]), name='z= -100 cm', mode='lines')
# fig1<- fig1 %>% add_trace(Conc200,y=rowMeans(Conc200[,-1]), name='z= -200 cm', mode='lines')
# fig1<- fig1 %>% layout(xaxis=list(title = 'days'),yaxis=list(title = 'EC, dS/m'))
# fig1
# #
# # Plot of EC time series at different depths
# fig2<-plot_ly(Conc10,x=Conc10RZ[,1],y=rowMeans(Conc10RZ[,-1]), name='z= -10 cm',type='scatter',mode='lines')
# fig2<- fig2 %>% add_trace(Conc50RZ,y=rowMeans(Conc50RZ[,-1]), name='z= -50 cm', mode='lines')
# fig2<- fig2 %>% add_trace(Conc100RZ,y=rowMeans(Conc100RZ[,-1]), name='z= -100 cm', mode='lines')
# fig2<- fig2 %>% add_trace(Conc200RZ,y=rowMeans(Conc200RZ[,-1]), name='z= -200 cm', mode='lines')
# fig2<- fig2 %>% layout(xaxis=list(title = 'days'),yaxis=list(title = 'EC, dS/m'))
# fig2
# #
# 
# NodeLocs[ObsNodesList==Nodes200]
# ObsNodesList

####################
#
