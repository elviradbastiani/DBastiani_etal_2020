############################################################################
#SCRIPT 01
############################################################################
#Mestranda Elvira D'Bastiani
#Programa de Ecologia e Conservacao - UFPR
#Orientadoras: Sabrina B. L. Araujo 
#Laboratorio de Interacoes Biologicas

#reference
#D'Bastiani, E., Campião, K. M., Boeger, W. A., & Araújo, S. B. (2020). The role 
#of ecological opportunity in shaping host–parasite networks. 
#Parasitology, 147(13), 1452-1460.

#clean 
rm(list=ls())
#directory
dir()

############################################################################
#Supplementary material S2
#Description of the models
############################################################################


############################################################################
#Filter i – NEUTRAL MODEL
############################################################################

#COMMENTS

#installing the package
install.packages("bipartite")

#open the directory
setwd(choose.dir())

#Loading the package
library(bipartite)

#inform how many networks to simulate:
samples =1000

#load interaction database matrix 
AS<-data.frame(read.csv2("matrix.database.csv",header=TRUE))

#load observed environment matrix 
observed <-data.frame(read.csv2("matrix.observed .csv",header=TRUE))

#creating lines 
row.names(AS)=AS[,1]
row.names(observed )=observed [,1]
AS=as.matrix(AS[,2:ncol(AS)])
observed =as.matrix(observed [,2:ncol(observed )])
observed.bin=observed.bin=((observed >0)*1)

#metrics for the observed
#number of hosts 
H_observed =nrow(observed.bin)
#number of parasites 
P_observed =ncol(observed.bin)
#nestedness
NODF_observed =nested(observed.bin,method="NODF2",rescale=F)
NODFarr_observed =round(NODF_observed ,digits=3) 
#connectance
Connectance_observed =sum(observed.bin)/(ncol(observed.bin)*nrow(observed.bin))
#modularity
Modular_observed =computeModules(observed .bin)
Modularity_observed =slot(Modular_observed ,"likelihood")

#number of rows and columns
nl = nrow (observed)
nc = ncol (observed)

#fazer para todas as redes amostradas  
for (n in 1:samples){
  sub=AS
  while(nrow(sub)> nrow(observed )){
    
    #return the hosts
    ls=floor(runif(1, 1,nrow(sub)+1)) 
    sub=sub[-ls,]
    
    #deletes parasites that do not interact
    sub=sub[,!!colSums(sub)] 
    sub=sub[!!rowSums(sub),]
  }
  
  #save network with name of parasite and host species
  name <-paste(n,"subrede_observed _aleatorio.txt", sep=" ")
  write.table(sub,name ,col.names = NA, row.names=T,sep="\t")
  
  #metrics for simulated network 
  #number of hosts 
  Hosts=nrow(sub)
  
  #number of parasites 
  Parasites=ncol(sub)
  
  #(0 e 1)
  sub.bin=((sub>0)*1)
  
  #nestedness
  NODF=nested(sub.bin,method="NODF2",rescale=F)
  NODFarr=round(NODF,digits=3)   
  #connectance
  Connectance=sum(sub.bin)/(ncol(sub.bin)*nrow(sub.bin))
  #modularity
  Modular=computeModules(sub.bin <- sub.bin + 1E-5)
  Modularity=slot(Modular,"likelihood")
  
  #number of rows and columns
  nl = nrow (observed)
  nc = ncol (observed)
  
  if (n == 1) {
    exit = c (n, Hosts, Parasites, NODFarr, Connectance, Modularity)
  } else {
    exit1 = c (n, Hosts, Parasites, NODFarr, Connectance, Modularity)
    exit = rbind (exit, exit1)
  }
}

#creating final random model output file
exit = as.data.frame (exit) 
colnames (exit) [1] <- "n"
colnames (exit) [2] <- "Hosts"
colnames (exit) [3] <- "Parasites"
colnames (exit) [4] <- "Nestedness"
colnames (exit) [5] <- "Connectance"
colnames (exit) [6] <- "Modularity"

write.table (exit, file = " exit.random.model.txt", col.names = TRUE, row.names = FALSE, sep = "\ t")

############################################################################
#End neutral model
############################################################################


############################################################################
#Filter ii – Phylogeny Model
############################################################################

#COMMENTS
#installing the package
install.packages("bipartite")

#open the directory
setwd(choose.dir())

#Loading the package
library(bipartite)

#inform how many networks to simulate:
samples =1000

#load interaction database matrix 
AS<-data.frame(read.csv2("matrix.database.csv",header=TRUE))

#load observed environment matrix 
observed <-data.frame(read.csv2("matrix.observed .csv",header=TRUE))

#load the spreadsheet with the host family analyzed
HOSTS<-data.frame(read.csv2("hosts_america_do_sul.csv",header=TRUE))

#inform how many hosts of each family to contain in the simulated network
nHy=6
nameHy="Hylidae"
nLep=5
nameLep="Leptodactylidae"

#creating lines
row.names(AS)=AS[,1]
row.names(observed)=observed[,1]
AS=as.matrix(AS[,2:ncol(AS)])
observed=as.matrix(observed[,2:ncol(observed)])
observed.bin=observed.bin=((observed>0)*1)

#selecting the name of the parasites
namesParasites=colnames(AS)

#considering the same number of rows and columns
nl=nrow(observed)
nc=ncol(observed)    

#selecting families
HHy=as.data.frame(HOSTS[HOSTS$family==nameHy,])
HLep=as.data.frame(HOSTS[HOSTS$family==nameLep,])  

# catch only hylidae hosts 
HHy=HHy$hosts
# catch only leptodactylidae
HLep=HLep$hosts

#selecting the parasite species that interact with hylidae
RedeHy1=AS[HHy,]
#selecting the parasite species that interact with leptodactylidae
RedeHLep1=AS[HLep, ]

#number of hosts
Hosts=nrow(observed.bin)
#number of parasites
Parasites=ncol(observed.bin)

#metrics for observed environment
#Nestedness
NODF_observed=nested(observed.bin,method="NODF2",rescale=F)
NODFarr_observed=round(NODF_observed,digits=3) 
#Connectance 
Connectance_observed=sum(observed.bin)/(ncol(observed.bin)*nrow(observed.bin))
#Modularity
Modular_observed=computeModules(observed.bin)
Modularity_observed=slot(Modular_observed,"likelihood")

#creating files to control the frequency of networks with> and <number of parasites of this model
#less
hostless=matrix(0,nrow (AS), dimnames = list (row.names(AS),"frequencia menor")) 
#major
hostmajor=matrix(0,nrow (AS), dimnames = list (row.names(AS),"frequencia maior")) 

#do for all sampled networks
for (n in 1:samples){  
  #selecting only the number of hosts in each family
  Sub_Hy=RedeHy1[sample.int(length(HHy),nHy),]  
  Sub_Lep=RedeHLep1[sample.int(length(HLep),nLep),] 
  
  #combining columns and rows
  sub=rbind(Sub_Lep,Sub_Hy)
  
  # Eliminating parasites that do not interact
  sub=sub[,!!colSums(sub)] 
  sub=sub[!!rowSums(sub),]
  
  #number of hosts
  Hosts=nrow(sub)
  #number of parasites 
  Parasites=ncol(sub)
  
  #save higher and lower frequency networks
  #Only host names that appear in the sub
  host.sub=rownames(sub)
  #major
  if (Parasites>72) { 
    name<-paste(n,"subrede_observed_F_maior.txt", sep=" ")
    write.table(sub, name, col.names=NA,row.names=TRUE,sep="\t") 
    #array greater than 72 Parasites adds one when generated      hostmajor[host.sub,]=1+hostmajor[host.sub,]
  } 
  #less
  else { 
    name<-paste(n,"subrede_observed_F_menor.txt", sep=" ")
    write.table(sub, name, col.names=NA,row.names=TRUE,sep="\t") 
    #array less than 72 Parasites adds one when generated      hostless[host.sub,]=1+hostless[host.sub,]
  }
  
  #metrics for simulated network
  #number of Hosts
  Hosts=nrow(sub)
  #number Parasites
  Parasites=ncol(sub)
  #binary (0 e 1)
  sub.bin=((sub>0)*1)
  
  #save network with name of Parasite and Host species
  name<-paste(n,"subrede_observed_F.txt", sep=" ")
  write.table(sub, name, col.names=FALSE,row.names=FALSE,sep="\t")
  
  #Nestedness
  NODF_AS=nested(sub.bin,method="NODF2",rescale=F)
  NODFarr_AS=round(NODF_AS,digits=3) 
  #Connectance 
  Connectance_AS=sum(sub.bin)/(ncol(sub.bin)*nrow(sub.bin))
  #Modularity
  Modular_AS=computeModules(sub.bin <- sub.bin + 1E-5)
  Modularity_AS=slot(Modular_AS,"likelihood")
  
  if (n==1) {
    exit=c(n,Hosts,Parasites,NODFarr_AS,Connectance_AS,Modularity_AS)
  } else {
    exit1=c(n,Hosts,Parasites,NODFarr_AS,Connectance_AS,Modularity_AS)
    exit=rbind(exit,exit1)
  }
}
#creating final output file of phylogeny model
exit=as.data.frame(exit)
colnames(exit)[1] <- "n"
colnames(exit)[2] <- "Hosts"
colnames(exit)[3] <- "Parasites"
colnames(exit)[4] <- "Nestedness"
colnames(exit)[5] <- "Connectance"
colnames(exit)[6] <- "Modularity"

write.table(exit,file = "exit.phylogeny.model.txt", col.names=TRUE,row.names=FALSE,sep="\t")


############################################################################
#End Phylogeny Model
############################################################################


############################################################################
#Filter iii – Body size Model
############################################################################

#installing the package
install.packages("bipartite")

#open the directory
setwd(choose.dir())

#Loading the package
library(bipartite)

#inform how many networks to simulate:
samples =1000

#load interaction database matrix 
AS<-data.frame(read.csv2("matrix.database.csv",header=TRUE))

#load observed environment matrix 
observed <-data.frame(read.csv2("matrix.observed .csv",header=TRUE))

#load interaction database matrix - AS = south america
AS_size <-data.frame(read.csv2("hosts_america_do_sul.csv",header=TRUE,dec="."))
#load observed environment matrix
observed_size<-data.frame(read.csv2("Hosts. environment .observado.csv",
                                    header=TRUE,dec="."))

#creating lines
row.names(AS)=AS[,1]
row.names(observed)=observed[,1]
AS=as.matrix(AS[,2:ncol(AS)])
observed=as.matrix(observed[,2:ncol(observed)])
observed.bin=observed.bin=((observed>0)*1)

#selecting the name of the parasites 
namesParasites=colnames(AS)
#selecting the name of the hosts
namesHospedeiro=rownames(AS)

#considering the same number of rows and columns 
nl=nrow(observed)
nc=ncol(observed)   

#number of hosts 
Hosts=nrow(observed.bin)
#number of parasites 
Parasites=ncol(observed.bin)

#metrics for the observed environment
#nestedness
NODF_observed=nested(observed.bin,method="NODF2",rescale=F)
NODFarr_observed=round(NODF_observed,digits=3) 
#connectance 
Connectance_observed=sum(observed.bin)/(ncol(observed.bin)*nrow(observed.bin))
#modularity
modular_observed=computeModules(observed.bin)
Modularity_observed=slot(modular_observed,"likelihood")

#creating network frequency file 
host=matrix(0,nrow (AS), dimnames = list (row.names(AS),"frequencia")) 

#do for all sampled networks
for (n in 1:samples){
  #list the chosen species
  shortlist =as.character(matrix(0,nl,1))
  #create a host source
  fonte_AS=as.data.frame(AS_size )
  #confirming that AS_source is character AS_source
  fonte_AS$hosts<-as.character(fonte_AS$hosts)
  
  for(j in 1: nl){
    #body size reference
    ref_size=as.data.frame(observed_size [,c(4,5)])
    #minimum size 
    minimumsize=observed_size $minimum[j]
    #maximum size 
    maximumsize=observed_size $maximum[j]
    #get only species within the minimum and maximum range 
    sample=fonte_AS[fonte_AS$size>=minimumsize & fonte_AS$size<=maximumsize,]
    #create a sample with these species
    Host_sample=sample$hosts
    #randomly pick a species
    sppchosen<-as.character(sample(Host_sample,1))
    #create a list of the chosen species
    shortlist[j]=sppchosen
    #removing source pick
    fonte_AS=fonte_AS[fonte_AS$hosts!=sppchosen,]
    #write the source
    write.table(fonte_AS,file = "fonte_AS.txt", col.names=TRUE,	
                row.names=F,sep="\t")
  }
  
  # get only lines from the chosen list in AS
  #creating the sub array
  sub=AS[shortlist,]
  
  #save network with parasite and host species
  name <-folder (n, "subrede_observed_T.txt", sep=" ")
    	write.table(sub,name, col.names = NA, row.names=T,sep="\t")
    
    #eliminate non-interacting parasites
    	sub=sub[,!!colSums(sub)] 
    	sub=sub[!!rowSums(sub),]
    
  #metrics for simulated network
    #number of hosts
    Hosts=nrow(sub)
    #number of Parasites
    	    Parasites=ncol(sub)
    #binary (0 and 1)    
    sub.bin=((sub>0)*1)
    
  #save networks
 #only hostnames that appear in sub
host.sub=rownames(sub)
name<-paste(n,"subrede_tamanho.txt", sep=" ")
 #write the network 
write.table(sub, name, col.names=FALSE,row.names=FALSE,sep="\t") 
  
    #n hosts
    	Hosts_AS=nrow(sub.bin)
    #n parasites
    	Parasites_AS=ncol(sub.bin)
    
    #nestedness
    	NODF_AS=nested(sub.bin,method="NODF2",rescale=F)
   	 NODFarr_AS=round(NODF_AS,digits=3) 
    #connectance 
    	Connectance_AS=sum(sub.bin)/(ncol(sub.bin)*nrow(sub.bin))
    #modularity
    	modular_AS=computeModules(sub.bin <- sub.bin + 1E-5)
    	Modularity_AS=slot(modular_AS,"likelihood")
       
        
    if (n==1) {
      exit=c(n,Hosts_AS,Parasites_AS,NODFarr_AS,Connectance_AS,Modularity_AS)
    } else {
exit1=c(n,Hosts_AS,Parasites_AS,NODFarr_AS,Connectance_AS,Modularity_AS)
      exit=rbind(exit,exit1)
    }
    
} 

#creating final output file of phylogeny model
exit=as.data.frame(exit)
colnames(exit)[1] <- "n"
colnames(exit)[2] <- "Hosts"
colnames(exit)[3] <- "Parasites"
colnames(exit)[4] <- "Nestedness"
colnames(exit)[5] <- "Connectance"
colnames(exit)[6] <- "Modularity"

write.table(exit,file = "exit.body.size.model.txt", col.names=TRUE,row.names=FALSE,sep="\t")
  
############################################################################
#End Body size Model
############################################################################
