#1. Preperations

#Load libraries
library(affy)
library(vctrs) #necessary? 
library(vsn)
library(limma)
library(AnnotationDbi)
library(highr)
library(mouse4302mmenstcdf) #download source files from brain array for "Affymetrix Mouse Genome 430 2.0 Array"
library(mouse4302mmenstprobe) #download source files from brain array for "Affymetrix Mouse Genome 430 2.0 Array"
library(hexbin, lib.loc = "C:/Users/lydiasteiner/Desktop/Datascience/project/rawdata")
library(hgu133plus2hsenstprobe, lib.loc = "C:/Users/lydiasteiner/Desktop/Datascience/project/rawdata")
library(hgu133plus2hsenstcdf, lib.loc = "C:/Users/lydiasteiner/Desktop/Datascience/project/rawdata")
library(tidyverse)

#Show me the version numbers
sessionInfo()

#Read in .CEL files
setwd("C:/Users/lydiasteiner/Desktop/Datascience/project/rawdata")
data.mouse=ReadAffy()
data.mouse@cdfName<-"HGU133Plus2_Hs_ENST"
setwd("/Users/lydiasteiner/Desktop/Datascience/project/sessions/rda")
save.image(file="rawdata.mouse.rda")

#2.   Quality Control
#-------------------------
#2.1  Single Chip control
#-------------------------
image(data.mouse, col=rainbow(100, start=0, end=0.75)[100:1])

ind=which(new%in%c())



#     2.2 Read in pheno data
#----------------------------------
#remove file-endings from sample names
filenames = rownames(pData(data.mouse))
samples = substr(filenames, 1,9)
rownames(pData(data.mouse))=filenames

#3    Normalization
#----------------------
breast.vsnrma<-vsnrma(data.mouse)

#4.1  Plot meanSdPlot
#--------------------
meanSdPlot(mouse.vsnrma)
setwd()
dev.copy2eps()


#4.2  Boxplot before and after normalization
#--------------------------------------------
#raw data
x11(w=10, h=7)
par(las=2)
mmi=c(1,0.7,1.0477939,0.5366749)
par(mai=mmi)
boxplot(data.mouse, col=rainbow(20), cex.axis=0.5, main="Gene expression in mouse")
setwd()
dev.copy2eps(file="boxplot_mouse_rawdata.eps")

#normalized data
x11(w=10, h=7)
par(las=2)
mmi=c(1,0.7,1.0477939,0.5366749)
par(mai=mmi)
boxplot(exprs(mouse.vsnrma), col=rainbow(20), cex.axis=0.5, main="Gene expression in mouse, normalized")
setwd()
dev.copy2eps(file="boxplot_mouse_vsnrma_normalized.eps")

#if necessary quantile normalization 



