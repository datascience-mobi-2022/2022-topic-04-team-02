#TRAs in human embryonic development
#---------------------------------
library(affy)
#Alina Aksianova
#date: 02.05.2022
#---------------------------------

#libraries
library(vsn)
library(limma)
library(AnnotationDbi)
library(hexbin, lib.loc = "C:/Users/amaks/AppData/Local/R/win-library/4.2")
library(hgu133plus2hsenstprobe, lib.loc = "C:/Users/amaks/AppData/Local/R/win-library/4.2")
library(hgu133plus2hsenstcdf, lib.loc = "C:/Users/amaks/AppData/Local/R/win-library/4.2")
library(tidyverse)


#read in the data

setwd("C:/Users/amaks/Data analysis project/GSE15744/rawdata")

data_raw=ReadAffy()
data_raw@cdfName <- "HGU133Plus2_Hs_ENST"

setwd("C:/Users/amaks/Data analysis project/GSE15744/Images")
save.image(file="rawdata_15744.rda")

#Quality Control
##single chip control

image(data_raw[,1], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,2], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,3], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,4], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,5], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,6], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,7], col=rainbow(100, start=0, end=0.75)[100:1]) 
image(data_raw[,8], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,9], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,10], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,11], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,12], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,13], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,14], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,15], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,16], col=rainbow(100, start=0, end=0.75)[100:1]) 
image(data_raw[,17], col=rainbow(100, start=0, end=0.75)[100:1])
image(data_raw[,18], col=rainbow(100, start=0, end=0.75)[100:1])

rownames(data_raw)
colnames(data_raw)

#Normalization
#-------------
data_vsnrma <- vsnrma(data_raw)

setwd("C:/Users/amaks/Data analysis project/GSE15744/data")
save.image(file="normalized_data.rda")

view(data_vsnrma)

#meanSdPlot
#----------

meanSdPlot(data_vsnrma)
setwd("C:/Users/amaks/Data analysis project/GSE15744/plots")
dev.copy2eps(file="meanSdPlot_vsnrma_normalized_data_15744.eps")
dev.copy2pdf(file="meanSdPlot_vsnrma_normalized_data_15744.pdf")


#boxplot
#-------

#before normalization
x11(w=10,h=7) #width and height of plotting window

par(las=2) #style of axis labels

mmi=c(1,0.7,1.0477939,0.5366749)
par(mai=mmi) #margin size
boxplot(data_raw, col=rainbow(20), cex.axis=0.5, main="Gene expression in human embryogenesis before normalization")
setwd("C:/Users/amaks/Data analysis project/GSE15744/plots")
dev.copy2eps(file="boxplot_data_raw.eps")
dev.copy2pdf(file="boxplot_data_raw.pdf")

#after normalization

boxplot(exprs(data_vsnrma), col=rainbow(20), cex.axis=0.5, main="Gene expression in human embryogenesis after normalization")
setwd("C:/Users/amaks/Data analysis project/GSE15744/plots")
dev.copy2eps(file="boxplot_data_vsnrma.eps")
dev.copy2pdf(file="boxplot_data_vsnrma.pdf")
dev.off()

#histogram
#----------------

#before normalization

hist(data_raw, col=rainbow(20), main="Density function of log intensity of human embryogenesis before normalization")
dev.copy2pdf(file="Histogram_raw_15774.pdf", width = 10, height = 7)

#after normalization

eset=exprs(data_vsnrma)

plot(density(eset[,1]), type="n", xlab="log Intensity", ylim=c(0,1), main="Density function of log intensity of human embryogenesis after normalization")
for (i in 1:ncol(eset)) {
  lines(density(eset[,i]), col=rainbow(20)[i])
}

dev.copy2pdf(file="Histogram_vsnrma_15774.pdf", width = 10, height = 7)


# RNA degradation plot

rnadeg.raw = AffyRNAdeg(data_raw)

plotAffyRNAdeg(rnadeg.raw, col=rainbow(20))
title(sub="GSE15774 rawdata")

dev.copy2pdf(file="rnadeg_rawdata_15774.pdf")
dev.copy2eps(file="rnadeg_rawdata_15774.eps")


plotAffyRNAdeg(rnadeg.raw, col=rainbow(20), transform = "shift.only")
title(sub="GSE15774 rawdata")

dev.copy2pdf(file="rnadeg_shiftonly_rawdata_15774.pdf")
dev.copy2eps(file="rnadeg_shiftonly_rawdata_15774.eps")

setwd("C:/Users/amaks/Data analysis project/GSE15744/data")
data_vsnrma=load(file="normalized_data.rda")
eset = exprs(data_vsnrma)

#Scatter plots

x11(w=9,h=5)
par(mfrow=c(1,2))

setwd("C:/Users/amaks/Data analysis project/GSE15744/data")
for (i in 1:17){
  
  plot(eset[,c(i,i+1)],pch='.')
  abline(0,1,col="red")
  
  title(main=paste("Scatterplot of probe", substr(colnames(data_vsnrma)[i],1,
                                                  nchar(colnames(data_vsnrma)[i])-4), "and", substr(colnames(data_vsnrma)[i+1],1,
                                                                                                    nchar(colnames(data_vsnrma)[i+1])-4), sep =" ", collapse = NULL))
  
  file.name= paste("scatterplot",as.character(substr(colnames(data_vsnrma)[i],1,
                                                     nchar(colnames(data_vsnrma)[i])-4)), "_and", as.character(substr(colnames(data_vsnrma)[i+1],1,
                                                                                                                      nchar(colnames(data_vsnrma)[i+1])-4)),".pdf", sep="")
  
  dev.copy2pdf(file= file.name, width=9, height=9)
  dev.off()
}

file.name