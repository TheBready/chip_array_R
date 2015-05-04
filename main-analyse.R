#########################################
# Einlesen und Analyse von .CEL-Dateien #
#                 von                   #
#       Nadine, Felix und Philipp       #
#               Gruppe 2                #
#########################################

##################
# Sub-Funktionen #
##################

######################################
# Installieren der benötigten Pakete #
######################################
installPackages <- function(){
  #package source location
  source("http://bioconductor.org/biocLite.R")
  biocLite("affy")
  biocLite("affyPLM")
  biocLite("hgu133plus2.db") # Chip-Datenbank
  biocLite("affyQCReport")
  biocLite("panp")
  biocLite("scatterplot3d")
  biocLite("AffyRNADegradation")
}

##############
# Write Info #
##############
writeInfo <- function(data){
  dir.create("info", showWarnings =FALSE)
  setwd("info")
  # Start writing to an output file
  sink("info.txt")
  print(data)
  print(phenoData(data))
  print(pData(data))
  # Stop writing to the file
  sink()
  setwd("..")
}

##################
# Detection Call #
##################
detectionCall <- function(data,PNGnames,colors,CELnames){
  print("Detection Call")
  dir.create("detection_call", showWarnings =FALSE)
  setwd("detection_call")
  data.exprsSet <- gcrma(data)
  data.PA <- pa.calls(data.exprsSet)
  data.PAcalls <- data.PA$Pcalls
  data.Pvalues <- data.PA$Pvals
  for(n in 1:length(CELnames)){
    png(filename= gsub('.{3}$', '.png', PNGnames[n]))
    hist(data.Pvalues[,n], breaks = 100,border = F, col=colors[n], main=CELnames[n],ylab="Anzahl",xlab="P-Wert")#, ylim = c(0,150000), xlim = c(3,10))
    dev.off() 
  }
  setwd("..")
  #Erstellen der PMA.txt
  print("Erstelle PMA.txt...")
  dir.create("PMA", showWarnings =FALSE)
  setwd("PMA")
  write.table(data.PAcalls, file = "PMA_Calls.txt", row.names = TRUE, quote = FALSE)
  setwd("..")
}

#######
# RMA #
#######
writeRMA <- function(data,dir,j){
  print("Normalisierung RMA")
  dir.create("RMA", showWarnings =FALSE)
  setwd("RMA")
  data.rma <- rma(data)
  write.exprs(data.rma,file = gsub('.{0}$', '_signals_RMA.txt', dir[j]))
  setwd("..")
  return(data.rma)
}


###########
# MAS 5.0 #
###########
writeMAS5 <- function(data,dir,j=1,scale=500){
  print("Normalisierung MAS 5.0")
  dir.create("MAS5", showWarnings =FALSE)
  setwd("MAS5")
  data.mas5 <- mas5(data, sc = scale)
  write.exprs(data.mas5,file = gsub('.{0}$', '_MAS5_500.txt', dir[j]))
  setwd("..")
  return(data.mas5)
}

############################
# Erstellen der Histogramme #
############################
histogramms <- function(data,PNGnames,CELnames,colors,data.mas5exp,data.rmaexp){
  print("Histogramme")
  dir.create("histograms", showWarnings =FALSE)
  setwd("histograms")
  #Einzel-Histogramme
  for(i in 1:length(CELnames)){
    
    # raw
    dir.create("raw", showWarnings =FALSE)
    setwd("raw") 
    png(filename= gsub('.{3}$', '_raw.png', PNGnames[i]))
    hist(log(intensity(data[, i])), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,150000), xlim = c(3,10))
    dev.off()
    setwd("..")
    
    # pm raw
    dir.create("pm_raw", showWarnings =FALSE)
    setwd("pm_raw") 
    png(filename= gsub('.{3}$', '_pm_raw.png', PNGnames[i]))
    hist(log(pm(data[, i])), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,40000), xlim = c(3,9))
    dev.off()
    setwd("..")
    
    
    # mm_raw
    dir.create("mm_raw", showWarnings =FALSE)
    setwd("mm_raw") 
    png(filename= gsub('.{3}$', '_mm_raw.png', PNGnames[i]))
    hist(log(mm(data[, i])), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,70000), xlim = c(3,8))
    dev.off()
    setwd("..")
    
    # mas5
    dir.create("mas5", showWarnings =FALSE)
    setwd("mas5") 
    png(filename= gsub('.{3}$', '_mas5.png', PNGnames[i]))
    hist(log(data.mas5exp[,i]), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)")
    dev.off()
    setwd("..")
    
    # rma
    dir.create("rma", showWarnings =FALSE)
    setwd("rma") 
    png(filename= gsub('.{3}$', '_rma.png', PNGnames[i]))
    hist(log(data.rmaexp[, i]), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)")
    dev.off()   
    setwd("..")
    
    
  }
  
  # Vergleich-Histogramme
  # raw
  setwd("raw") 
  png(filename= "Histogram-Vergleich-raw.png")
  histogramList <- vector('list', 6)
  histogramList[1] <- hist(log(intensity(data[, 1])), breaks = 100)  
  for(m in 2:length(CELnames)){
    histogramList[[m]] <- hist(log(intensity(data[, m])), breaks = 100)
  }
  plot(histogramList[[m]] ,border = F, col=colors[1], main="Histogramme der DNA-Arrays - raw",ylab="Anzahl",xlab="Intensität(log)",ylim=c(0,140000))   
  for(m in 2:length(CELnames)){
    plot( histogramList[[m]],col=colors[m], add=T,border = F)
  }
  legend('topright',colnames(data),fill = colors, bty = 'n',border = NA)
  dev.off()
  
  setwd("..")
  
  # pm raw
  setwd("pm_raw") 
  png(filename= "Histogram-Vergleich-pm_raw.png")
  histogramList <- vector('list', 6)
  histogramList[1] <- hist(log(pm(data[, 1])), breaks = 100)  
  for(m in 2:length(CELnames)){
    histogramList[[m]] <- hist(log(pm(data[, m])), breaks = 100)
  }
  plot(histogramList[[m]] ,border = F, col=colors[1], main="Histogramme der DNA-Arrays - pm_raw",ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,40000), xlim = c(3,8))   
  for(m in 2:length(CELnames)){
    plot( histogramList[[m]],col=colors[m], add=T,border = F)
  }
  legend('topright',colnames(data),fill = colors, bty = 'n',border = NA)
  dev.off()
  
  setwd("..")
  
  
  # mm_raw
  setwd("mm_raw") 
  png(filename= "Histogram-Vergleich-mm_raw.png")
  histogramList <- vector('list', 6)
  histogramList[1] <- hist(log(mm(data[, 1])), breaks = 100)  
  for(m in 2:length(CELnames)){
    histogramList[[m]] <- hist(log(mm(data[, m])), breaks = 100)
  }
  plot(histogramList[[m]] ,border = F, col=colors[1], main="Histogramme der DNA-Arrays - mm_raw",ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,70000), xlim = c(3,8))   
  for(m in 2:length(CELnames)){
    plot( histogramList[[m]],col=colors[m], add=T,border = F)
  }
  legend('topright',colnames(data),fill = colors, bty = 'n',border = NA)
  dev.off()
  
  setwd("..")
  
  # mas5
  setwd("mas5") 
  png(filename= "Histogram-Vergleich-MAS5.0.png")
  histogramList <- vector('list', 6)
  histogramList[1] <- hist(log(data.mas5exp[,1]), breaks = 100)  
  for(m in 2:length(CELnames)){
    histogramList[[m]] <- hist(log(data.mas5exp[,m]), breaks = 100)
  }
  plot(histogramList[[m]] ,border = F, col=colors[1], main="Histogramme der DNA-Arrays - MAS 5.0",ylab="Anzahl",xlab="Intensität(log)")   
  for(m in 2:length(CELnames)){
    plot( histogramList[[m]],col=colors[m], add=T,border = F)
  }
  legend('topright',colnames(data),fill = colors, bty = 'n',border = NA)
  dev.off()
  
  setwd("..")
  
  # rma
  setwd("rma") 
  png(filename= "Histogram-Vergleich-RMA.png")
  histogramList <- vector('list', 6)
  histogramList[1] <- hist(log(data.rmaexp[, 1]), breaks = 100)  
  for(m in 2:length(CELnames)){
    histogramList[[m]] <- hist(log(data.rmaexp[, m]), breaks = 100)
  }
  plot(histogramList[[m]] ,border = F, col=colors[1], main="Histogramme der DNA-Arrays - RMA",ylab="Anzahl",xlab="Intensität(log)")   
  for(m in 2:length(CELnames)){
    plot( histogramList[[m]],col=colors[m], add=T,border = F)
  }
  legend('topright',colnames(data),fill = colors, bty = 'n',border = NA)
  dev.off()
  
  setwd("..")
  setwd("..")
}


##########################
# Erstellen der BoxPlots #
##########################
chipBoxplot <- function(data,data.mas5exp,data.rmaexp){
  print("Erstelle Boxplot")
  dir.create("boxplot", showWarnings =FALSE)
  setwd("boxplot")
  
  #raw
  dir.create("raw", showWarnings =FALSE)
  setwd("raw") 
  png(filename="boxplot_raw.png")
  boxplot(data, col="red")
  dev.off() 
  setwd("..")
  
  # mas5
  dir.create("mas5", showWarnings =FALSE)
  setwd("mas5") 
  png(filename="boxplot_mas5.png")
  boxplot(data.mas5exp, col="red")
  dev.off() 
  setwd("..")
  
  # rma
  dir.create("rma", showWarnings =FALSE)
  setwd("rma") 
  png(filename="boxplot_rma.png")
  boxplot(data.rmaexp, col="red")
  dev.off() 
  setwd("../..")
}


################################
# Ausgabe der Rohdaten in .txt #
################################
rawdata <- function(data,dir,j){
  print("Export Rohdaten")
  dir.create("exprs", showWarnings =FALSE)
  setwd("exprs")
  data.exp <- probes(data)
  write.table(data.exp, file = gsub('.{0}$', '_signals.txt', dir[j]), row.names=TRUE, quote = FALSE)
  setwd("..")
}

##############################################
# Ausgabe der perfect match Rohdaten in .txt #
##############################################
pmdata <- function(data,dir,j,data.proGen){
  print("Export perfect match Daten")
  dir.create("pm", showWarnings =FALSE)
  setwd("pm")
  data.pm <- pm(data,data.proGen[,1])
  write.table(data.pm, file = gsub('.{0}$', '_signals_PM.txt', dir[j]), row.names=TRUE, quote = FALSE)
  setwd("..")
}

#########################################
# Ausgabe der mismatch Rohdaten in .txt #
#########################################
mmdata <- function(data,dir,j,data.proGen){
  print("Export mismatch Daten")
  dir.create("mm", showWarnings =FALSE)
  setwd("mm")
  data.mm <- mm(data,data.proGen[,1])
  write.table(data.mm, file = gsub('.{0}$', '_signals_MM.txt', dir[j]), row.names=TRUE, quote = FALSE)
  setwd("..")
}


###########
# Density #
###########
chipDensity <- function(data,CELnames,j,dir){
  print("Plot Density")
  dir.create("density_plot", showWarnings =FALSE)
  setwd("density_plot") 
  png(filename= "density_plot.png")
  plotDensity.AffyBatch(data, col = 1:length(CELnames), log = TRUE, which=c("pm","mm","both"),ylab = "density", main = dir[j])
  legend("topright",col=1:length(CELnames),lwd=1,legend=CELnames, bty="n")
  dev.off()
  setwd("..")
}

##############
# Clustering #
##############
chipCluster <- function(data,data.rmaexp,data.mas5exp){
  
  # raw data
  print("Erstelle hieraisches Clustering - raw data")
  dir.create("hiera_clust", showWarnings =FALSE)
  setwd("hiera_clust") 
  data.dist = as.matrix(t(exprs(data)))
  data.dist = dist(data.dist,method="euclidean")
  data.cluster = hclust(data.dist, method="average" )
  png(filename="hc_raw.png")
  plot(data.cluster, main= "hieraisches Clustering der Daten - Rohdaten", xlab="Distanz", ylab="Höhe")
  dev.off()
  
  # RMA data
  print("Erstelle hieraisches Clustering - RMA data")
  data.dist = as.matrix(t(data.rmaexp))
  data.dist = dist(data.dist,method="euclidean")
  data.cluster = hclust(data.dist, method="average" )
  png(filename="hc_rma.png")
  plot(data.cluster, main= "hieraisches Clustering der Daten - RMA-Daten", xlab="Distanz", ylab="Höhe")
  dev.off()
  
  # MaAS data
  print("Erstelle hieraisches Clustering - MAS 5.0 data")
  data.dist = as.matrix(t(data.mas5exp))
  data.dist = dist(data.dist,method="euclidean")
  data.cluster = hclust(data.dist, method="average" )
  png(filename="hc_mas.png")
  plot(data.cluster, main= "hieraisches Clustering der Daten - MAS 5.0-Daten", xlab="Distanz", ylab="Höhe")
  dev.off()
  setwd("..")    
}

####################
# Correlation plot #
####################
correlplot <- function(data,data.mas5){
  print("Correlation plot")
  dir.create("correlation_plot", showWarnings = FALSE)
  setwd("correlation_plot")  
  QCReport(data)
  png(filename = "correlation_plot_notNormalized.png")
  correlationPlot(data)
  dev.off()
  png(filename = "correlation_plot_normalized.png")
  correlationPlot(data.mas5)
  dev.off()
  setwd("..")
}

##########################
## Background Intensity ##
##########################

backgroundPlot <- function(data)
{
  print("Erstelle background intensity plot")
  dir.create("background_Intensity", showWarnings =FALSE)
  setwd("background_Intensity") 
  quality <- qc(data)
  
  chips <- colnames(exprs(data))
  num.chips <- length(chips)
  
  #savepath<-paste(savepath, "Background_Intensity_Plot", "background_intensity_plot.jpg", sep="/")
  
  backgroundMean<-mean(quality@average.background)
  backgroundMin <- min(quality@average.background)
  backgroundMax <- max(quality@average.background)
  
  reference<-cbind(c(backgroundMean-10,backgroundMin,backgroundMax-20),c(backgroundMean+10,backgroundMin+20,backgroundMax))
  reftext<-c("[mean-10 ; mean+10]", "[min ; min+20]", "[max-20 ; max]")
  
  x1 <- (quality@average.background >= (backgroundMean-10) & quality@average.background <= (backgroundMean+10))
  x2 <- (quality@average.background >= backgroundMin & quality@average.background <= (backgroundMin+20))  
  x3 <- (quality@average.background >= (backgroundMax-20) & quality@average.background <= backgroundMax)  
  
  outliers <- c(length(x1[x1==FALSE]),length(x2[x2==FALSE]), length(x3[x3==FALSE]))
  if(outliers[1] == 0)
  {
    reference<- reference[1,]
    reftext<-reftext[1]
  }
  else
  {
    reference<- reference[outliers==min(outliers),] # less outliers
    if(length(reference) > 3) { reference<-reference[1,] }
    reftext<- reftext[outliers==min(outliers)] # less outliers
    if(length(reftext) > 1) { reftext<-reftext[1] }  			
  }
  
  testMinimum <- reference[1]
  testMaximum <- reference[2]
  ylimit = c(0, max(max(quality@maximum.background), testMaximum)+5)
  
  png(filename = "background_Intensity_plot.png",width=2000,height=2500,pointsize=40)
  #jpeg(file = savepath,width=2000,height=2500,pointsize=40)
  par(mfrow=c(1,2),oma=c(13,0.5,3,0.5),cex.axis=0.8)
  
  plot(c(quality@minimum.background, quality@maximum.background), type = 'n', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', cex=10)
  par(new=TRUE)
  plot(quality@maximum.background, type = 'h', col=c(2:num.chips+1), ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', lwd=3, ylim = ylimit)
  par(new=TRUE)
  plot(quality@minimum.background, type = 'h', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', lwd=4, col='white', ylim = ylimit)        
  par(new=TRUE)
  plot(quality@maximum.background, type = 'p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=6, cex=0.5, ylim = ylimit)
  par(new=TRUE)
  plot(quality@minimum.background, type = 'p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=2, cex=0.5, ylim=ylimit)
  par(new=TRUE)
  plot(quality@average.background, type='p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=3, cex=1,  ylim=ylimit)
  abline(h=testMinimum, col="grey", ylim=ylimit, lwd=3)
  abline(h=testMaximum, col="grey", ylim=ylimit, lwd=3)	
  title(main="Plot of background intensity",xlab="",ylab="background intensity")
  axis(2)
  par(cex.axis=0.65)
  
  if((num.chips)<20)
  {
    axis(1,at=1:length(quality@average.background),las=2,labels=(substr(chips,1,nchar(chips)-4)))  
  }
  
  legend("bottomright", c("max bg","average bg", "min bg"), col = c(1, 1, 1), pch = c(6,3,2), cex=0.7)
  x1 <- (quality@average.background >= testMinimum & quality@average.background <= testMaximum)
  
  mtext(paste("Data should be between the two lines representing a spread of 20"), side=4, font=1, cex=0.7)
  par(cex.axis=0.8, cex.lab=0.8) 
  
  boxplot(quality@average.background)
  title(main="Average background intensity")
  
  if(backgroundMax - backgroundMin <= 20)
  {
    title(xlab="Background QC: OK (spread <= 20)")
  } 
  else
  {
    title(xlab="Background QC: not OK (spread > 20)")
  }     
  
  legend("bottomright", c(paste("min = ", round(backgroundMin,2), sep="" ), paste("max = ", round(backgroundMax,2), sep="" ), paste("max-min = ", round(backgroundMax-backgroundMin,2), sep="")), cex=0.7)        
  
  mtext("Background Intensity Plot", side=3, outer=TRUE, cex=1.2, font=2)
  
  dev.off()
  
  setwd("..")
  
} 


###################
# RNA Degradation #
###################
RNADegrad <- function(data, colors){
  print("RNA_Degradation")
  dir.create("RNA_Degradation", showWarnings = FALSE)
  setwd("RNA_Degradation")  
  data.rnadeg <- AffyRNAdeg(data)
  png(filename = "RNA_Degradation.png")
  plotAffyRNAdeg(data.rnadeg, cols = colors)
  legend('topleft',colnames(data),fill = colors, bty = 'n',border = NA)
  dev.off()
  setwd("..")
}


###########################
# one gene over all chips #
###########################
geneOverAll <- function(data,data.rmaexp,data.mas5exp){
  print("Es wird ein Gen über alle Chips geplottet")
  dir.create("one_gene_plot", showWarnings = FALSE)
  setwd("one_gene_plot")  
  data.probeset <- probeset(data)
  png(filename = "one_gene_plot_rma.png")
  plot(data.rmaexp["1553602_at",],type ="l", main = "Ein Gen über alle Chips - RMA",ylab = "Intesität", xlab= "Micro-Chip")
  dev.off()
  png(filename = "one_gene_plot_mas5.png")
  plot(data.mas5exp["1553602_at",],type ="l", main = "Ein Gen über alle Chips - MAS 5.0",ylab = "Intesität", xlab= "Micro-Chip")
  dev.off()
  data.ps <- probeset(data, genenames="1553602_at")
  data.psmean <- colMeans(pm(data.ps[[1]]))
  png(filename = "one_gene_plot_raw_mean.png")
  plot(data.psmean,type ="l", main = "Ein Gen über alle Chips - raw mean",ylab = "Intesität", xlab= "Micro-Chip")
  dev.off()
  setwd("..")
}


###############
# Scatterplot #
###############
chipScatter <- function(data,data.rmaexp,data.mas5exp,CELnames){
  print("Scatterplot")
  dir.create("scatterplot", showWarnings = FALSE)
  setwd("scatterplot") 
  names <- gsub('.{4}$', '', CELnames)
  
  #raw
  dir.create("raw", showWarnings = FALSE)
  setwd("raw") 
  for(i in 1:(length(names))){
    for(j in i:length(names)){
      name <- paste(names[i],names[j],sep = " vs ")
      print(name)
      png(filename = paste(name,"_raw.png",sep = ""),width = 1024, height = 1024, units = "px", pointsize = 18)
      trad.scatter.plot(log(exprs(data)[,i]),log(exprs(data)[,j]), main = name, xlab=names[i],ylab=names[j],fc.line.col="red")
      dev.off()
    }
  }
  setwd("..")
  
  #rma
  dir.create("rma", showWarnings = FALSE)
  setwd("rma") 
  for(i in 1:(length(names))){
    for(j in i:length(names)){
      name <- paste(names[i],names[j],sep = " vs ")
      print(name)
      png(filename = paste(name,"_rma.png",sep = ""),width = 1024, height = 1024, units = "px", pointsize = 18)
      trad.scatter.plot(log(data.rmaexp[,i]),log(data.rmaexp[,j]), main = name, xlab=names[i],ylab=names[j],fc.line.col="red")
      dev.off()
    }
  }
  setwd("..")
  
  #mas5
  dir.create("mas5", showWarnings = FALSE)
  setwd("mas5") 
  for(i in 1:(length(names))){
    for(j in i:length(names)){
      name <- paste(names[i],names[j],sep = " vs ")
      print(name)
      png(filename = paste(name,"_mas5.png",sep = ""),width = 1024, height = 1024, units = "px", pointsize = 18)
      trad.scatter.plot(log(data.mas5exp[,i]),log(data.mas5exp[,j]), main = name, xlab=names[i],ylab=names[j],fc.line.col="red")
      dev.off()
    }
  }
  setwd("..")
  setwd("..")
}

############
# QC-stats #
############
qc_stats_plot<-function(data){
  print("QC Stats")
  dir.create("QC Stats", showWarnings = FALSE)
  setwd("QC Stats") 
  qc_stats<-qc(data)
  png(filename = "QC_Stats-test.png",width = 1920, height = 1080, units = "px", pointsize = 24)
  plot(qc_stats)
  dev.off()
  setwd("..")
}




#######
# PCA #
#######
chipPCA <- function(data,data.rmaexp, data.mas5exp, CELnames){
  print("PCA")
  dir.create("PCA", showWarnings = FALSE)
  setwd("PCA")
  
  #Rohdaten
  #raw
  dir.create("raw", showWarnings = FALSE)
  setwd("raw")
  PCA<-prcomp(t(exprs(data))) #pca Rohdaten
  png(filename = "PCA.png") #Übersicht Hauptkomponenten
  plot(PCA,main = "Hauptkomponentenanalyse - erklärende Varianz" )
  dev.off()
  png(filename = "PCA_data.png")
  plot(PCA$x, col=1, pch=c(1:length(CELnames)), las=1, cex=2, main = "Hauptkomponentenanalyse")
  grid()
  legend("top",legend=CELnames, pch=c(1:length(CELnames)),pt.cex=1.5)
  dev.off()
  setwd("..")
  
  
  #skalierte Daten
  dir.create("scale", showWarnings = FALSE)
  setwd("scale")
  PCA_Scale<-prcomp(t(exprs(data)), scale=TRUE)# pca Rohdaten skaliert
  png(filename = "PCA_Scale.png")
  plot(PCA_Scale,main = "Hauptkomponentenanalyse - erklärende Varianz (Skalierte Daten)")
  dev.off()
  png(filename = "PCA_data_Scale.png")
  plot(PCA_Scale$x, col=1, pch=c(1:length(CELnames)), las=1, cex=2, main = "Hauptkomponentenanalyse skalierte Daten")
  grid()
  legend("top",legend=CELnames, pch=c(1:length(CELnames)),pt.cex=1.5)
  dev.off()
  setwd("..")
  
  #MAS5 Daten
  dir.create("mas5", showWarnings = FALSE)
  setwd("mas5") 
  PCA_MAS5<-prcomp(t(data.mas5exp)) #pca MAS5
  png(filename = "PCA_MAS5.png")
  plot(PCA_MAS5,main = "Hauptkomponentenanalyse - erklärende Varianz (MAS5-Daten)")
  dev.off()
  png(filename = "PCA_data_MAS5.png")
  plot(PCA_MAS5$x, col=1, pch=c(1:length(CELnames)), las=1, cex=2, main = "Hauptkomponentenanalyse MAS5-Daten")
  grid()
  legend("top",legend=CELnames, pch=c(1:length(CELnames)),pt.cex=1.5)
  dev.off()
  setwd("..")
  
  #RMA DAten
  dir.create("rma", showWarnings = FALSE)
  setwd("rma") 
  PCA_RMA<-prcomp(t(data.rmaexp)) #pca RMA
  png(filename = "PCA_RMA.png")  #Übersicht HAuptkomponenten
  plot(PCA_RMA,main = "Hauptkomponentenanalyse - erklärende Varianz (RMA-Daten)" )
  dev.off()
  png(filename = "PCA_data_RMA.png") 
  plot(PCA_RMA$x, col=1, pch=c(1:length(CELnames)), las=1, cex=2, main = "Hauptkomponentenanalyse RMA-Daten")
  grid()
  legend("top",legend=CELnames, pch=c(1:length(CELnames)),pt.cex=1.5)
  dev.off()
  setwd("..")
  
  #Zusammenfassung
  png(filename = "Übersicht erklärende Varianz")
  par(mfrow=c(2,2))
  plot(PCA,main = "Rohdaten-erklärende Varianz")
  plot(PCA_Scale,main = "skalierte Daten-erklärende Varianz")
  plot(PCA_MAS5,main = "MAS5-Daten-erklärende Varianz")
  plot(PCA_RMA,main = "RMA-erklärende Varianz" )
  dev.off()
  
  
  png(filename = "PCA Übersicht")
  par(mfrow=c(2,2))
  plot(PCA$x, col=1, pch=c(1:length(CELnames)), las=1, cex=2, main = "PCA Rohdaten")
  plot(PCA_Scale$x, col=1, pch=c(1:length(CELnames)), las=1, cex=2, main = "PCA skalierte Daten")
  plot(PCA_MAS5$x, col=1, pch=c(1:length(CELnames)), las=1, cex=2, main = "PCA MAS5-Daten")
  plot(PCA_RMA$x, col=1, pch=c(1:length(CELnames)), las=1, cex=2, main = "PCA RMA-Daten")
  dev.off()
  
  setwd("..")
} 





###################################################################################################
#*************************************************************************************************#
###################################################################################################


#################
# Main-Funktion #
#################
mainAnalyse<- function(resolution = 7500,scale = 500){
  ###########################
  # Installieren der Pakete #
  ###########################
  #installPackages()
  
  #####################
  # Laden von Paketen #
  #####################
  library("affy")
  library("hgu133plus2.db")
  library("affyPLM")
  library("affyQCReport")
  library("panp")
  library("gcrma")
  library("simpleaffy")
  library("scatterplot3d") 
  library("AffyRNADegradation")
  
  
  
  
  #########################
  # Einstellen des Pfades #
  #########################
  setwd("input")    # wenn Java genutzt wird
  #setwd("../input") # wenn mit R-Studio gearbeitet wird
  dir <- dir()
  setwd("..")
  
  
  ##################### 
  # globale Variablen #
  #####################                              
  data.proGen <- toTable(hgu133plus2SYMBOL) # erstelle Liste mit allen Probe ids und Gen ids
  
  
  # For-Schleife um alle Experimente im Input Ordner ab zu arbeiten
  for(j in 1:length(dir)){
    
    # Erstellen des Output Ordners  
    print("Erstellen der Output-Odner")
    dir.create("output", showWarnings =FALSE)
    setwd("output")
    dir.create(dir[j], showWarnings =FALSE)
    setwd(dir[j])
    
    
    # Erstellen der LogDatei 
    dir.create("log", showWarnings =FALSE)
    setwd("log")
    sink("log.txt")
    setwd("../../..")
    print(paste("Bearbeite Experiment:", dir[j], sep=" "))
    
    # Laden der .cel Files als Batch (alle in einem Ordner)
    setwd("input")
    print("Laden der Daten")
    data <- ReadAffy(celfile.path=dir[j],verbose = TRUE)
    setwd("../output")
    setwd(dir[j])  
    
    # Namen der verschiedenen Samples
    CELnames <- colnames(data)                  # pure CEL-Datei Namen
    PNGnames <- gsub('.{3}$', 'png', CELnames)  # Namen für Bilder
    colors <- rainbow(length(CELnames), alpha =0.5)         # Farben für Plots
    
    # Aufrufen der Funktionen
    writeInfo(data)
    
    detectionCall(data,PNGnames,colors,CELnames)
    
    data.rma <- writeRMA(data,dir,j)
    data.rmaexp <- exprs(data.rma)
    
    data.mas5 <- writeMAS5(data,dir,j,scale)
    data.mas5exp <- exprs(data.mas5)
    data.mas5calls <- mas5calls(data)
    
    histogramms(data,PNGnames,CELnames,colors,data.mas5exp,data.rmaexp)
        
    chipBoxplot(data,data.mas5exp,data.rmaexp)
    
    rawdata(data,dir,j)
    
    pmdata(data,dir,j,data.proGen)
    
    mmdata(data,dir,j,data.proGen)
    
    chipDensity(data,CELnames,j,dir)
    
    chipCluster(data,data.rmaexp,data.mas5exp)
    
    correlplot(data,data.mas5)
    
    geneOverAll(data,data.rmaexp,data.mas5exp)
    
    RNADegrad(data,colors)
    
    chipScatter(data,data.rmaexp,data.mas5exp,CELnames)
    
    qc_stats_plot(data)
      
    backgroundPlot(data)
    
    chipPCA(data, data.rmaexp, data.mas5exp, CELnames)
    
    
    # Ende eines Experiment -> Verlasse Ordner
    print("Bearbeiten des Experimentes beendet")
    sink()
    setwd("../..")
  }
}

mainAnalyse()
