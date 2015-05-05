#########################################
# Einlesen und Analyse von .CEL-Dateien #
#           Normalisierungen            #
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

#######
# RMA #
#######
writeRMA <- function(data,dir,j){
  print("Normalisierung RMA")
  dir.create("RMA", showWarnings =FALSE)
  setwd("RMA")
  data.rma <- rma(data)
  write.exprs(data.rma,file = gsub('.{0}$', '_signals_RMA.txt', dir[j]), sep = " ")
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
  write.exprs(data.mas5,file = gsub('.{0}$', '_MAS5_500.txt', dir[j]), sep = " ")
  setwd("..")
  return(data.mas5)
}


#############
# MVA-Plots #
#############
mvaPlot<-function(data,data.rmaexp,data.mas5exp){
  print("MVA")
  dir.create("MVA", showWarnings = FALSE)
  setwd("MVA") 
  
  #raw
  dir.create("raw", showWarnings = FALSE)
  setwd("raw") 
  png(filename = "mva-plot_raw.png",width = 1920, height = 1080, units = "px", pointsize = 24)
  mva.pairs(exprs(data), log.it = FALSE)
  dev.off()
  png(filename = "mva-plot_log_raw.png",width = 1920, height = 1080, units = "px", pointsize = 24)
  mva.pairs(exprs(data), log.it = TRUE)
  dev.off()
  setwd("..")
  
  #rma
  dir.create("rma", showWarnings = FALSE)
  setwd("rma") 
  png(filename = "mva-plot_rma.png",width = 1920, height = 1080, units = "px", pointsize = 24)
  mva.pairs(data.rmaexp, log.it = FALSE)
  dev.off()
  png(filename = "mva-plot_log_rma.png",width = 1920, height = 1080, units = "px", pointsize = 24)
  mva.pairs(data.rmaexp, log.it = TRUE)
  dev.off()
  setwd("..")
  
  #mas5
  dir.create("mas5", showWarnings = FALSE)
  setwd("mas5") 
  png(filename = "mva-plot_mas5.png",width = 1920, height = 1080, units = "px", pointsize = 24)
  mva.pairs(data.rmaexp, log.it = FALSE)
  dev.off()
  png(filename = "mva-plot_log_mas5.png",width = 1920, height = 1080, units = "px", pointsize = 24)
  mva.pairs(data.mas5exp, log.it = TRUE)
  dev.off()
  setwd("..")
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
    sink("Normalisation-log.txt")
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

    
    data.rma <- writeRMA(data,dir,j)
    data.rmaexp <- exprs(data.rma)
    
    data.mas5 <- writeMAS5(data,dir,j,scale)
    data.mas5exp <- exprs(data.mas5)
    
    mvaPlot(data,data.rmaexp,data.mas5exp)
    
    chipPCA(data, data.rmaexp, data.mas5exp, CELnames)
    
    chipScatter(data,data.rmaexp,data.mas5exp,CELnames)
    
    geneOverAll(data,data.rmaexp,data.mas5exp)
    
    chipCluster(data,data.rmaexp,data.mas5exp)
    
    correlplot(data,data.mas5)
    
    chipBoxplot(data,data.mas5exp,data.rmaexp)
    
    histogramms(data,PNGnames,CELnames,colors,data.mas5exp,data.rmaexp)
    
    
    # Ende eines Experiment -> Verlasse Ordner
    print("Bearbeiten des Experimentes beendet")
    sink()
    setwd("../..")
  }
}

mainAnalyse()
