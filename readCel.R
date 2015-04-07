######################################
#       Einlesen von  .CEL-Dateien   #
#               von                  #
#     Nadine, Felix und Philipp      #
#             Gruppe 2               #
######################################

##################
# Sub-Funktionen #
##################

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
detectionCall <- function(data,PNGnames,colors){
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
    hist(log(pm(data[, i])), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,150000), xlim = c(3,10))
    dev.off()
    setwd("..")
    
    # mm_raw
    dir.create("mm_raw", showWarnings =FALSE)
    setwd("mm_raw") 
    png(filename= gsub('.{3}$', '_mm_raw.png', PNGnames[i]))
    hist(log(mm(data[, i])), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,150000), xlim = c(3,10))
    dev.off()
    setwd("..")
    
    # mas5
    dir.create("mas5", showWarnings =FALSE)
    setwd("mas5") 
    png(filename= gsub('.{3}$', '_mas5.png', PNGnames[i]))
    hist(log(data.mas5exp[,i]), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,2000), xlim = c(0,25))
    dev.off()
    setwd("..")
    
    # rma
    dir.create("rma", showWarnings =FALSE)
    setwd("rma") 
    png(filename= gsub('.{3}$', '_rma.png', PNGnames[i]))
    hist(log(data.rmaexp[, i]), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,1000), xlim = c(0,5))
    dev.off()
    setwd("..")
  }
  setwd("..")
}

########################
# Erstellen der Bilder #
########################
chipImages <- function(data,PNGnames,resolution){
  print("Erstelle Bilder")
  dir.create("images", showWarnings =FALSE)
  setwd("images")
  data.Pset <- fitPLM(data)
  for(i in 1:length(PNGnames)){
    
    print(PNGnames[i])
    #raw
    dir.create("raw", showWarnings =FALSE)
    setwd("raw") 
    png(filename=PNGnames[i], width = resolution, height = resolution, units = "px")
    image(data[,i])
    dev.off()  
    setwd("..")
    
    #topo
    dir.create("topo", showWarnings =FALSE)
    setwd("topo") 
    png(filename=gsub('.{3}$', '_topo.png', PNGnames[i]), width = resolution, height = resolution, units = "px")
    image(data.Pset,which=i)
    dev.off()
    setwd("..")
    
    #heat
    dir.create("heat", showWarnings =FALSE)
    setwd("heat") 
    png(filename=gsub('.{3}$', '_heat.png', PNGnames[i]), width = resolution, height = resolution, units = "px")
    image(data.Pset,which=i,col=pseudoPalette(low="yellow",high="red"))
    dev.off()
    setwd("..")
    
    #palm
    dir.create("palm", showWarnings =FALSE)
    setwd("palm") 
    png(filename=gsub('.{3}$', '_palm.png', PNGnames[i]), width = resolution, height = resolution, units = "px")
    image(data.Pset,which=i,col=pseudoPalette(low="green",high="blue"))
    dev.off() 
    setwd("..")
    
    #resids
    dir.create("resids", showWarnings =FALSE)
    setwd("resids") 
    png(filename=gsub('.{3}$', '_resids.png', PNGnames[i]), width = resolution, height = resolution, units = "px")
    image(data.Pset,which=i, type="resids")
    dev.off()
    setwd("..")
    
    #pos.resids
    dir.create("posResids", showWarnings =FALSE)
    setwd("posResids") 
    png(filename=gsub('.{3}$', '_pos.resids.png', PNGnames[i]), width = resolution, height = resolution, units = "px")
    image(data.Pset,which=i, type="pos.resids",col=pseudoPalette(low="yellow",high="darkblue"))
    dev.off()
    setwd("..")
    
    #neg.resids
    dir.create("negResids", showWarnings =FALSE)
    setwd("negResids") 
    png(filename=gsub('.{3}$', '_neg.resids.png', PNGnames[i]), width = resolution, height = resolution, units = "px")
    image(data.Pset,which=i, type="neg.resids")
    dev.off()
    setwd("..")
    
    #sign.resids
    dir.create("signResids", showWarnings =FALSE)
    setwd("signResids") 
    png(filename=gsub('.{3}$', '_sign.resids.png', PNGnames[i]), width = resolution, height = resolution, units = "px")
    image(data.Pset,which=i, type="sign.resids")
    dev.off()
    setwd("..")
    
    #sign.resids
    dir.create("weight", showWarnings =FALSE)
    setwd("weight") 
    png(filename=gsub('.{3}$', '_weight.png', PNGnames[i]), width = resolution, height = resolution, units = "px")
    image(data.Pset,which=i, type="weight")
    dev.off()
    setwd("..")
    
  }
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
  write.table(data.exp, file = gsub('.{0}$', '_signals.txt', dir[j]), row.names=TRUE)
  setwd("..")
}

##############################################
# Ausgabe der perfect match Rohdaten in .txt #
##############################################
pmdata <- function(data,dir,j){
  dir.create("pm", showWarnings =FALSE)
  setwd("pm")
  data.pm <- pm(data,data.proGen[,1])
  write.table(data.pm, file = gsub('.{0}$', '_signals_PM.txt', dir[j]), row.names=TRUE)
  setwd("..")
}

#########################################
# Ausgabe der mismatch Rohdaten in .txt #
#########################################
mmdata <- function(data,dir,j){
  dir.create("mm", showWarnings =FALSE)
  setwd("mm")
  data.mm <- mm(data,data.proGen[,1])
  write.table(data.mm, file = gsub('.{0}$', '_signals_MM.txt', dir[j]), row.names=TRUE)
  setwd("..")
}


###########
# Density #
###########
chipDensity <- function(data,CELnames){
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
chipCluster <- function(data.exp,data.rmaexp,data.mas5exp){
  
  # raw data
  print("Erstelle hieraisches Clustering")
  dir.create("hiera_clust", showWarnings =FALSE)
  setwd("hiera_clust") 
  data.dist = as.matrix(t(exprs(data)))
  data.dist = dist(data.dist,method="euclidean")
  data.cluster = hclust(data.dist, method="average" )
  png(filename="hc_raw.png")
  plot(data.cluster, main= "hieraisches Clustering der Daten - Rohdaten", xlab="Distanz", ylab="Höhe")
  dev.off()
  
  # RMA data
  data.dist = as.matrix(t(data.rmaexp))
  data.dist = dist(data.dist,method="euclidean")
  data.cluster = hclust(data.dist, method="average" )
  png(filename="hc_rma.png")
  plot(data.cluster, main= "hieraisches Clustering der Daten - RMA-Daten", xlab="Distanz", ylab="Höhe")
  dev.off()
  
  # MaAS data
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
  png(filename = "one_gene_plot_mas5.png")
  plot(data.psmean,type ="l", main = "Ein Gen über alle Chips - MAS 5.0",ylab = "Intesität", xlab= "Micro-Chip")
  dev.off()
  setwd("..")
}


#######
# PCA #
#######
# chipPCA <- function(data){
#   PCA<-prcomp(exprs(data))
#   summary(PCA)                  # Prints variance summary for all principal components.
#   scatterplot3d(PCA$x[,1:3])
#   scatterplot3d(PCA$x[,2:4])
#   scatterplot3d(PCA$x[,3:5])
#   scatterplot3d(PCA$x[,4:6])
#   }


###################################################################################################
#*************************************************************************************************#
###################################################################################################


#################
# Main-Funktion #
#################
mainAnalyse<- function(resolution = 7500,scale = 500){

  ######################################
  # Installieren der benötigten Pakete #
  ######################################
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("affy")
  #biocLite("affyPLM")
  #biocLite("hgu133plus2.db") # Chip-Datenbank
  #biocLite("affyQCReport")
  #biocLite("panp")
  #biocLite("scatterplot3d")

  ##################
  # Laden von affy #
  ##################
  library("affy")
  library("hgu133plus2.db")
  library("affyPLM")
  library("affyQCReport")
  library("panp")
  library("gcrma")
  library("simpleaffy")
  library("scatterplot3d")      



  #########################
  # Einstellen des Pfades #
  #########################
  setwd("../input")
  dir <- dir()
  setwd("..")


  ##################### 
  # globale Variablen #
  #####################                              
  data.proGen <- toTable(hgu133plus2SYMBOL) # erstelle Liste mit allen Probe ids und Gen ids


  # For-Schleife um alle Experimente im Input Ordner ab zu arbeiten
  for(j in 1:length(dir)){
    print(paste("Bearbeite Experiment:", dir[j], sep=" "))
  
  # Laden der .cel Files als Batch (alle in einem Ordner)
    setwd("input")
    print("Laden der Daten")
    data <- ReadAffy(celfile.path=dir[j],verbose = TRUE)
    setwd("..")

  # Namen der verschiedenen Samples
    CELnames <- colnames(data)                  # pure CEL-Datei Namen
    PNGnames <- gsub('.{3}$', 'png', CELnames)  # Namen für Bilder
    colors <- rainbow(length(CELnames))         # Farben für Plots

  # Erstellen des Output Ordners
    print("Erstellen der Output-Odners")
    dir.create("output", showWarnings =FALSE)
    setwd("output")
    dir.create(dir[j], showWarnings =FALSE)
    setwd(dir[j])

  # Aufrufen der Funktionen
    #writeInfo(data)

    #detectionCall(data,PNGnames,colors)

    #data.rma <- writeRMA(data,dir,j)
    #data.rmaexp <- exprs(data.rma)

    #data.mas5 <- writeMAS5(data,dir,j,scale)
    #data.mas5exp <- exprs(data.mas5)
    #data.mas5calls <- mas5calls(data)

    #histogramms(data,PNGnames,CELnames,colors,data.mas5exp,data.rmaexp)

    chipImages(data,PNGnames,resolution)

    #chipBoxplot(data,data.mas5exp,data.rmaexp)

    #rawdata(data,dir,j)

    #pmdata(data,dir,j)

    #mmdata(data,dir,j)

    #chipDensity(data,CELnames)

    #chipCluster(data.exp,data.rmaexp,data.mas5exp)
  
    #correlplot(data,data.mas5)

    #geneOverAll(data,data.rmaexp)

  

  # Ende eines Experiment -> Verlasse Ordner
    print("Bearbeiten des Experimentes beendet")
    setwd("../..")
  }
}