#########################################
# Einlesen und Analyse von .CEL-Dateien #
#             Rohdaten-Verarbeitung     #
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
  write.table(data.PAcalls, file = "PMA_Calls.txt", row.names = TRUE, quote = FALSE, sep = " ")
  setwd("..")
}


################################
# Ausgabe der Rohdaten in .txt #
################################
rawdata <- function(data,dir,j){
  print("Export Rohdaten")
  dir.create("exprs", showWarnings =FALSE)
  setwd("exprs")
  data.exp <- probes(data)
  write.table(data.exp, file = gsub('.{0}$', '_signals.txt', dir[j]), row.names=TRUE, quote = FALSE, sep = " ")
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
  write.table(data.pm, file = gsub('.{0}$', '_signals_PM.txt', dir[j]), row.names=TRUE, quote = FALSE, sep = " ")
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
  write.table(data.mm, file = gsub('.{0}$', '_signals_MM.txt', dir[j]), row.names=TRUE, quote = FALSE, sep = " ")
  setwd("..")
}

############################
# Ausgabe der Gene Symbols #
############################
writeSymbols <- function(data.proGen,dir,j){
  print("Export Gene Symbols")
  dir.create("symbols", showWarnings =FALSE)
  setwd("symbols")
  write.table(data.proGen, file = gsub('.{0}$', '_gene_symbols.txt', dir[j]), quote = FALSE, sep = " ")
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
    sink("raw_analysis_log.txt")
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
    
    writeSymbols(data.proGen,dir,j)
    
    detectionCall(data,PNGnames,colors,CELnames)
                
    rawdata(data,dir,j)
    
    pmdata(data,dir,j,data.proGen)
    
    mmdata(data,dir,j,data.proGen)
    
    chipDensity(data,CELnames,j,dir)
    
    RNADegrad(data,colors)
    
    qc_stats_plot(data)
      
    backgroundPlot(data)
     
    # Ende eines Experiment -> Verlasse Ordner
    print("Bearbeiten des Experimentes beendet")
    sink()
    print("Bearbeiten des Experimentes beendet")
    setwd("../..")
  }
}

mainAnalyse()
