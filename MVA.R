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
    sink("MVA-log.txt")
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

    
    
    # Ende eines Experiment -> Verlasse Ordner
    print("Bearbeiten des Experimentes beendet")
    sink()
    setwd("../..")
  }
}

mainAnalyse()
