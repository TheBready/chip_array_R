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
  library("affyPLM")
  
  
  #########################
  # Einstellen des Pfades #
  #########################
  setwd("input")    # wenn Java genutzt wird
  #setwd("../input") # wenn mir R-Studio gearbeitet wird
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
    sink("images-log.txt")
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
    
    # Aufrufen der Funktionen
    
    
    chipImages(data,PNGnames,resolution)
    
    
    # Ende eines Experiment -> Verlasse Ordner
    print("Bearbeiten des Experimentes beendet")
    sink()
    setwd("../..")
  }
}

mainAnalyse()
