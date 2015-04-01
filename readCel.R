######################################
#       Einlesen von  .CEL-Dateien   #
#               von                  #
#     Nadine, Felix und Philipp      #
#             Gruppe 2               #
######################################
# Installieren der benötigten Pakete
#source("http://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("affyPLM")
#biocLite("hgu133plus2.db") # Chip-Datenbank
#biocLite("affyQCReport")


# Laden von affy
library("affy")
library("hgu133plus2.db")
library("affyPLM")
library("affyQCReport")

#Einstellen des Pfades
setwd("../input")
dir <- dir()
setwd("..")

# For-Schleife um alle Experimente im Input Ordner ab zu arbeiten
for(j in 1:length(dir)){
  print(paste("Bearbeite Experiment:", dir[j], sep=" "))
  
# Laden der .cel Files als Batch (alle in einem Ordner)
  setwd("input")
  print("Laden der Daten")
  data <- ReadAffy(celfile.path=dir[j],verbose = TRUE)

# Namen der verschiedenen Samples
  CELnames <- colnames(data)

# Erstellen des Output Ordners
  print("Erstellen der Output-Odners")
  setwd("..")
  dir.create("output", showWarnings =FALSE)
  setwd("output")
  dir.create(dir[j], showWarnings =FALSE)
  setwd(dir[j])

# RMA
  print("Normalisierung RMA")
  dir.create("RMA", showWarnings =FALSE)
  setwd("RMA")
  data.rma <- rma(data)
  data.rmaexp <- exprs(data.rma)
  write.exprs(data.rma,file = gsub('.{0}$', '_signals_RMA.txt', dir[j]))

# MAS 5.0
  print("Normalisierung MAS 5.0")
  setwd("..")
  dir.create("MAS5", showWarnings =FALSE)
  setwd("MAS5")
  data.mas5 <- mas5(data, sc = 500)
  data.mas5exp <- exprs(data.mas5)
  data.mas5calls <- mas5calls(data)
  write.exprs(data.mas5,file = gsub('.{0}$', '_MAS5_500.txt', dir[j]))

# Erstellen der Histgramme
  print("Histogramme")
  setwd("..")
  dir.create("histograms", showWarnings =FALSE)
  setwd("histograms")
  PNGnames <- gsub('.{3}$', 'png', CELnames)
  colors <- rainbow(length(CELnames))
  for(i in 1:length(CELnames)){
    # raw
    png(filename= gsub('.{3}$', '_raw.png', PNGnames[i]))
    hist(log(intensity(data[, i])), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,150000), xlim = c(3,10))
    dev.off()
    # mas5
    png(filename= gsub('.{3}$', '_mas5.png', PNGnames[i]))
    hist(log(data.mas5exp[,i]), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,2000), xlim = c(0,25))
    dev.off()
    # rma
    png(filename= gsub('.{3}$', '_rma.png', PNGnames[i]))
    hist(log(data.rmaexp[, i]), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)", ylim = c(0,1000), xlim = c(0,5))
    dev.off()
  }



# Erstellen der Bilder
  print("Erstelle Bilder")
  setwd("..")
  dir.create("images", showWarnings =FALSE)
  setwd("images")
  for(i in 1:length(PNGnames)){
    png(filename=PNGnames[i], width = 7500, height = 7500, units = "px")
    image(data[,i])
    dev.off()  
    data.Pset <- fitPLM(data)
    png(filename=gsub('.{3}$', '_topo.png', PNGnames[i]), width = 7500, height = 7500, units = "px")
    image(data.Pset,which=i)
    dev.off()
    png(filename=gsub('.{3}$', '_heat.png', PNGnames[i]), width = 7500, height = 7500, units = "px")
    image(data.Pset,which=i,col=pseudoPalette(low="yellow",high="red"))
    dev.off()
    png(filename=gsub('.{3}$', '_palm.png', PNGnames[i]), width = 7500, height = 7500, units = "px")
    image(data.Pset,which=i,col=pseudoPalette(low="green",high="blue"))
    dev.off() 
    png(filename=gsub('.{3}$', '_resids.png', PNGnames[i]), width = 7500, height = 7500, units = "px")
    image(data.Pset,which=2, type="resids")
    dev.off()
    png(filename=gsub('.{3}$', '_pos.resids.png', PNGnames[i]), width = 7500, height = 7500, units = "px")
    image(data.Pset,which=2, type="pos.resids",col=pseudoPalette(low="yellow",high="darkblue"))
    dev.off()
    png(filename=gsub('.{3}$', '_neg.resids.png', PNGnames[i]), width = 7500, height = 7500, units = "px")
    image(data.Pset,which=2, type="neg.resids")
    dev.off()
    png(filename=gsub('.{3}$', '_sign.resids.png', PNGnames[i]), width = 7500, height = 7500, units = "px")
    image(data.Pset,which=2, type="sign.resids")
    dev.off()
  }

# Erstellen der BoxPlots
  print("Erstelle Boxplot")
  setwd("..")
  dir.create("boxplot", showWarnings =FALSE)
  setwd("boxplot")
  png(filename="boxplot.png")
  boxplot(data, col="red")
  dev.off()  


# erstelle Liste mit allen Probe ids und Gen ids
  data.proGen <- toTable(hgu133plus2SYMBOL)

#Ausgabe der Rohdaten in .txt
  print("Export Rohdaten")
  setwd("..")
  dir.create("exprs", showWarnings =FALSE)
  setwd("exprs")
  data.exp <- probes(data)
  write.table(data.exp, file = gsub('.{0}$', '_signals.txt', dir[j]), row.names=TRUE)

# Ausgabe der perfect match Rohdaten in .txt
  setwd("..")
  dir.create("pm", showWarnings =FALSE)
  setwd("pm")
  data.pm <- pm(data,data.proGen[,1])
  write.table(data.pm, file = gsub('.{0}$', '_signals_PM.txt', dir[j]), row.names=TRUE)

# Ausgabe der mismatch Rohdaten in .txt
  setwd("..")
  dir.create("mm", showWarnings =FALSE)
  setwd("mm")
  data.mm <- mm(data,data.proGen[,1])
  write.table(data.mm, file = gsub('.{0}$', '_signals_MM.txt', dir[j]), row.names=TRUE)


# RMA
  print("Normalisierung RMA")
  setwd("..")
  dir.create("RMA", showWarnings =FALSE)
  setwd("RMA")
  data.rma <- rma(data)
  write.exprs(data.rma,file = gsub('.{0}$', '_signals_RMA.txt', dir[j]))

# MAS 5.0
  print("Normalisierung MAS 5.0")
  setwd("..")
  dir.create("MAS5", showWarnings =FALSE)
  setwd("MAS5")
  data.mas5 <- mas5(data, sc = 500)
  data.mas5exp <- exprs(data.mas5)
  data.mas5calls <- mas5calls(data)
  write.exprs(data.mas5,file = gsub('.{0}$', '_MAS5_500.txt', dir[j]))

# Density
  print("Plot Density")
  setwd("..")
  dir.create("density_plot", showWarnings =FALSE)
  setwd("density_plot") 
  png(filename= "density_plot.png")
  plotDensity.AffyBatch(data, col = 1:length(CELnames), log = TRUE, which=c("pm","mm","both"),ylab = "density", main = dir[j])
  legend("topright",col=1:length(CELnames),lwd=1,legend=CELnames, bty="n")
  dev.off()

# QC
  data.qc <- qc(data)
  data.back <- avbg(data.qc)
  plot(data.back)

# Clustering 
  # raw data
  print("Erstelle hieraisches Clustering")
  setwd("..")
  dir.create("hiera_clust", showWarnings =FALSE)
  setwd("hiera_clust") 
  data.dist = as.matrix(t(data.exp))
  data.dist = dist(data.dist,method="euclidean")
  data.cluster = hclust(data.dist, method="average" )
  png(filename="hc_raw.png")
  plot(data.cluster, main= "hieraisches Clustering der Daten - Rohdaten", xlab="Distanz", ylab="Höhe")
  dev.off()
  # RMA data
  data.dist = as.matrix(t(exprs(data.rma)))
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

# Correlation plot
  setwd("..")
  dir.create("correlation_plot", showWarnings = FALSE)
  setwd("correlation_plot")  
  QCReport(data)
  png(filename = "correlation_plot_notNormalized.png")
  correlationPlot(data)
  dev.off()
  png(filename = "correlation_plot_normalized.png")
  correlationPlot(data.mas5)
  dev.off()
  #correlationPlot(as.matrix(exprs(data.rma)))



# Ende eines Experiment -> Verlasse Ordner
  print("Bearbeiten des Experimentes beendet")
  setwd("../../..")
}
