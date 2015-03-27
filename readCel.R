######################################
#       Einlesen von  .CEL-Dateien   #
#               von                  #
#     Nadine, Felix und Philipp      #
#             Gruppe 2               #
######################################
# Installieren der benötigten Pakete
#source("http://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("hgu133plus2.db") # Chip-Datenbank


# Laden von affy
library("affy")
library("hgu133plus2.db")

#Einstellen des Pfades
setwd("../input")
dir <- dir()

# Laden der .cel Files als Batch (alle in einem Ordner)
data <- ReadAffy(celfile.path=dir,verbose = TRUE)

#Namen der verschiedenen Samples
CELnames <- colnames(data)

#Erstellen des Output Ordners
setwd("..")
dir.create("output")
setwd("output")
dir.create(dir)
setwd(dir)
dir.create("histograms")
setwd("histograms")

#Erstellen der Histgramme
PNGnames <- gsub('.{3}$', 'png', CELnames)
colors <- rainbow(length(CELnames))
for(i in 1:length(CELnames)){
  png(filename= PNGnames[i])
  hist(log(intensity(data[, i])), breaks = 100,border = F, col=colors[i], main=CELnames[i],ylab="Anzahl",xlab="Intensität(log)")
  dev.off()
}



# Erstellen der Bilder
setwd("..")
dir.create("images")
setwd("images")
for(i in 1:length(PNGnames)){
  png(filename=PNGnames[i], width = 7500, height = 7500, units = "px")
  image(data[,i])
  dev.off()  
}

# Erstellen der BoxPlots
setwd("..")
dir.create("boxplot")
setwd("boxplot")
png(filename="boxplot.png")
boxplot(data, col="red")
dev.off()  


# erstelle Liste mit allen Probe ids und Gen ids
data.proGen <- toTable(hgu133plus2SYMBOL)

#Ausgabe der Rohdaten in .txt
setwd("..")
dir.create("exprs")
setwd("exprs")
data.exp <- probes(data)
write.table(data.exp, file = gsub('.{0}$', '_signals.txt', dir), row.names=TRUE)

#Ausgabe der perfect match Rohdaten in .txt
setwd("..")
dir.create("pm")
setwd("pm")
data.pm <- pm(data,data.proGen[,1])
write.table(data.pm, file = gsub('.{0}$', '_signals_PM.txt', dir), row.names=TRUE)

#Ausgabe der mismatch Rohdaten in .txt
setwd("..")
dir.create("mm")
setwd("mm")
data.mm <- mm(data,data.proGen[,1])
write.table(data.mm, file = gsub('.{0}$', '_signals_MM.txt', dir), row.names=TRUE)


#RMA
setwd("..")
dir.create("RMA")
setwd("RMA")
data.rma <- rma(data)
write.exprs(data.rma,file = gsub('.{0}$', '_signals_RMA.txt', dir))

#MAS 5.0
setwd("..")
dir.create("MAS5")
setwd("MAS5")
data.mas5 <- mas5(data, sc = 500)
data.mas5exp <- exprs(data.mas5)
data.mas5calls <- mas5calls(data)
write.exprs(data.mas5,file = gsub('.{0}$', '_MAS5_500.txt', dir))


