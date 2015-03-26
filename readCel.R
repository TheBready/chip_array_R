######################################
#       Einlesen .CEL-Dateien        #
#               von                  #
#     Nadine, Felix und Philipp      #
#             Gruppe 2               #
######################################
# Installieren der benötigten Pakete
#source("http://bioconductor.org/biocLite.R")
#biocLite("affy")

# Laden von affy
library("affy")

#Einstellen des Pfades
setwd("../input")
dir <- dir()

# Laden der .cel Files als Batch (alle in einem Ordner)
data <- ReadAffy(celfile.path=dir,verbose = TRUE)

#Namen der verschiedenen Samples
CELnames <- colnames(data)

#Erstellen der Bilder der Chips
PNGnames <- gsub('.{3}$', 'png', CELnames)
setwd("../output")
for(i in 1:length(PNGnames)){
  png(filename=PNGnames[i], width = 1164, height = 1164, units = "px")
  image(data[,i])
  dev.off()  
}


#Ausgabe der Rohdaten in .txt
data.exp <- exprs(data)
data.exp <- cbind(rownames(data.exp),data.exp)
write(t(data.exp), file = gsub('.{0}$', '_signals.txt', dir),ncolumns = length(CELnames))


#Ausgabe der perfect match Rohdaten in .txt
data.pm <- pm(data)
data.pm <- cbind(rownames(data.pm),data.pm)
data.pm <- data.pm [order(as.numeric(data.pm[,1])),]
write(t(data.pm), file = gsub('.{0}$', '_signals_PM.txt', dir),ncolumns = length(CELnames))

#Ausgabe der mismatch Rohdaten in .txt
data.mm <- mm(data)
data.mm <- cbind(rownames(data.mm),data.mm)
data.mm <- data.mm [order(as.numeric(data.mm[,1])),]
write(t(data.mm), file = gsub('.{0}$', '_signals_MM.txt', dir),ncolumns = length(CELnames))


#RMA
data.rma <- rma(data)
write.exprs(data.rma,file = gsub('.{0}$', '_signals_RMA.txt', dir))

#MAS 5.0
data.mas5 <- mas5(data, sc = 500)
data.mas5exp <- exprs(data.mas5)
data.mas5calls <- mas5calls(data)
write.exprs(data.mas5,file = gsub('.{0}$', '_MAS5_500.txt', dir))



