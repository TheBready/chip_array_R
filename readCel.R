
# Installieren der benötigten Pakete
#test
#source("http://bioconductor.org/biocLite.R")
#biocLite("affy")

# Laden von affy
library("affy", lib.loc="~/R/win-library/3.1")

setwd("C:/Users/Felix/OwnCloud/Studium/7. Fachsemester/Software-Praktikum/R")

# Anezigen der vorandenden .cel Dateien
list.celfiles(path="input/ND_Group2_133Plus_2",full.names=TRUE)

# Laden der .cel Files als Batch (alle in einem Ordner)
data <- ReadAffy(celfile.path="input/ND_Group2_133Plus_2",verbose = TRUE)
data.probes <- probes(data)




# Histogramme der einzelnen Chips
p1 <-hist(log(intensity(data[, 1])), breaks = 100)
p2 <-hist(log(intensity(data[, 2])), breaks = 100)
p3 <-hist(log(intensity(data[, 3])), breaks = 100)
p4 <-hist(log(intensity(data[, 4])), breaks = 100)
p5 <-hist(log(intensity(data[, 5])), breaks = 100)
p6 <-hist(log(intensity(data[, 6])), breaks = 100)
#p7 <-hist(log(intensity(data[, 7])), breaks = 100)

png(filename="output/hist.png")
plot( p1 ,border = F, col=rgb(0, 1, 0, 0.5), main="Histogramme der DNA-Arrays",ylab="Anzahl",xlab="Intensität(log)",ylim=c(0,130000)) 
plot( p2,col=rgb(1, 0, 0, 0.5), add=T,border = F)  
plot( p3,col=rgb(0, 0, 1, 0.5), add=T,border = F) 
plot( p4,col=rgb(1, 1, 0, 0.5), add=T,border = F) 
plot( p5,col=rgb(0, 1, 1, 0.5), add=T,border = F) 
plot( p6,col=rgb(1, 0, 1, 0.5), add=T,border = F) 
#plot( p7,col=rgb(0, 0, 0, 0.5), add=T,border = F) 
plocol <- c(rgb(0, 1, 0, 0.5),rgb(1, 0, 0, 0.5),rgb(0, 0, 1, 0.5),rgb(1, 1, 0, 0.5),rgb(0, 1, 1, 0.5),rgb(1, 0, 1, 0.5))#,rgb(0, 0, 0, 0.5)
legend('topright',colnames(data),fill = plocol, bty = 'n',border = NA)
dev.off()


# Bilder zum Vergleich 
png(filename="output/Image_1.png", width = 1164, height = 1164, units = "px")
image(data[,1])
dev.off()
png(filename="output/Image_2.png", width = 1164, height = 1164, units = "px")
image(data[,2])
dev.off()
png(filename="output/Image_3.png", width = 1164, height = 1164, units = "px")
image(data[,3])
dev.off()
png(filename="output/Image_4.png", width = 1164, height = 1164, units = "px")
image(data[,4])
dev.off()
png(filename="output/Image_5.png", width = 1164, height = 1164, units = "px")
image(data[,5])
dev.off()
png(filename="output/Image_6.png", width = 1164, height = 1164, units = "px")
image(data[,6])
dev.off()
#png(filename="output/Image_7.png", width = 1164, height = 1164, units = "px")
#image(data[,7])
#dev.off()

# Boxplot der Daten
png(filename="output/boxplot.png")
boxplot(data, col="red")
dev.off()

#Ausgabe der Rohdaten in .txt
data.exp <- exprs(data)
data.exp <- cbind(rownames(data.exp),data.exp)
write(t(data.exp), file = "output/ND_Group_133Plus_2_signals.txt",ncolumns = 7)


#Ausgabe der perfect match Rohdaten in .txt
data.pm <- pm(data)
data.pm <- cbind(rownames(data.pm),data.pm)
data.pm <- data.pm [order(as.numeric(data.pm[,1])),]
write(t(data.pm), file = "output/ND_Group_133Plus_2_signals_PM.txt",ncolumns = 7)

#Ausgabe der mismatch Rohdaten in .txt
data.mm <- mm(data)
data.mm <- cbind(rownames(data.mm),data.mm)
data.mm <- data.mm [order(as.numeric(data.mm[,1])),]
write(t(data.mm), file = "output/ND_Group_133Plus_2_signals_MM.txt",ncolumns = 7)


#hieraisches clustering exp.head exp.d
data.dist = as.matrix(t(data.exp[,2:7]))
rownames(data.dist)= c("1","2","3","4","5","6")#,"7"
data.dist = dist(data.dist,method="euclidean")
data.cluster = hclust(data.dist, method="average" )
png(filename="output/hc_tree.png")
plot(data.cluster, main= "hieraisches Clustering der Daten", xlab="Distanz", ylab="Höhe")
dev.off()

#RMA

data.rma <- rma(data)
write.exprs(data.rma,file="output/ND_Group_133Plus_2_signals_RMA.txt")

#MAS 5.0

data.mas5 <- mas5(data, sc = 500)
data.mas5exp <- exprs(data.mas5)
data.mas5calls <- mas5calls(data)
write.exprs(data.mas5,file="output/ND_Group_133Plus_2_signals_MAS5_500.txt")

# Histogramme der einzelnen Chips
p1 <-hist(log(data.mas5exp[, 1]), breaks = 100)
p2 <-hist(log(data.mas5exp[, 2]), breaks = 100)
p3 <-hist(log(data.mas5exp[, 3]), breaks = 100)
p4 <-hist(log(data.mas5exp[, 4]), breaks = 100)
p5 <-hist(log(data.mas5exp[, 5]), breaks = 100)
p6 <-hist(log(data.mas5exp[, 6]), breaks = 100)
#p7 <-hist(log(intensity(data[, 7])), breaks = 100)

png(filename="output/histmas5.png")
plot( p1 ,border = F, col=rgb(0, 1, 0, 0.5), main="Histogramme der DNA-Arrays",ylab="Anzahl",xlab="Intensität(log)") 
plot( p2,col=rgb(1, 0, 0, 0.5), add=T,border = F)  
plot( p3,col=rgb(0, 0, 1, 0.5), add=T,border = F) 
plot( p4,col=rgb(1, 1, 0, 0.5), add=T,border = F) 
plot( p5,col=rgb(0, 1, 1, 0.5), add=T,border = F) 
plot( p6,col=rgb(1, 0, 1, 0.5), add=T,border = F) 
#plot( p7,col=rgb(0, 0, 0, 0.5), add=T,border = F) 
plocol <- c(rgb(0, 1, 0, 0.5),rgb(1, 0, 0, 0.5),rgb(0, 0, 1, 0.5),rgb(1, 1, 0, 0.5),rgb(0, 1, 1, 0.5),rgb(1, 0, 1, 0.5))#,rgb(0, 0, 0, 0.5)
legend('topright',colnames(data),fill = plocol, bty = 'n',border = NA)
dev.off()

# Unser Gen: 1554934_at  72.814280355092	60.2198997878917	6.00507904239215	33.030183674161	44.3790478387515	79.6205949766544


gn <- featureNames(data)
ps <- probeset(data, genenames=c("1554934_at"))
pms <- pm(ps[[1]])
mms <- mm(ps[[1]])
pmnames <- rownames(pms)
mmnames <- rownames(mms)

biocLite("hgu133plus2.db") # Chip-Datenbank
seq <- read.table("input/1554934_at.txt")
