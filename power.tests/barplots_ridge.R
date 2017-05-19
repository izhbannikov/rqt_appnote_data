# Power bar plots #
n <- 1000

pdir <- "~/Projects/rqt.dev/tests/power.tests/power.ld.dich_ridge/"

data <- read.table(paste(pdir, "res.table.pow.txt", sep=''), header=T)
data <- t(1-data[,1:8])
rownames(data) <- c("rQT1", "rQT2", "rQT3", "QT1", "QT2", "QT3", "SKAT", "SKAT-O")
colnames(data) <- c("10", "25")

png(width = 1024, height = 600, filename = paste(pdir, "barplots.png", sep=''))
bp<-barplot(data, col=rainbow(dim(data)[1]), 
            beside=TRUE,cex.names = T,
            ylim = c(0,1.15), ylab= "Power", 
            xlab="% of causal SNPs", main = "Power estimates for dichotomous phenotype with presence of LD\nMethod: rigde regression")

text(bp, c(data), 
     c("RQT1", "RQT2", "RQT3", "QT1", "QT2", "QT3", "SKAT", "SKAT-O",
       "RQT1", "RQT2", "RQT3", "QT1", "QT2", "QT3", "SKAT", "SKAT-O"),
     cex=1.2, pos=4, srt = 90)
dev.off()        

pdir <- "~/Projects/rqt.dev/tests/power.tests/power.ld.cont_ridge/"
data <- read.table(paste(pdir, "res.table.pow.txt", sep=''), header=T)
data <- t(1-data[,1:8])
rownames(data) <- c("rQT1", "rQT2", "rQT3", "QT1", "QT2", "QT3", "SKAT", "SKAT-O")
colnames(data) <- c("10", "25")

png(width = 1024, height = 600, filename = paste(pdir, "barplots.png", sep=''))
bp<-barplot(data, col=rainbow(dim(data)[1]), 
            beside=TRUE,cex.names = T,
            ylim = c(0,1.15), ylab= "Power", 
            xlab="% of causal SNPs", main = "Power estimates for continuous phenotype with presence of LD\nMethod: ridge regression")

text(bp, c(data), 
     c("RQT1", "RQT2", "RQT3", "QT1", "QT2", "QT3", "SKAT", "SKAT-O",
       "RQT1", "RQT2", "RQT3", "QT1", "QT2", "QT3", "SKAT", "SKAT-O"),
     cex=1.2, pos=4, srt = 90)
dev.off()


