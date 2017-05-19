# QQ Plots under Null Hyphotesis

source("~/Dropbox/rqt/AppNote/MinorRevision/tests/type1err.tests/QQ_plots/qq.plot.null.R")

alpha.levels <- c(0.01, 0.001, 1e-4)
methods <- c("pca", "pls", "lasso", "ridge")
outcomes <- c("cont", "dich")

dir <- "~/Projects/rqt.dev/tests/type1error.tests/"


for(method in methods) {
    for(outcome in outcomes) {
        if (outcome == "dich" & method == "pls") 
                next
        
        png(filename=paste(dir, "QQ_plot_", outcome, "_", method, ".png", sep=""), width = 2700, height = 1024)
        par(mfrow=c(3,8), mgp=c(2, 0, 0))
        
        # Alpha 0.01
        for(i in 1:3) {
            d <- read.table(paste(dir, outcome, "_", method, "/pval.table.rqt.0.01.txt", sep=""), header = TRUE)
            ggd.qqplot(d[,i], main=paste("RQT",i,sep=""), cex.lab=2, cex.main=2)
        }
        for(i in 1:3) {
            d <- read.table(paste(dir, outcome, "_", method, "/pval.table.qt.0.01.txt", sep=""), header = TRUE)
            ggd.qqplot(d[,i], main=paste("QT",i,sep=""), cex.lab=2, cex.main=2)
        }
        d <- read.table(paste(dir, outcome, "_", method, "/pval.table.skat.0.01.txt", sep=""), header = TRUE)
        ggd.qqplot(d[,1], main="SKAT", cex.lab=2, cex.main=2)
        d <- read.table(paste(dir, outcome, "_", method, "/pval.table.skato.0.01.txt", sep=""), header = TRUE)
        ggd.qqplot(d[,1], main="SKAT-O", cex.lab=2, cex.main=2)
        
        mtext('Alpha level = 0.01', side=3, line=-2, outer=TRUE, cex=2)
        
        # Alpha 0.001
        for(i in 1:3) {
          d <- read.table(paste(dir, outcome, "_", method, "/pval.table.rqt.0.001.txt", sep=""), header = TRUE)
          ggd.qqplot(d[,i], main=paste("RQT",i,sep=""), cex.lab=2, cex.main=2)
        }
        for(i in 1:3) {
          d <- read.table(paste(dir, outcome, "_", method, "/pval.table.qt.0.001.txt", sep=""), header = TRUE)
          ggd.qqplot(d[,i], main=paste("QT",i,sep=""), cex.lab=2, cex.main=2)
        }
        d <- read.table(paste(dir, outcome, "_", method, "/pval.table.skat.0.001.txt", sep=""), header = TRUE)
        ggd.qqplot(d[,1], main="SKAT", cex.lab=2, cex.main=2)
        d <- read.table(paste(dir, outcome, "_", method, "/pval.table.skato.0.001.txt", sep=""), header = TRUE)
        ggd.qqplot(d[,1], main="SKAT-O", cex.lab=2, cex.main=2)
        
        mtext('Alpha level = 0.001', side=3, line=-38, outer=TRUE, cex=2) 
        
        # Alpha 1e-04
        for(i in 1:3) {
          d <- read.table(paste(dir, outcome, "_", method, "/pval.table.rqt.1e-04.txt", sep=""), header = TRUE)
          ggd.qqplot(d[,i], main=paste("RQT",i,sep=""), cex.lab=2, cex.main=2)
        }
        for(i in 1:3) {
          d <- read.table(paste(dir, outcome, "_", method, "/pval.table.qt.1e-04.txt", sep=""), header = TRUE)
          ggd.qqplot(d[,i], main=paste("QT",i,sep=""), cex.lab=2, cex.main=2)
        }
        d <- read.table(paste(dir, outcome, "_", method, "/pval.table.skat.1e-04.txt", sep=""), header = TRUE)
        ggd.qqplot(d[,1], main="SKAT", cex.lab=2, cex.main=2)
        d <- read.table(paste(dir, outcome, "_", method, "/pval.table.skato.1e-04.txt", sep=""), header = TRUE)
        ggd.qqplot(d[,1], main="SKAT-O", cex.lab=2, cex.main=2)
        
        mtext('Alpha level = 1E-04', side=3, line=-73, outer=TRUE, cex=2)
        
        #mtext(paste("QQ plots under null hyphotesis for alpha level", alpha, "method", method, "and", outcome, "outcome"), side = 3, line=-1, outer = TRUE)
        dev.off()
    }
}
