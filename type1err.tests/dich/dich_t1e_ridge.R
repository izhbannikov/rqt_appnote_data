# Power test for rqt #
library(rqt)
library(SKAT)

# how to run: source("/Users/ilya/Dropbox/rqt/AppNote/MinorRevision/tests/type1err.tests/dich/dich_t1e_ridge.R")

proj.dir <- "/Users/ilya/Projects/rqt.dev/tests/type1error.tests/dich_ridge/" # project directory, change it if you need.
n.snp <- 50
n <- 20000 # testing n times
alpha.level <- c(1e-2, 1e-3, 1e-4)
res.table <- data.frame(matrix(nrow=3, ncol=20))
colnames(res.table) <- c("rQT1", "rQT2", "rQT3", 
                         "QT1", "QT2", "QT3", 
                         "p.skat", "p.skat.o", 
                         "var.pool.rqt", "var.pool.qt", 
                         "vif.rqt", "vif.qt", 
                         "prqt1", "prqt2", "prqt3", "pqt1", "pqt2", "pqt3", "pskat", "pskato")
rownames(res.table) <- as.character(alpha.level)

res.table.case <- data.frame(matrix(nrow=3, ncol=1, 0))
colnames(res.table.case) <- c("case")
rownames(res.table.case) <- as.character(alpha.level)

ncs <- 1500 # Number of cases
nct <- 1500 # Number of controls


for(alpha in alpha.level) { 
    res.t1e <- list(rQT1=0, rQT2=0, rQT3=0, 
                    QT1=0, QT2=0, QT3=0, 
                    p.skat=0, p.skat.o=0,
                    var.pool.rqt=0, var.pool.qt=0,
                    vif.rqt=0, vif.qt=0,
                    prqt1=0, prqt2=0, prqt3=0, pqt1=0, pqt2=0, pqt3=0, pskat=0, pskato=0)
    
    pval.table.rqt <- matrix(nrow=n, ncol=3)
    pval.table.qt <- matrix(nrow=n, ncol=3)
    pval.table.skat <- matrix(nrow=n, ncol=1)
    pval.table.skato <- matrix(nrow=n, ncol=1)
    
    
    for (i in 1:n){
    
        print(paste("alpha:", alpha, "Iteration:", i))
        # Data simulation #
        # Prepare raw files with cophesim () and plink #
        system(paste("plink --simulate-ncases ", ncs, " --simulate-ncontrols ", nct, " --simulate ", proj.dir, "wgas.sim ", " --out ", proj.dir, "sim.plink ", " --make-bed >/dev/null", sep=''))
        system(paste("python /Users/ilya/Projects/cophesim_stable/cophesim.py -i ", proj.dir, "sim.plink -o ", proj.dir, "testout" , " >/dev/null", sep=''))
        system(paste("plink --file ", proj.dir, "testout_pheno_bin.txt", " --recodeA --out ", proj.dir, "testout_pheno_bin.txt >/dev/null", sep=''))
  
        # Combining Phenotype and Genotype #
        ## Dichotomous ##
        p <- read.table(paste(proj.dir, "testout_pheno_bin.txt", sep=''), header=TRUE)
        p <- p[["pheno"]]
        g <- read.table(paste(proj.dir, "testout_pheno_bin.txt.raw", sep=''), header=TRUE)
        g <- g[,7:dim(g)[2]]
        colnames(g) <- paste("snp", 1:n.snp, sep="")
        d <- cbind(pheno=p, g)
        write.table(x = d, file=paste(proj.dir, "test.bin", i, ".dat", sep=''), quote = FALSE, row.names = FALSE)
  
        res.table.case[as.character(alpha), 1] <- res.table.case[as.character(alpha), 1] + length(which(p==1))
        
        # Tests #
        ## RQT ##
        data <- data.matrix(read.table(paste(proj.dir, "test.bin", i, ".dat", sep=''), header=TRUE))
        pheno <- data[,1]
        geno <- data[, 2:dim(data)[2]]
        colnames(geno) <- paste(seq(1, dim(geno)[2]))
        geno.obj <- SummarizedExperiment(geno)
        obj <- rqt(phenotype=pheno, genotype=geno.obj)
        res <- geneTest(obj, method="ridge", out.type = "D", cumvar.threshold = 70, weight=F)
        
        ## QT ##
        res.qt <- geneTest(obj, method="none", out.type = "D", weight=F)
        
        ## SKAT ##
        obj.b<-SKAT_Null_Model(pheno ~ 1, out_type="D")
        res.skat <- SKAT(geno, obj.b)$p.value
        ## SKAT-O ##
        obj.b<-SKAT_Null_Model(pheno ~ 1, out_type="D")
        res.skat.o <- SKAT(geno, obj.b, method = "optimal.adj")$p.value
        
        
        
        ######################### Calclulating Type-I error #########################
        ## rqt ##
        if(!is.na(results(res)$pValue$pVal1)) {
          if (results(res)$pValue$pVal1 < alpha) (res.t1e$rQT1 <- res.t1e$rQT1 + 1) 
          pval.table.rqt[i,1] <- results(res)$pValue$pVal1
        } else {
          res.t1e$rQT1 <- res.t1e$rQT1 + 1
        }
        
        if(!is.na(results(res)$pValue$pVal2)) {
          if (results(res)$pValue$pVal2 < alpha) (res.t1e$rQT2 <- res.t1e$rQT2 + 1) 
          pval.table.rqt[i,2] <- results(res)$pValue$pVal2
        } else {
          res.t1e$rQT2 <- res.t1e$rQT2 + 1
        }
        
        if(!is.na(results(res)$pValue$pVal3)) {
          if (results(res)$pValue$pVal3 < alpha) (res.t1e$rQT3 <- res.t1e$rQT3 + 1) 
          pval.table.rqt[i,3] <- results(res)$pValue$pVal3
        } else {
          res.t1e$rQT3 <- res.t1e$rQT3 + 1
        }
        
        ## qt ##
        if(!is.na(results(res.qt)$pValue$pVal1)) {
          if (results(res.qt)$pValue$pVal1 < alpha) (res.t1e$QT1 <- res.t1e$QT1 + 1) 
          pval.table.qt[i,1] <- results(res.qt)$pValue$pVal1
        } else {
          res.t1e$QT1 <- res.t1e$QT1 + 1
        }
        
        if(!is.na(results(res.qt)$pValue$pVal2)) {
          if (results(res.qt)$pValue$pVal2 < alpha) (res.t1e$QT2 <- res.t1e$QT2 + 1) 
          pval.table.qt[i,2] <- results(res.qt)$pValue$pVal2
        } else {
          res.t1e$QT2 <- res.t1e$QT2 + 1
        }
        
        if(!is.na(results(res.qt)$pValue$pVal3)) {
          if (results(res.qt)$pValue$pVal3 < alpha) (res.t1e$QT3 <- res.t1e$QT3 + 1) 
          pval.table.qt[i,3] <- results(res.qt)$pValue$pVal3
        } else {
          res.t1e$QT3 <- res.t1e$QT3 + 1
        }
        
  
        ## SKAT, SKAT-O ##
        if (res.skat < alpha) (res.t1e$p.skat <- res.t1e$p.skat + 1)
        if (res.skat.o < alpha) (res.t1e$p.skat.o <- res.t1e$p.skat.o + 1)
        
        pval.table.skat[i,1] <- res.skat
        pval.table.skato[i,1] <- res.skat.o
        
        ### Pooled variance ###
        res.t1e$var.pool.rqt <- res.t1e$var.pool.rqt + results(res)$var.pooled
        res.t1e$var.pool.qt <- res.t1e$var.pool.qt + results(res.qt)$var.pooled
        
        ### Mean VIF ###
        res.t1e$vif.rqt <- res.t1e$vif.rqt + results(res)$mean.vif
        res.t1e$vif.qt <- res.t1e$vif.qt + results(res.qt)$mean.vif
        
        
        ######################### ####################### #########################
        
        # Printing results #
        print(paste("Type I error rate for p.rQT1 in percentage is", (res.t1e$rQT1/i)*100,"%"))
        print(paste("Type I error rate for p.rQT2 in percentage is", (res.t1e$rQT2/i)*100,"%"))
        print(paste("Type I error rate for p.rQT3 in percentage is", (res.t1e$rQT3/i)*100,"%"))
        ###
        print(paste("Type I error rate for p.QT1 in percentage is", (res.t1e$QT1/i)*100,"%"))
        print(paste("Type I error rate for p.QT2 in percentage is", (res.t1e$QT2/i)*100,"%"))
        print(paste("Type I error rate for p.QT3 in percentage is", (res.t1e$QT3/i)*100,"%"))
        ###
        print(paste("Type I error rate for p.skat in percentage is", (res.t1e$p.skat/i)*100,"%"))
        print(paste("Type I error rate for p.skat.o in percentage is", (res.t1e$p.skat.o/i)*100,"%"))
        ###
        
        print(paste("Avg pooled variance rQT", (res.t1e$var.pool.rqt/i)))
        print(paste("Avg pooled variance QT", (res.t1e$var.pool.qt/i)))
        ###
        print(paste("Avg mean vif rQT", (res.t1e$vif.rqt/i)))
        print(paste("Avg mean vif QT", (res.t1e$vif.qt/i)))
        ### P-values ###
        print(paste("Avg pvalue rqt1", (res.t1e$prqt1/i)))
        print(paste("Avg pvalue rqt2", (res.t1e$prqt2/i)))
        print(paste("Avg pvalue rqt3", (res.t1e$prqt3/i)))
        ###
        print(paste("Avg pvalue qt1", (res.t1e$pqt1/i)))
        print(paste("Avg pvalue qt2", (res.t1e$pqt2/i)))
        print(paste("Avg pvalue qt3", (res.t1e$pqt3/i)))
        ###
        print(paste("Avg pvalue pskat", (res.t1e$pskat/i)))
        print(paste("Avg pvalue pskato", (res.t1e$pskato/i)))
        
        
        # Cleaning up #
        rm(p, d, obj.b, res.skat.o, data, pheno, geno, obj, res, res.qt, geno.obj)
        system(paste("rm", paste(proj.dir, "test.bin", i, ".dat", sep='')))
        system(paste("rm", paste(proj.dir, "testout_pheno_bin.txt*", sep='')))
        system(paste("rm", paste(proj.dir, "sim.plink.*", sep='')))
        
    }
    
    # Printing results #
    print(paste("Type II error rate for p.rQT1 in percentage is", (res.t1e$rQT1/n)*100,"%"))
    print(paste("Type II error rate for p.rQT2 in percentage is", (res.t1e$rQT2/n)*100,"%"))
    print(paste("Type II error rate for p.rQT3 in percentage is", (res.t1e$rQT3/n)*100,"%"))
    
    print(paste("Type II error rate for p.QT1 in percentage is", (res.t1e$QT1/n)*100,"%"))
    print(paste("Type II error rate for p.QT2 in percentage is", (res.t1e$QT2/n)*100,"%"))
    print(paste("Type II error rate for p.QT3 in percentage is", (res.t1e$QT3/n)*100,"%"))
    
    print(paste("Type II error rate for p.skat in percentage is", (res.t1e$p.skat/n)*100,"%"))
    print(paste("Type II error rate for p.skat.o in percentage is", (res.t1e$p.skat.o/n)*100,"%"))
    ###
    ###
    print(paste("Avg pooled variance rQT", (res.t1e$var.pool.rqt/n)))
    print(paste("Avg pooled variance QT", (res.t1e$var.pool.qt/n)))
    ###
    print(paste("Avg mean vif rQT", (res.t1e$vif.rqt/n)))
    print(paste("Avg mean vif QT", (res.t1e$vif.qt/n)))
    print(paste("Avg pvalue rqt1", (res.t1e$prqt1/n)))
    print(paste("Avg pvalue rqt2", (res.t1e$prqt2/n)))
    print(paste("Avg pvalue rqt3", (res.t1e$prqt3/n)))
    ###
    print(paste("Avg pvalue qt1", (res.t1e$pqt1/n)))
    print(paste("Avg pvalue qt2", (res.t1e$pqt2/n)))
    print(paste("Avg pvalue qt3", (res.t1e$pqt3/n)))
    ###
    print(paste("Avg pvalue pskat", (res.t1e$pskat/n)))
    print(paste("Avg pvalue pskato", (res.t1e$pskato/n)))
    
    res.table[as.character(alpha),] <- unlist(res.t1e)
    write.table(x = res.t1e, file = paste(proj.dir, "res.t1e", alpha, ".txt", sep=""), row.names = F, quote=F)
    
    res.table.case[as.character(alpha), 1] <- res.table.case[as.character(alpha), 1]/n
    res.table[as.character(alpha), ] <- res.table[as.character(alpha), ]/n
    print(res.table)
    
    write.table(x = pval.table.rqt, file = paste(proj.dir, "pval.table.rqt.", alpha, ".txt", sep=""), row.names = F, quote=F)
    write.table(x = pval.table.qt, file = paste(proj.dir, "pval.table.qt.", alpha, ".txt", sep=""), row.names = F, quote=F)
    write.table(x = pval.table.skat, file = paste(proj.dir, "pval.table.skat.", alpha, ".txt", sep=""), row.names = F, quote=F)
    write.table(x = pval.table.skato, file = paste(proj.dir, "pval.table.skato.", alpha, ".txt", sep=""), row.names = F, quote=F)
    
}

print("!!! Final results !!!")
print(res.table)
write.table(x = res.table, file = paste(proj.dir, "res.table.t1e.txt", sep=""))
write.table(x = res.table.case, file = paste(proj.dir, "res.casecont.txt", sep=""))