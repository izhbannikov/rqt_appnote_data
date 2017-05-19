# Power test for rqt #
library(rqt)
library(SKAT)

# source("~/Dropbox/rqt/AppNote/MajorRevision/tests/power.tests/power.ld.cont/power.ld.cont_pls.R")

proj.dir <- "/Users/ilya/Projects/rqt.dev/tests/power.tests/power.ld.cont_pls/" # project directory, change it if you need.
n.snp <- 50
n <- 1000 # testing 1,000 times
alpha.level <- 1e-3
res.table <- data.frame(matrix(nrow=2, ncol=20))
colnames(res.table) <- c("rQT1", "rQT2", "rQT3", 
                         "QT1", "QT2", "QT3", 
                         "p.skat", "p.skat.o", 
                         "var.pool.rqt", "var.pool.qt", 
                         "vif.rqt", "vif.qt", 
                         "prqt1", "prqt2", "prqt3", "pqt1", "pqt2", "pqt3", "pskat", "pskato")

rownames(res.table) <- c("0.1","0.25")


res.table.case <- data.frame(matrix(nrow=2, ncol=1, 0))
colnames(res.table.case) <- c("case")
rownames(res.table.case) <- c("0.1","0.25")

ncs <- 1500 # Number of cases
nct <- 1500 # Number of controls

caus.var <- c(0.1, 0.25)

set.seed(100500)

for(s in caus.var) {
    res.pow <- list(rQT1=0, rQT2=0, rQT3=0, 
                    QT1=0, QT2=0, QT3=0, 
                    p.skat=0, p.skat.o=0,
                    var.pool.rqt=0, var.pool.qt=0,
                    vif.rqt=0, vif.qt=0,
                    prqt1=0, prqt2=0, prqt3=0, pqt1=0, pqt2=0, pqt3=0, pskat=0, pskato=0)
  
   
    n.caus <- round(n.snp*s,digits = 0)
    
    eff <- sapply(0:(n.caus-1), function(n) { paste(n, ':', ifelse( (runif(1) > 0.5), -0.2, 0.2 ), sep='') })
    write(x = eff, paste(proj.dir, "effects.txt", sep=""))
    
    
    ld <- c() 
    for(j in 0:(n.snp-2)) {
      if(j %% 4) {
        ld <- c(ld, paste( j, ',', j+1, ',', 0.95, sep=''))
      }
    }
    write(x = ld, paste(proj.dir, "ldfile.txt", sep=""))
  
    for (i in 1:n){
        
        print(paste("n.caus.snp:", n.caus, "Iteration:", i))
        # Data simulation #
        # Prepare raw files with cophesim () and plink #
        system(paste("plink --simulate-ncases ", ncs, " --simulate-ncontrols ", nct, " --simulate ", proj.dir, "wgas.sim ", " --out ", proj.dir, "sim.plink ", " --make-bed >/dev/null", sep=''))
        system(paste("python /Users/ilya/Projects/cophesim_stable/cophesim.py -i ", proj.dir, "sim.plink -o ", proj.dir, "testout" , " -c -ce ",  proj.dir, "effects.txt", " -LD ", proj.dir, "ldfile.txt", " >/dev/null", sep=''))
        system(paste("plink --file ", proj.dir, "testout_pheno_cont.txt", " --recodeA --out ", proj.dir, "testout_pheno_cont.txt >/dev/null", sep=''))
  
        # Combining Phenotype and Genotype #
        ## Dichotomous ##
        p <- read.table(paste(proj.dir, "testout_pheno_cont.txt", sep=''), header=TRUE)
        p <- p[["pheno"]]
        g <- read.table(paste(proj.dir, "testout_pheno_cont.txt.raw", sep=''), header=TRUE)
        g <- g[,7:dim(g)[2]]
        colnames(g) <- paste("snp", 1:n.snp, sep="")
        d <- cbind(pheno=p, g)
        write.table(x = d, file=paste(proj.dir, "test.cont", i, ".dat", sep=''), quote = FALSE, row.names = FALSE)
  
        res.table.case[as.character(s), 1] <- res.table.case[as.character(s), 1] + length(which(p==1))
        
        # Tests #
        ## RQT ##
        data <- data.matrix(read.table(paste(proj.dir, "test.cont", i, ".dat", sep=''), header=TRUE))
        pheno <- data[,1]
        geno <- data[, 2:dim(data)[2]]
        colnames(geno) <- paste(seq(1, dim(geno)[2]))
        geno.obj <- SummarizedExperiment(geno)
        obj <- rqt(phenotype=pheno, genotype=geno.obj)
        res <- geneTest(obj, method="pls", out.type = "C", cumvar.threshold = 10, weight=F)
        print(paste("p.rQT1", results(res)$pValue$pVal1))
        print(paste("p.rQT2", results(res)$pValue$pVal2))
        print(paste("p.rQT3", results(res)$pValue$pVal3))
        print(paste("Var pooled rQT: ", results(res)$var.pooled))
        print(paste("Mean vif rQT: ", results(res)$mean.vif))
        ## QT ##
        res.qt <- geneTest(obj, method="none", out.type = "C", weight=F)
        print(paste("p.QT1", results(res.qt)$pValue$pVal1))
        print(paste("p.QT2", results(res.qt)$pValue$pVal2))
        print(paste("p.QT3", results(res.qt)$pValue$pVal3))
        print(paste("Var pooled QT: ", results(res.qt)$var.pooled))
        print(paste("Mean vif QT: ", results(res.qt)$mean.vif))
        ## SKAT ##
        obj.b<-SKAT_Null_Model(pheno ~ 1, out_type="C")
        res.skat <- SKAT(geno, obj.b)$p.value
        print(paste("skat", res.skat))
        ## SKAT-O ##
        obj.b<-SKAT_Null_Model(pheno ~ 1, out_type="C")
        res.skat.o <- SKAT(geno, obj.b, method = "optimal.adj")$p.value
        print(paste("skat.o", res.skat.o))
        
        # Calculating power #
        ### rQT ###
        if(!is.na(results(res)$pValue$pVal1)) {
          if (results(res)$pValue$pVal1 > alpha.level) (res.pow$rQT1 <- res.pow$rQT1 + 1) 
          res.pow$prqt1 <- res.pow$prqt1 + results(res)$pValue$pVal1
        } else {
          res.pow$rQT1 <- res.pow$rQT1 + 1
          res.pow$prqt1 <- res.pow$prqt1 + 1
        }
        
        if(!is.na(results(res)$pValue$pVal2)) {
          if (results(res)$pValue$pVal2 > alpha.level) (res.pow$rQT2 <- res.pow$rQT2 + 1) 
          res.pow$prqt2 <- res.pow$prqt2 + results(res)$pValue$pVal2
        } else {
          res.pow$rQT2 <- res.pow$rQT2 + 1
          res.pow$prqt2 <- res.pow$prqt2 + 1
        }
        
        if(!is.na(results(res)$pValue$pVal3)) {
          if (results(res)$pValue$pVal3 > alpha.level) (res.pow$rQT3 <- res.pow$rQT3 + 1) 
          res.pow$prqt3 <- res.pow$prqt3 + results(res)$pValue$pVal3
        } else {
          res.pow$rQT3 <- res.pow$rQT3 + 1
          res.pow$prqt3 <- res.pow$prqt3 + 1
        }
        
        ### QT ###
        if(!is.na(results(res.qt)$pValue$pVal1)) {
          if (results(res.qt)$pValue$pVal1 > alpha.level) (res.pow$QT1 <- res.pow$QT1 + 1) 
          res.pow$pqt1 <- res.pow$pqt1 + results(res.qt)$pValue$pVal1
        } else {
          res.pow$QT1 <- res.pow$QT1 + 1
          res.pow$pqt1 <- res.pow$pqt1 + 1
        }
        
        if(!is.na(results(res.qt)$pValue$pVal2)) {
          if (results(res.qt)$pValue$pVal2 > alpha.level) (res.pow$QT2 <- res.pow$QT2 + 1) 
          res.pow$pqt2 <- res.pow$pqt2 + results(res.qt)$pValue$pVal2
        } else {
          res.pow$QT2 <- res.pow$QT2 + 1
          res.pow$pqt2 <- res.pow$pqt2 + 1
        }
        
        if(!is.na(results(res.qt)$pValue$pVal3)) {
          if (results(res.qt)$pValue$pVal3 > alpha.level) (res.pow$QT3 <- res.pow$QT3 + 1) 
          res.pow$pqt3 <- res.pow$pqt3 + results(res.qt)$pValue$pVal3
        } else {
          res.pow$QT3 <- res.pow$QT3 + 1
          res.pow$pqt3 <- res.pow$pqt3 + 1
        }
        
        
        ### SKAT/SKAT-O ###
        if (res.skat > alpha.level) (res.pow$p.skat <- res.pow$p.skat + 1)
        res.pow$pskat <- res.pow$pskat + res.skat
        if (res.skat.o > alpha.level) (res.pow$p.skat.o <- res.pow$p.skat.o + 1)
        res.pow$pskato <- res.pow$pskato + res.skat.o
        
        ### Pooled variance ###
        res.pow$var.pool.rqt <- res.pow$var.pool.rqt + results(res)$var.pooled
        res.pow$var.pool.qt <- res.pow$var.pool.qt + results(res.qt)$var.pooled
        
        ### Mean VIF ###
        res.pow$vif.rqt <- res.pow$vif.rqt + results(res)$mean.vif
        res.pow$vif.qt <- res.pow$vif.qt + results(res.qt)$mean.vif
        
        
        # Printing results #
        print(paste("Type II error rate for p.rQT1 in percentage is", (res.pow$rQT1/i)*100,"%"))
        print(paste("Type II error rate for p.rQT2 in percentage is", (res.pow$rQT2/i)*100,"%"))
        print(paste("Type II error rate for p.rQT3 in percentage is", (res.pow$rQT3/i)*100,"%"))
        ###
        print(paste("Type II error rate for p.QT1 in percentage is", (res.pow$QT1/i)*100,"%"))
        print(paste("Type II error rate for p.QT2 in percentage is", (res.pow$QT2/i)*100,"%"))
        print(paste("Type II error rate for p.QT3 in percentage is", (res.pow$QT3/i)*100,"%"))
        ###
        print(paste("Type II error rate for p.skat in percentage is", (res.pow$p.skat/i)*100,"%"))
        print(paste("Type II error rate for p.skat.o in percentage is", (res.pow$p.skat.o/i)*100,"%"))
        ###
        ###
        
        print(paste("Avg pooled variance rQT", (res.pow$var.pool.rqt/i)))
        print(paste("Avg pooled variance QT", (res.pow$var.pool.qt/i)))
        ###
        print(paste("Avg mean vif rQT", (res.pow$vif.rqt/i)))
        print(paste("Avg mean vif QT", (res.pow$vif.qt/i)))
        ### P-values ###
        print(paste("Avg pvalue rqt1", (res.pow$prqt1/i)))
        print(paste("Avg pvalue rqt2", (res.pow$prqt2/i)))
        print(paste("Avg pvalue rqt3", (res.pow$prqt3/i)))
        ###
        print(paste("Avg pvalue qt1", (res.pow$pqt1/i)))
        print(paste("Avg pvalue qt2", (res.pow$pqt2/i)))
        print(paste("Avg pvalue qt3", (res.pow$pqt3/i)))
        ###
        print(paste("Avg pvalue pskat", (res.pow$pskat/i)))
        print(paste("Avg pvalue pskato", (res.pow$pskato/i)))
        
        
        # Cleaning up #
        rm(p, d, obj.b, res.skat.o, data, pheno, geno, obj, res, res.qt, geno.obj)
        system(paste("rm", paste(proj.dir, "test.cont", i, ".dat", sep='')))
        system(paste("rm", paste(proj.dir, "testout_pheno_cont.txt*", sep='')))
        system(paste("rm", paste(proj.dir, "sim.plink.*", sep='')))
        
    }
    
    # Printing results #
    print(paste("Type II error rate for p.rQT1 in percentage is", (res.pow$rQT1/n)*100,"%"))
    print(paste("Type II error rate for p.rQT2 in percentage is", (res.pow$rQT2/n)*100,"%"))
    print(paste("Type II error rate for p.rQT3 in percentage is", (res.pow$rQT3/n)*100,"%"))
    
    print(paste("Type II error rate for p.QT1 in percentage is", (res.pow$QT1/n)*100,"%"))
    print(paste("Type II error rate for p.QT2 in percentage is", (res.pow$QT2/n)*100,"%"))
    print(paste("Type II error rate for p.QT3 in percentage is", (res.pow$QT3/n)*100,"%"))
    
    print(paste("Type II error rate for p.skat in percentage is", (res.pow$p.skat/n)*100,"%"))
    print(paste("Type II error rate for p.skat.o in percentage is", (res.pow$p.skat.o/n)*100,"%"))
    ###
    print(paste("Avg pooled variance rQT", (res.pow$var.pool.rqt/n)))
    print(paste("Avg pooled variance QT", (res.pow$var.pool.qt/n)))
    ###
    print(paste("Avg mean vif rQT", (res.pow$vif.rqt/n)))
    print(paste("Avg mean vif QT", (res.pow$vif.qt/n)))
    print(paste("Avg pvalue rqt1", (res.pow$prqt1/n)))
    print(paste("Avg pvalue rqt2", (res.pow$prqt2/n)))
    print(paste("Avg pvalue rqt3", (res.pow$prqt3/n)))
    ###
    print(paste("Avg pvalue qt1", (res.pow$pqt1/n)))
    print(paste("Avg pvalue qt2", (res.pow$pqt2/n)))
    print(paste("Avg pvalue qt3", (res.pow$pqt3/n)))
    ###
    print(paste("Avg pvalue pskat", (res.pow$pskat/n)))
    print(paste("Avg pvalue pskato", (res.pow$pskato/n)))
    
    res.table[as.character(s),] <- unlist(res.pow)
    write.table(x = res.pow, file = paste(proj.dir, "res.pow", s, ".txt", sep=""), row.names = F, quote=F)
    
    res.table.case[as.character(s), 1] <- res.table.case[as.character(s), 1]/n
    res.table[as.character(s), ] <- res.table[as.character(s), ]/n
    print(res.table)
}

print("!!! Final results !!!")
print(res.table)
write.table(x = res.table, file = paste(proj.dir, "res.table.pow.txt", sep=""))
write.table(x = res.table.case, file = paste(proj.dir, "res.casecont.txt", sep=""))