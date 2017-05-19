#### Analysis for GSA 2016 ####

library(openxlsx)
library(qvalue)
library(metap)
library(plyr)
library(rqt)

#==========================================================================================================#
table_filename <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/AD_genes_GWAS_catalog_NHGRI_June2016_Ilya.xlsx"
genes_to_write_filename <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/genes.txt"
gnamesfile <- "W:/data/baru/reference-data/preprocessed/Ensembl/release-82/other/egid_to_symb.csv"
#==========================================================================================================#
gene_dir_chs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/genes_chs"
gene_dir_fhs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/genes_fhs"
gene_dir_loadfs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/genes_loadfs"
#==========================================================================================================#
pheno.data.dir <- "W:/data/processed"
#==========================================================================================================#
indiv_chs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/individuals_chs.txt"
pheno_bin_chs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/pheno_bin_chs.txt"
indiv_fhs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/individuals_fhs.txt"
pheno_bin_fhs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/pheno_bin_fhs.txt"
indiv_loadfs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/individuals_loadfs.txt"
pheno_bin_loadfs <- "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/pheno_bin_loadfs.txt"
  #==========================================================================================================#

pheno.data.dir_chs <- "/CHS/dbGaP_current/pheno_CARe/consent_ALL_ALL/current/other/diseases/ICD9/chs_AD.csv"
pheno.data.dir_llfs <- "/LLFS/LLFS_current/pheno/consent_ALL_ALL/current/other/diseases/AD.csv"
pheno.data.dir_fhs <- "/FRAM/dbGaP_current/pheno_CARe/consent_ALL_ALL/current/other/diseases/AD_dementia/FHS_Dementia_mod_AD.csv"
pheno.data.dir_loadfs <- "/LOADFS/dbGaP_current/pheno/consent_ALZ_ALL/current/LOADFS_ALZ_ALL_wide.csv"


###################################################### FUNCTIONS ############################################################

trim.leading <- function (x)  sub("^\\s+", "", x)
# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)
# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


loadGeneNames <- function(gnamesfile) { #{{{
  cat(sprintf("-- Loading gene names from %s...",gnamesfile))
  gn_df <- read.delim(gnamesfile, header=F, stringsAsFactors=F, col.names=c("EGID","GSymb"))
  genenames <- gn_df$GSymb
  names(genenames) <- gn_df$EGID
  genenames
}#}}}

loadInputFileNames <- function(geno_datadir) { #{{{
  cat(sprintf("-- Reading names of gene data files from %s...",geno_datadir))
  ifilenames <-  list.files(geno_datadir, pattern="^ENSG\\d*\\.txt$")
  names(ifilenames) <- sub(".txt", "", ifilenames, fixed=T)
  cat(sprintf(" Done, %d files.\n", length(ifilenames)))
  ifilenames
}#}}}


loadGenTable <- function(geno_datadir, fn) {
  ffn <- file.path(geno_datadir, fn)
  lines <- readLines(ffn)
  nlines <- length(lines)
  snp_ids <- unlist(strsplit(lines[1], ",", fixed=T))
  
  data <- matrix(unlist(strsplit(lines[2:nlines], "\t", fixed=T)), byrow=T, nrow=nlines-1, ncol=3)
  ids <- data.frame(data.frame(x=data[,1:2]))
  snps <- data[,3]
  nsnps <- length(snp_ids)
  
  data_matrix <- matrix(
    suppressWarnings(as.numeric(unlist(strsplit(snps[1:(nlines-1)], "", fixed=T)))),
    byrow=T, nrow=nlines-1, ncol=nsnps 
  )
  
  dm <- cbind(ids,data_matrix)
  colnames(dm) <- c("FAM_ID", "ID", snp_ids)
  dm <- dm[complete.cases(dm),] # Return a logical vector indicating which cases are complete, i.e., have no missing values.
  
  list(dm, snp_ids, fn)
  
}

calcrQTest <- function(mat, d) {
  
  #mat <- loadGenTable(dd$gene_dir, dd$ifilenames[27]) 
  snps <- mat[[1]]
  
  if(!is.null(d$covariates)) {
    data <- join_all(list(snps, d$covariates, d$pheno), by=c("FAM_ID","ID"), type="inner")
    g <- as.matrix(data[,3:(dim(data)[2]-2)])
    p <- as.matrix(data$IsIncid)
    cv <- data.frame(data$cov1)
    geno.obj <- SummarizedExperiment(g)
    obj <- rqtClass(phenotype=p, genotype=geno.obj, covariates = cv)
  } else {
    data <- join(snps, d$pheno, by=c("FAM_ID","ID"), type="inner")
    g <- as.matrix(data[,3:(dim(data)[2]-1)])
    p <- as.matrix(data$IsIncid)
    geno.obj <- SummarizedExperiment(g)
    obj <- rqtClass(phenotype=p, genotype=geno.obj)
  }
  
  #### Association test ####
  tmp.res <- rQTest(obj, method="pca", out.type = "D", cumvar.threshold=75, scale=T)
  if(!is.na(tmp.res@results$p.value$p.Q3)) {
    res <- list(pValue=tmp.res@results$p.value$p.Q3, testStat=tmp.res@results$Qstatistic$Q3, beta=tmp.res@results$beta)
  } else {
    res <- list(pValue=NA, testStat=NA, beta=NA)
  }
  
  #obj.b <- SKAT_Null_Model(p ~ data$cov1, out_type="D")
  #res.skat <- SKAT(g, obj.b)
  #res <- list(pValue=res.skat$p.value, testStat=res.skat$Q, method="optimal.adj")
  
  res
}


runTests <- function(d) { #{{{
  d$tmp_results <- mclapply(d$ifilenames, function(fn) {
    mat <- loadGenTable(d$gene_dir, fn)
    calcrQTest(mat, d)
  })
  
  # Only for debug purposes:
  #d$tmp_results <- lapply(d$ifilenames, function(fn) {
  #                        print(fn)
  #                        mat <- loadGenTable(d$gene_dir, fn)
  #                        calcrQTest(mat, d)
  #                    })
  cat(" Done.\n")
}#}}}


makeTestResults <- function(dd) { #{{{
  perm <- dd$perm
  correctPvalues <- function(pv,dd) {
    pv[is.na(pv)] <- 1
    min_pv <- min(pv[pv>0])
    pv[pv==0] <- min_pv / 50000 # min_pv/perm
    pv[pv>1] <- 1
    pv
  }
  
  cat("-- Preparing test results...")
  rmat <- t(simplify2array(dd$tmp_results))   # rows -- genes, columns -- tests
    
  dd$test_results <- new.env()
  
  tpv <- correctPvalues( unlist(rmat[,"pValue"]), dd )
  
  #if (dd$gcontrol) {
  #  print("Genomic control...")
  #  res <- genomic_control(tpv,method = "median")
  #  tpv <- res$adjusted_p
  #  dd$test_results[[tname]]$lambda <- res$lambda
  #  print("Done genomic control.")
  #}
    
  o <- order(tpv)
  tqv <- qvalue(tpv,lambda=0.1,robust=TRUE)$qvalues
  ignames <- names(tpv)
  gnames <- dd$genenames[ignames]
  tstat <- unlist(rmat[,"testStat"])
  dd$test_results$data <- data.frame( iGene   = ignames[o],
                                               Gene    = gnames[o],
                                               pValue  = tpv[o],
                                               qValue  =tqv[o],
                                               testStat = tstat[o],
                                               row.names=NULL, stringsAsFactors=F )
    
  
  cat(" Done.\n")
  
}#}}}

######################################################################################################################################################### 
###################################################### ANALYSIS #########################################################################################
######################################################################################################################################################### 

####### Prepare list of genes ########
table <- read.xlsx(table_filename, startRow = 3)
genes <- unique(sapply(1:dim(table)[1], function(n) {trim(unlist(strsplit(x = table$Reportedgenesbyauthor[n], split = ','))[1])}))
genes <- genes[which(!is.na(genes))]
write(genes, genes_to_write_filename)

################################################################ ################################################################
###################################################### CHS CARe ################################################################
################################################################ ################################################################

dd.chs <- new.env()
dd.chs$gnamesfile <- gnamesfile
dd.chs$genenames <- loadGeneNames(dd.chs$gnamesfile)
dd.chs$gene_dir <- gene_dir_chs
dd.chs$perm <- 0


#### Prepare pheno ####
columns <- c("fid", "SubjID", "IsIncid")
pheno <- read.csv(paste(pheno.data.dir, pheno.data.dir_chs, sep=''), header = TRUE)
pheno <- pheno[, columns]
colnames(pheno) <- c("FAM_ID", "ID", "IsIncid")
dd.chs$pheno <- pheno
# Save phenotype #
write.table(pheno[,1:2],indiv_chs, row.names = F, col.names=F, quote = F)
write.table(pheno, pheno_bin_chs, row.names = F, col.names=F, quote = F)

print(paste("Total:", dim(pheno)[1], "Cases:", length(which(pheno$IsIncid == 1))))


#### Prepare geno ####
######### Shell scripts ##########
#vpython select_candidate_genes.py -i genes.txt -a  /data/work/iz12/pipeline/ANNOTATION/CHS/ann5Gene.csv -o out_chs.txt
#plink --keep-allele-order --bfile /data/work/iz12/pipeline/ANNOTATION/CHS/CHS --extract out_chs.txt --keep individuals_chs.txt --make-bed --out geno_chs
#plink --keep-allele-order --bfile geno_chs --nonfounders --mind 0.05 --hwe 0.001 --maf 0.05 --geno 0.05 --make-bed --out geno_chs_qc
#plink   --bfile geno_chs_qc --keep-allele-order --recodeA --out geno_chs_qc_recorded
#sed -i -e '1s|_[A-Z]||g' -e '2,$s| NA| .|g' geno_chs_qc_recorded.raw
#mkdir genes_chs
#vpython prepSKAT.py genes_chs  geno_chs_qc_recorded.raw  pheno_bin_chs.txt /data/work/iz12/pipeline/ANNOTATION/CHS/ann5Gene.csv


#### Run tests ####
dd.chs$ifilenames <- loadInputFileNames(dd.chs$gene_dir)
runTests(dd.chs)

# Conduct meta-analysis by summarizing p-values using Fisher's method #
makeTestResults(dd.chs)
colnames(dd.chs$test_results$data) <- c("iGene", "Gene", "chs.pValue", "chs.qValue", "chs.testStat")
#hist(dd.chs$test_results$data$pValue)
dd.chs$test_results$data <- data.frame(dd.chs$test_results$data[with(dd.chs$test_results$data, order(chs.pValue)), ])
head(dd.chs$test_results$data)
write.table(x = dd.chs$test_results$data[which(dd.chs$test_results$data$Gene %in% genes),], 
            file="W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/chs.results.txt", row.names = F, quote=F)


################################################################ ################################################################
###################################################### FHS CARe ################################################################
################################################################ ################################################################
dd.fhs <- new.env()
dd.fhs$gnamesfile <- gnamesfile
dd.fhs$genenames <- loadGeneNames(dd.fhs$gnamesfile)
dd.fhs$gene_dir <- gene_dir_fhs
dd.fhs$perm <- 0


#### Prepare pheno ####
pheno <- read.csv(paste(pheno.data.dir, pheno.data.dir_fhs, sep=''), header = TRUE)
columns <- c("FID", "ID", "IsIncid")
pheno <- pheno[, columns]
colnames(pheno) <- c("FAM_ID", "ID", "IsIncid")
dd.fhs$pheno <- pheno
# Save phenotype #
write.table(pheno[,1:2], indiv_fhs, row.names = F, col.names=F, quote = F)
write.table(pheno, pheno_bin_fhs, row.names = F, col.names=F, quote = F)
print(paste("Total:", dim(pheno)[1], "Cases:", length(which(pheno$IsIncid == 1))))

#### Prepare geno ####
######### Shell scripts ##########
# To execute the following commands below, open the terminal and paste them, then press enter:
#vpython select_candidate_genes.py -i genes.txt -a  /data/work/iz12/pipeline/ANNOTATION/FRAM_CARe/ann5Gene.csv -o out_fhs.txt
#plink --keep-allele-order --bfile /data/work/iz12/pipeline/ANNOTATION/FRAM_CARe/FRAM_CARe --extract out_fhs.txt --keep individuals_fhs.txt --make-bed --out geno_fhs
#plink --keep-allele-order --bfile geno_fhs --nonfounders --mind 0.05 --hwe 0.001 --maf 0.05 --geno 0.05 --make-bed --out geno_fhs_qc
#plink   --bfile geno_fhs_qc --keep-allele-order --recodeA --out geno_fhs_qc_recorded
#sed -i -e '1s|_[A-Z]||g' -e '2,$s| NA| .|g' geno_fhs_qc_recorded.raw
#mkdir genes_fhs
#vpython prepSKAT.py genes_fhs  geno_fhs_qc_recorded.raw  pheno_bin_fhs.txt /data/work/iz12/pipeline/ANNOTATION/FRAM_CARe/ann5Gene.csv
######### End of shell scripts ##########


dd.fhs$ifilenames <- loadInputFileNames(dd.fhs$gene_dir)
runTests(dd.fhs)
makeTestResults(dd.fhs)
colnames(dd.fhs$test_results$data) <- c("iGene", "Gene", "fhs.pValue", "fhs.qValue", "fhs.testStat")
#hist(dd.fhs$test_results$data$pValue)
dd.fhs$test_results$data <- data.frame(dd.fhs$test_results$data[with(dd.fhs$test_results$data, order(fhs.pValue)), ])
print(head(dd.fhs$test_results$data))
write.table(x = dd.fhs$test_results$data[which(dd.fhs$test_results$data$Gene %in% genes),], 
            file="W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/fhs.results.txt", row.names = F, quote=F)

################################################################################################################################
###################################################### LOADFS ################################################################
################################################################################################################################

dd.loadfs <- new.env()
dd.loadfs$gnamesfile <- gnamesfile
dd.loadfs$genenames <- loadGeneNames(dd.loadfs$gnamesfile)
dd.loadfs$gene_dir <- gene_dir_loadfs
dd.loadfs$perm <- 0


#### Prepare pheno ####
pheno <- read.csv(paste(pheno.data.dir, pheno.data.dir_loadfs, sep=''), header = TRUE)
columns <- c("FID", "Subj_NO", "Case_Control")
pheno <- pheno[, columns]
colnames(pheno) <- c("FAM_ID", "ID", "IsIncid")
dd.loadfs$pheno <- pheno
# Save phenotype #
write.table(pheno[,1:2], indiv_loadfs, row.names = F, col.names=F, quote = F)
write.table(pheno, pheno_bin_loadfs, row.names = F, col.names=F, quote = F)
print(paste("Total:", dim(pheno)[1], "Cases:", length(which(pheno$IsIncid == 1))))
#### Prepare geno ####
######### Shell scripts ##########
# To execute the following commands below, open the terminal and paste them, then press enter:
#vpython select_candidate_genes.py -i genes.txt -a  /data/work/iz12/pipeline/ANNOTATION/LOADFS/ann5Gene.csv -o out_loadfs.txt
#plink --keep-allele-order --bfile /data/work/iz12/pipeline/ANNOTATION/LOADFS/LOADFS --extract out_loadfs.txt --keep individuals_loadfs.txt --make-bed --out geno_loadfs
#plink --keep-allele-order --bfile geno_loadfs --nonfounders --mind 0.05 --hwe 0.001 --maf 0.05 --geno 0.05 --make-bed --out geno_loadfs_qc
#plink   --bfile geno_loadfs_qc --keep-allele-order --recodeA --out geno_loadfs_qc_recorded
#sed -i -e '1s|_[A-Z]||g' -e '2,$s| NA| .|g' geno_loadfs_qc_recorded.raw
#mkdir genes_loadfs
#vpython prepSKAT.py genes_loadfs  geno_loadfs_qc_recorded.raw  pheno_bin_loadfs.txt /data/work/iz12/pipeline/ANNOTATION/LOADFS/ann5Gene.csv
######### End of shell scripts ##########


dd.loadfs$ifilenames <- loadInputFileNames(dd.loadfs$gene_dir)
runTests(dd.loadfs)

# Conduct meta-analysis by summarizing p-values using Fisher's method #
makeTestResults(dd.loadfs)
colnames(dd.loadfs$test_results$data) <- c("iGene", "Gene", "loadfs.pValue", "loadfs.qValue", "loadfs.testStat")
#hist(dd.loadfs$test_results$data$pValue)
dd.loadfs$test_results$data <- data.frame(dd.loadfs$test_results$data[with(dd.loadfs$test_results$data, order(loadfs.pValue)), ])
head(dd.loadfs$test_results$data)
write.table(x = dd.loadfs$test_results$data[which(dd.loadfs$test_results$data$Gene %in% genes),], 
            file="W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/loadfs.results.txt", row.names = F, quote=F)


#################################################### #################################################### ######################### 
#################################################### Summarizing p-values #########################################################
#################################################### #################################################### #########################

dat1 <- read.table("W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/fhs.results.txt", header=T)
dat2 <- read.table("W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/chs.results.txt", header=T)
dat3 <- read.table("W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/loadfs.results.txt", header=T)

d1 <- join(dat2, dat1, by=c("iGene"), type="full")
d2 <- join(d1, dat3, by=c("iGene"), type="full")

#d1 <- join(dd.chs$test_results$data, dd.fhs$test_results$data, by=c("iGene"), type="full")
#d2 <- join(d1, dd.loadfs$test_results$data, by=c("iGene"), type="full")


pvv <- d2[,seq(3,dim(d2)[2],by=3)]

tests <- c("wilkinson", "fisher", "minimump", "sump", "sumlog", "logitp")

final.pvv <- data.frame(matrix(nrow=dim(pvv), ncol=0))

for(test in tests) {
  print(test)
  comb.pvv <- sapply(1:dim(pvv)[1], function(n) {
     if(length(which(!is.na(pvv[n,]))) >= 2) {
       if(test == "wilkinson") {
         wilkinsonp(pvv[n, which(!is.na(pvv[n,]))])$p
       } else if(test == "fisher") {
         chi.comb <- sum(-2*log(pvv[n, which(!is.na(pvv[n,]))]))
         df <- 2*length(which(!is.na(pvv[n,])))
         ans <- 1-pchisq(q=chi.comb, df=df)
       } else if(test == "minimump") {
         minimump(pvv[n, which(!is.na(pvv[n,]))])$p
       } else if(test == "sump") {
         sump(pvv[n, which(!is.na(pvv[n,]))])$p
       } else if(test == "sumlog") {
         sumlog(pvv[n, which(!is.na(pvv[n,]))])$p
       } else if(test == "logitp") {
         logitp(pvv[n, which(!is.na(pvv[n,]))])$p
       }
    } else {
       pvv[n, which(!is.na(pvv[n,]))]
    }
    })
  
  final.pvv <- cbind(final.pvv, comb.pvv)
  

}

colnames(final.pvv) <- tests
rownames(final.pvv) <- d2[,1]
final.pvv <- data.frame(cbind(gene=d2[,2], final.pvv))
final.pvv <- data.frame(final.pvv[with(final.pvv, order(wilkinson, fisher, minimump, sump, sumlog, logitp)), ])

pt <- 1
write.table(x = final.pvv[which(final.pvv$gene %in% genes), ], 
            file="W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/meta.p-value.txt", row.names = F, quote=F)

write.xlsx(dd.chs$test_results$data, "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/chs.results.xlsx")
write.xlsx(dd.fhs$test_results$data, "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/fhs.results.xlsx")
write.xlsx(dd.loadfs$test_results$data, "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/loadfs.results.xlsx")
write.xlsx(final.pvv[which(final.pvv$gene %in% genes), ], 
           "W:/data/work/iz12/rqt/tests/app_note/ad.tests/rqt.analysis/meta.p-value.xlsx")
