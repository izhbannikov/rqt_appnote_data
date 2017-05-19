# Tables #
library(openxlsx)

n <- 1000

### With LD ###
data.power.epi.d <- read.table("~/Projects/rqt.dev/tests/power.tests/power.ld.dich/res.table.pow.txt", header=T)
data.power.epi.c <- read.table("~/Projects/rqt.dev/tests/power.tests/power.ld.cont/res.table.pow.txt", header=T)

### Without LD ###
data.power.d <- read.table("~/Projects/rqt.dev/tests/power.tests/power.nold.dich/res.table.pow.txt", header=T)
data.power.c <- read.table("~/Projects/rqt.dev/tests/power.tests/power.nold.cont/res.table.pow.txt", header=T)

## Average VIF ##
#### LD ####
vif.ld <- cbind(data.power.epi.d[,c("vif.rqt", "vif.qt")], 
                data.power.epi.c[,c("vif.rqt", "vif.qt")])
colnames(vif.ld) <- c("VIF.RQT.D", "VIF.QT.D", "VIF.RQT.C", "VIF.QT.C")
write.xlsx(x = vif.ld, file = "~/Projects/rqt.dev/tests/power.tests/vif.ld.xlsx", row.names=T)
#### No LD ####
vif.nold <- cbind(data.power.d[,c("vif.rqt", "vif.qt")], 
                  data.power.c[,c("vif.rqt", "vif.qt")])
colnames(vif.nold) <- c("VIF.RQT.D", "VIF.QT.D", "VIF.RQT.C", "VIF.QT.C")
write.xlsx(x = vif.nold, file = "~/Projects/rqt.dev/tests/power.tests/vif.nold.xlsx", row.names=T)

## Average pooled variance ##
#### LD ####
var.pooled.ld <- cbind(data.power.epi.d[,c("var.pool.rqt", "var.pool.qt")], 
                data.power.epi.c[,c("var.pool.rqt", "var.pool.qt")])
colnames(var.pooled.ld) <- c("Var.pooled.RQT.D", "Var.pooled.QT.D", "Var.pooled.RQT.C", "Var.pooled.QT.C")
write.xlsx(x = var.pooled.ld, file = "~/Projects/rqt.dev/tests/power.tests/var.pooled.ld.xlsx", row.names=T)
#### No LD ####
var.pooled.nold <- cbind(data.power.d[,c("var.pool.rqt", "var.pool.qt")], 
                  data.power.c[,c("var.pool.rqt", "var.pool.qt")])
colnames(var.pooled.nold) <- c("Var.pooled.RQT.D", "Var.pooled.QT.D", "Var.pooled.RQT.C", "Var.pooled.QT.C")
write.xlsx(x = var.pooled.nold, file = "~/Projects/rqt.dev/tests/power.tests/var.pooled.nold.xlsx", row.names=T)

## Average p-values ##
### LD ###
pval.ld <- cbind(data.power.epi.d[,c("prqt1", "prqt2", "prqt3", "pqt1", "pqt2", "pqt3", "pskat", "pskato")], 
                       data.power.epi.c[,c("prqt1", "prqt2", "prqt3", "pqt1", "pqt2", "pqt3", "pskat", "pskato")])
colnames(pval.ld) <- c("p.RQT1.D", "p.RQT2.D", "p.RQT3.D", 
                       "p.QT1.D", "p.QT2.D", "p.QT3.D", 
                       "p.SKAT.D", "p.SKAT-O.D",
                       "p.RQT1.C", "p.RQT2.C", "p.RQT3.C", 
                       "p.QT1.C", "p.QT2.C", "p.QT3.C", 
                       "p.SKAT.C", "p.SKAT-O.C")
write.xlsx(x = pval.ld, file = "~/Projects/rqt.dev/tests/power.tests/pvalue.ld.xlsx", row.names=T)
#### No LD ####
pval.nold <- cbind(data.power.d[,c("prqt1", "prqt2", "prqt3", "pqt1", "pqt2", "pqt3", "pskat", "pskato")], 
                         data.power.c[,c("prqt1", "prqt2", "prqt3", "pqt1", "pqt2", "pqt3", "pskat", "pskato")])
colnames(pval.nold) <- c("p.RQT1.D", "p.RQT2.D", "p.RQT3.D", 
                         "p.QT1.D", "p.QT2.D", "p.QT3.D", 
                         "p.SKAT.D", "p.SKAT-O.D",
                         "p.RQT1.C", "p.RQT2.C", "p.RQT3.C", 
                         "p.QT1.C", "p.QT2.C", "p.QT3.C", 
                         "p.SKAT.C", "p.SKAT-O.C")
write.xlsx(x = pval.nold, file = "~/Projects/rqt.dev/tests/power.tests/pvalue.nold.xlsx", row.names=T)
