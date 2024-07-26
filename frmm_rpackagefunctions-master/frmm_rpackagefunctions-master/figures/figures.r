## Susana, Daphne: add a line for your absolute path to the git folder.
## To run this file locally, uncomment your line and comment the other ones
gitPath <- c('~/git/frmm_rpackagefunctions/') # SHAHIN
#gitPath = c("https://gitlab.com/sconde778/frmm_rpackagefunctions/") ## Susana # CodeFiguresPaper/GitlabRepository/frmm_rpackagefunctions/
#gitPath = ## Daphne

## relative paths
dataPath <- 'data/'
figPath <- 'figures/'

## figure specific parameters
figFullWidth <- 6
figHeight <- 3
figFullHeight = 9
figFullWidth2 <- 10
figHeight2 <- 9.5
figHeight3 <- 10

## Name of our method
ourMethodName = 'FRECL'

library(tidyr)
library(dplyr)
library(reshape)
library(ggplot2)


##############################
## Figures 1 & 2
##############################

### Loading the data and saving as tibble

## N=500
varNames = load(paste0(gitPath, dataPath, "Res_50simulations_48seasonR_K3n500L2iids2_1_FuRMi_Alt1_15-8-20.RData"))
varNames = c(varNames, load(paste0(gitPath, dataPath, "Res_Alt2-3_K3L2n500iid_15-8-20.RData")))
varNames = c(varNames, load(paste0(gitPath, dataPath, "Res_Alt4-5_K3L2n500iid_15-8-20.RData")))
#dim(ari.tpr.our.algorithm) # ari.tpr.our.algorithm[, c(1:3)]

allData = cbind(500, rbind( 
                cbind(ourMethodName, ari.tpr.our.algorithm %>% t),
                cbind('FPCA.oracle', ari.tpr.altern.fpc %>% t),
                cbind('HDDC.Cattell', ari.tpr.alt2.hddc.default %>% t),
                cbind('HDDC.BIC', ari.tpr.alt2.hddc.bic %>% t),
                cbind('FunHDDC.Cattell', ari.tpr.alt4.default.mat %>% t),
                cbind('FunHDDC.BIC', ari.tpr.alt5.bic.mat %>% t)
))

rm(list=varNames)

### N=1000
varNames = load(paste0(gitPath, dataPath, "Res_50simulations_50seasonR_K3n1000L2iids2_1_FuRMi_Alt1_15-8-20.RData")) 
varNames = c(varNames, load(paste0(gitPath, dataPath, "Res_50simulations_50seasonR_K3n1000L2iids2_1_Alt2-3_15-8-20.RData")) )
varNames = c(varNames, load(paste0(gitPath, dataPath, "Res_Alt4-5_K3L2n1000iid_15-8-20.RData")) )


allData = rbind(allData, cbind(1000, rbind( 
                cbind(ourMethodName, ari.tpr.our.algorithm %>% t),
                cbind('FPCA.oracle', ari.tpr.altern.fpc %>% t),
                cbind('HDDC.Cattell', ari.tpr.alt2.hddc.default %>% t),
                cbind('HDDC.BIC', ari.tpr.alt2.hddc.bic %>% t),
                cbind('FunHDDC.Cattell', ari.tpr.alt4.default.mat %>% t),
                cbind('FunHDDC.BIC', ari.tpr.alt5.bic.mat %>% t)
)))

colnames(allData)[1:2] = c('N', 'Method')
allData = as_tibble(allData)
allData = mutate(allData, SimNumber=seq(nrow(allData)))
library(readr)
allData = type_convert(allData)
allData$Method = factor(allData$Method, levels=c(ourMethodName, "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC"))
allData$N = as.factor(allData$N)

#### Settings for getting base-r type plots
theme_set(theme_bw())
theme_update(text = element_text(size=12),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
strip.background = element_blank()
)

#### Saving the boxplots
pdf(file=paste0(gitPath, figPath, 'fig1_ARI.pdf'), width=figFullWidth, height=figHeight)
allData %>% ggplot(aes(x=Method, y=`%ARI`, fill=N)) +
    geom_boxplot(lwd=.1, outlier.size=.1) +
    coord_cartesian(ylim=c(0,100)) + # limit y-axis without removing data points
    labs(title='Adjusted Rand Index (%)', y='', x='', fill='Sample size') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file=paste0(gitPath, figPath, 'fig1_RI.pdf'), width=figFullWidth, height=figHeight)
allData %>% ggplot(aes(x=Method, y=`%RI`, fill=N)) +
    geom_boxplot(lwd=.1, outlier.size=.1) +
    coord_cartesian(ylim=c(0,100)) + # limit y-axis without removing data points
    labs(title='Rand Index (%)', y='', x='', fill='Sample size') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
dev.off()


pdf(file=paste0(gitPath, figPath, 'fig2_TPR.pdf'), width=figFullWidth, height=figHeight)
allData %>% ggplot(aes(x=Method, y=`%TPR`, fill=N)) +
    geom_boxplot(lwd=.1, outlier.size=.1) +
    coord_cartesian(ylim=c(0,100)) + # limit y-axis without removing data points
    labs(title='True Positive Rate (%)', y='', x='', fill='Sample size') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file=paste0(gitPath, figPath, 'fig2_TNR.pdf'), width=figFullWidth, height=figHeight)
allData %>% ggplot(aes(x=Method, y=`%TNR`, fill=N)) +
    geom_boxplot(lwd=.1, outlier.size=.1) +
    coord_cartesian(ylim=c(0,100)) + # limit y-axis without removing data points
    labs(title='True Negative Rate (%)', y='', x='', fill='Sample size') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
dev.off()


##############################
## Figure 3
##############################


load(paste0(gitPath, dataPath, "Data_Figure1b_ARI50_AR1vsIID_n500-1000AllMeth.RData") ) # "ari50.ar1.n1000" "ari50.ar1.n500"  "ari50.iid.n1000" "ari50.iid.n500"  "pathtemp3

n1 <- 500
n2 <- 1000
results.ar1 <- rbind(colMeans(ari50.ar1.n500), colMeans(ari50.ar1.n1000))
rownames(results.ar1) <- paste0("n", c(n1, n2))
colnames(results.ar1) = c(ourMethodName, "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC") ####### NEEDS CHECKING by SUSANA!!
results.iid <- rbind(colMeans(ari50.iid.n500), colMeans(ari50.iid.n1000))
rownames(results.iid) <- paste0("n", c(n1, n2))
colnames(results.iid) = c(ourMethodName, "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC") ####### NEEDS CHECKING by SUSANA!!

## TODO: change the name and order of the columns to correspond to name you want to be displayed and their desired order
# all.equal(colnames(results.ar1), colnames(results.iid)) ## just a sanity check

## set the colors and line types
cols <- c("darkblue", "chocolate", "chartreuse4", "brown1", "orange", "mediumorchid") # like violet or perhaps magenta
pchs <- 1:6
cex = .7
lty.iid <- 1
lty.ar1 <- 2

n.mat <- matrix(c(n1, n2), nrow=2, ncol=ncol(results.iid))
set.seed(1) # for reproducing the jittering
pdf(paste0(gitPath, figPath, 'fig3.pdf'), width=figFullWidth, height=figHeight, pointsize=10)
par(mar=c(5.1, 4.1, 1.1, 11.1), xpd=TRUE)
matplot(jitter(n.mat, factor=.1), results.iid, type='b', lty=lty.iid, ylim=c(0,100), col=cols, pch=pchs, xlab ="n", ylab = "Adjusted Rand Index (%)", xaxt='n', cex=cex)
axis(side=1, at=c(n1,n2), labels=c(n1,n2))
matlines(jitter(n.mat, factor=.1), results.ar1,  lty=lty.ar1, col=cols, pch=pchs, type='b', cex=cex)
legend('topright', inset=c(-.5,0), legend=c(colnames(results.iid), "iid errors", 'AR(1) errors'), pch=c(pchs, NA,NA),  col=c(cols, 1,1), lty=c(rep(NA,6), lty.iid, lty.ar1))
dev.off()


##############################
## Figure 4
##############################

### This is the data I found in the old file
### it is NOT REPRODUCIBLE and the data NEEDS TO BE LOADED DIRECTLY FROM A .RData file - *** THIS IS ALREADY DONE; PLEASE SEE UPDATED CODE BELOW. ***
### ARIperIterationInFRMM_... files (4): see box.com, FuRMi.
###
filenames.vec <- c("SolsBriefARIandFinalPartition_FuRMi_Alt1-5_K3n500L2AR1phi0_5s2_0_1_15-9-20_58seasonR.RData", "SolsBriefARIandFinalPartition_FuRMi_Alt1-5_K6n500L2AR1phi0_5s2_0_1_15-9-20_69seasonR.RData", 
 "SolsBriefARIandFinalPartition_FuRMi_Alt1-5_K9n500L2AR1phi0_5s2_0_1_15-9-20_67seasonR.RData",
 "SolsBriefARIandFinalPartition_FuRMi_Alt1-5_K12n500L2AR1phi0_5s2_0_1_15-9-20_65seasonR.RData") # K = 12, n=1000, L2 distance, AR(1), r=1st simulated data set. #
filenames1000.vec <- c("SolsBriefARIandFinalPartition_FuRMi_Alt1-5_K3n1000L2AR1phi0_5s2_0_1_15-9-20_60seasonR.RData", "SolsBriefARIandFinalPartition_FuRMi_Alt1-5_K6n1000L2AR1phi0_5s2_0_1_15-9-20_70seasonR.RData", 
 "SolsBriefARIandFinalPartition_FuRMi_Alt1-5_K9n1000L2AR1phi0_5s2_0_1_15-9-20_68seasonR.RData", 
 "SolsBriefARIandFinalPartition_FuRMi_Alt1-5_K12n1000L2AR1phi0_5s2_0_1_15-9-20_66seasonR.RData") 
data500 <- matrix(NA, nrow=4, ncol=6)
rownames(data500) <- paste0("K", c(3, 6, 9, 12))
method.lab <- c("FRECL", "FPC.oracle",  "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
colnames(data500) <- method.lab
data1000 <- data500

for(i in 1:4) # n=500 #
{
 load(filenames.vec[i])
 data500[i, ] <- ari.tpr.furmi.alt1.5["%ARI", ]
}

for(i in 1:4) # n=1000 #
{
 load(filenames1000.vec[i])
 data1000[i, ] <- ari.tpr.furmi.alt1.5["%ARI", ] 
}

data500
data1000
# ######################################################################## #
#        FRECL FPC.oracle HDDC.Cattell HDDC.BIC FunHDDC.Cattell FunHDDC.BIC
# K3  87.77833  17.052288     9.152097 9.712343        8.637560   16.811067
# K6  69.52058   8.174384     5.145261 5.723684        4.208810    6.061212
# K9  30.94162  11.303351     1.802568 7.409953        3.170575    7.024331
# K12 16.16378   8.160929     5.881035 5.428198        3.800409    3.539659

#        FRECL FPC.oracle HDDC.Cattell  HDDC.BIC FunHDDC.Cattell FunHDDC.BIC
# K3  93.49990  18.656469     4.900801 12.864430       13.847023   20.663555
# K6  84.07439   7.265435     2.919836  3.621690        7.167989    5.881535
# K9  77.04684   9.482990     1.400790  4.989690        4.788361    6.939552
# K12 36.61812   9.051169     5.505215  6.338536        3.528132    5.891347
# ######################################################################### #

Ks <- seq(3, 12, by=3)

## graphical parameters
cols <- c("darkblue", "chocolate", "chartreuse4", "brown1", "orange", "mediumorchid") 
pchs <- 1:6
lty500=1
lty1000=2
cex = .7

set.seed(1) # for reproducing the jittering
pdf(paste0(gitPath, figPath, 'fig4.pdf'), width=figFullWidth, height=figHeight, pointsize=10)
par(mar=c(5.1, 4.1, 1.1, 11.1), xpd=TRUE)
matplot(jitter(Ks, factor=.1), data500, type='b', lty=lty500, ylim=c(0,100), col=cols, pch=pchs, xlab ="K", ylab = "Adjusted Rand Index (%)", xaxt='n', cex=cex)
matlines(jitter(Ks, factor=.1), data1000, type='b', lty=lty1000, ylim=c(0,100), col=cols, pch=pchs, xaxt='n', cex=cex)
axis(side=1, at=Ks, labels=Ks)
legend('topright', inset=c(-.5,0), legend=c(colnames(data500), "N=500", 'N=1000'), pch=c(pchs, NA,NA),  col=c(cols, 1,1), lty=c(rep(NA,6), lty500, lty1000))
dev.off()



##############################
## Figure 5
##############################
# ARI versus iteration number, K=3, 12, n=500, 1000.

setwd(path1)  
load(gitPath, dataPath, "data_Figure5_ARIvsIteration_K3-12n500-1000.RData") # ari.mat.list.K3n500, ari.mat.list.K3n1000, ari.mat.list.K12n500, ari.mat.list.K12n1000 #
load(gitPath, dataPath, "data_Figure5_2_FinalARI_FuRMi_78seasonR.RData")

ari.mat.list <- ari.mat.list.K3n500
mtemp <- 51
min.y <- -2
max.y <- 100
xtemp <- c(2:mtemp)
colvec <- rep(NA, M.run)
for(l in 1:M.run) 
 colvec[l] <- rgb(0, 0, 255, max = 255, alpha = 50, names = paste0("blue", l))
cextemp <- 0.5
lwdtemp <- 2

pchtemp <- 19

l <- 1
ytemp <- ari.mat.list[[l]]["%ARI", c(2:mtemp)]
pdf(paste0(gitPath, figPath, "Fig5aARI_versus_IterationNo_K3n500L2AR1.pdf"), width=10, height=9.5)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) 
plot(xtemp, ytemp, type='o', xlim=c(1, mtemp), xlab="", ylab="Adjusted Rand Index (*100)", pch=19, cex=cextemp,
 lwd=lwdtemp, col=colvec[l], ylim=c(min.y, max.y))
for(l in 2:M.run)
{
 points(xtemp, ari.mat.list[[l]]["%ARI", c(2:mtemp)], pch=19, col=colvec[l], cex=cextemp, lwd=lwdtemp, type='o')
} # StudyARIperIteration_78seasonR
abline(h=final.ari.vec["K3n500"], lty=1, lwd=2, col=2)
dev.off()

ari.mat.list <- ari.mat.list.K3n1000
l <- 1
ytemp <- ari.mat.list[[l]]["%ARI", c(2:mtemp)]
pdf(paste0(gitPath, figPath, "Fig5bARI_versus_IterationNo_K3n1000L2AR1.pdf"), width=10, height=9.5)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) 
plot(xtemp, ytemp, type='o', xlim=c(1, mtemp), xlab="", ylab="Adjusted Rand Index (*100)", pch=19, cex=cextemp,
 lwd=lwdtemp, col=colvec[l], ylim=c(min.y, max.y))
for(l in 2:M.run)
 points(xtemp, ari.mat.list[[l]]["%ARI", c(2:mtemp)], pch=19, col=colvec[l], cex=cextemp, lwd=lwdtemp, type='o')
abline(h=final.ari.vec["K3n1000"], lty=1, lwd=2, col=2)
dev.off()

ari.mat.list <- ari.mat.list.K12n500
mtemp <- 70
min.y <- 0
max.y <- 42
xtemp <- c(2:mtemp)

l <- 1
ytemp <- ari.mat.list[[l]]["%ARI", c(2:mtemp)]
pdf(paste0(gitPath, figPath,  "Fig5cARI_versus_IterationNo_K12n500L2AR1.pdf"), width=10, height=9.5)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) 
plot(xtemp, ytemp, type='o', xlim=c(1, mtemp), xlab="", ylab="Adjusted Rand Index (*100)", pch=19, cex=cextemp,
 lwd=lwdtemp, col=colvec[l], ylim=c(min.y, max.y))
for(l in 2:M.run)
 points(xtemp, ari.mat.list[[l]]["%ARI", c(2:mtemp)], pch=19, col=colvec[l], cex=cextemp, lwd=lwdtemp, type='o')
abline(h=final.ari.vec["K12n500"], lty=1, lwd=2, col=2)
dev.off()

ari.mat.list <- ari.mat.list.K12n1000
l <- 1
ytemp <- ari.mat.list[[l]]["%ARI", c(2:mtemp)]
pdf(paste0(gitPath, figPath, "Fig5dARI_versus_IterationNo_K12n1000L2AR1.pdf"), width=10, height=9.5)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) 
plot(xtemp, ytemp, type='o', xlim=c(1, mtemp), xlab="", ylab="Adjusted Rand Index (*100)", pch=19, cex=cextemp,
 lwd=lwdtemp, col=colvec[l], ylim=c(min.y, max.y))
for(l in 2:M.run)
 points(xtemp, ari.mat.list[[l]]["%ARI", c(2:mtemp)], pch=19, col=colvec[l], cex=cextemp, lwd=lwdtemp, type='o')
abline(h=final.ari.vec["K12n1000"], lty=1, lwd=2, col=2)
dev.off()

# Final ARI: consensus matrix, K-means - construct data_Figure5_2_FinalARI_FuRMi_78seasonR.RData #
load(paste0(gitPath, dataPath, "ARIperIterationInFRMM_78seasonR_K3L2n500AR1.RData"))
final.ari.vec <- rep(NA, 4)
names(final.ari.vec) <- c("K3n500", "K3n1000", "K12n500", "K12n1000")
binary.mat.true.partition <- binary.matrix.truth.function(K1=K, ngen1=ngen, true.tildetildeP1=true.tildetildeP) 
genelab.ordered.true.partition <- rownames(binary.mat.true.partition) 
summary(achieved.convergence.run.vec) 
system.time(
prop.pairs <- prop.correct.pairs.lth.partition.and.truth.function(M1=M.run,
 K6=K, ngen6=ngen, final.clme.per.run.mat1=final.clme.per.run.mat, genelab.m6=genelab.m,
 true.tildetildeP6=true.tildetildeP, m.gen.clu6=as.numeric(table(true.tildetildeP)),
 genelab.ordered.true.partition6=genelab.ordered.true.partition,
 achieved.convergence.run.vec1=parent.frame()$achieved.convergence.run.vec) 
) 
consensus.mat <- prop.pairs$B.consensus.mat 
set.seed(21)
clme.kmeans <- kmeans(consensus.mat, centers=K, iter.max=100, nstart=50) 
subject.lab.ordered.true.partition <- genelab.ordered.true.partition
Q <- 1 
part.frmm.mat <- matrix(clme.kmeans$cluster) 
rownames(part.frmm.mat) <- genelab.m 
colnames(part.frmm.mat) <- c("FuRMi")
system.time(
res.list <- prop.correct.pairs.lth.partition.and.truth.function(M1=Q, K6=K, ngen6=ngen,
 final.clme.per.run.mat1=part.frmm.mat, genelab.m6=genelab.m, true.tildetildeP6=true.tildetildeP, 
 m.gen.clu6=as.numeric(table(true.tildetildeP)),
 genelab.ordered.true.partition=subject.lab.ordered.true.partition,
 achieved.convergence.run.vec1=1) 
) 
ari.tpr.furmi.mat <- res.list$prop.correct.ones.zeros.acc.FP.FN.mat.ari[c("ARI", "Acc", "TPR", "TNR"), , drop=FALSE]*100
rownames(ari.tpr.furmi.mat) <- c("%ARI", "%RI", "%TPR", "%TNR")
colnames(ari.tpr.furmi.mat)
ari.tpr.furmi.mat
final.ari.vec["K3n500"] <- ari.tpr.furmi.mat["%ARI", ]
#            l1
# %ARI 87.77833
# %RI  94.57876
# %TPR 91.84634
# %TNR 95.93678
# ######################################################################### #
# ############################################################################################# #
load(paste0(gitPath, dataPath, "ARIperIterationInFRMM_78seasonR_K3L2n1000AR1.RData"))
binary.mat.true.partition <- binary.matrix.truth.function(K1=K, ngen1=ngen, true.tildetildeP1=true.tildetildeP) 
genelab.ordered.true.partition <- rownames(binary.mat.true.partition) #"Ahg493925" ...
summary(achieved.convergence.run.vec) # 6-40 no NAs #
system.time(
prop.pairs <- prop.correct.pairs.lth.partition.and.truth.function(M1=M.run, K6=K, ngen6=ngen, final.clme.per.run.mat1=final.clme.per.run.mat, genelab.m6=genelab.m,
 true.tildetildeP6=true.tildetildeP, m.gen.clu6=as.numeric(table(true.tildetildeP)),
 genelab.ordered.true.partition6=genelab.ordered.true.partition,
 achieved.convergence.run.vec1=parent.frame()$achieved.convergence.run.vec) 
) # Next: K-means with the consensus matrix plus final ARI
consensus.mat <- prop.pairs$B.consensus.mat 
set.seed(21)
clme.kmeans <- kmeans(consensus.mat, centers=K, iter.max=100, nstart=50)
clme.kmeans$cluster[1:5] 
subject.lab.ordered.true.partition <- genelab.ordered.true.partition
Q <- 1 
part.frmm.mat <- matrix(clme.kmeans$cluster) 
rownames(part.frmm.mat) <- genelab.m 
colnames(part.frmm.mat) <- c("FuRMi")
system.time(
res.list <- prop.correct.pairs.lth.partition.and.truth.function(M1=Q, K6=K, ngen6=ngen, final.clme.per.run.mat1=part.frmm.mat, genelab.m6=genelab.m, true.tildetildeP6=true.tildetildeP, 
 m.gen.clu6=as.numeric(table(true.tildetildeP)),
 genelab.ordered.true.partition=subject.lab.ordered.true.partition,
 achieved.convergence.run.vec1=1) 
) 
ari.tpr.furmi.mat <- res.list$prop.correct.ones.zeros.acc.FP.FN.mat.ari[c("ARI", "Acc", "TPR", "TNR"), , drop=FALSE]*100
rownames(ari.tpr.furmi.mat) <- c("%ARI", "%RI", "%TPR", "%TNR")
colnames(ari.tpr.furmi.mat)
ari.tpr.furmi.mat
final.ari.vec["K3n1000"] <- ari.tpr.furmi.mat["%ARI", ]
#            l1
# %ARI 93.49990
# %RI  97.11371
# %TPR 95.67905
# %TNR 97.82890
# ##
load(gitPath, dataPath, "ARIperIterationInFRMM_78seasonR_K12L2n500AR1.RData")
binary.mat.true.partition <- binary.matrix.truth.function(K1=K, ngen1=ngen, true.tildetildeP1=true.tildetildeP) 
genelab.ordered.true.partition <- rownames(binary.mat.true.partition) #"Ahg493925" ...
summary(achieved.convergence.run.vec) # 9-20 no NAs #
system.time(
prop.pairs <- prop.correct.pairs.lth.partition.and.truth.function(M1=M.run, K6=K, ngen6=ngen, final.clme.per.run.mat1=final.clme.per.run.mat, genelab.m6=genelab.m,
 true.tildetildeP6=true.tildetildeP, m.gen.clu6=as.numeric(table(true.tildetildeP)),
 genelab.ordered.true.partition6=genelab.ordered.true.partition,
 achieved.convergence.run.vec1=parent.frame()$achieved.convergence.run.vec) 
) 
consensus.mat <- prop.pairs$B.consensus.mat 
set.seed(22)
clme.kmeans <- kmeans(consensus.mat, centers=K, iter.max=100, nstart=50)
subject.lab.ordered.true.partition <- genelab.ordered.true.partition
Q <- 1 
part.frmm.mat <- matrix(clme.kmeans$cluster) 
rownames(part.frmm.mat) <- genelab.m 
colnames(part.frmm.mat) <- c("FuRMi")
system.time(
res.list <- prop.correct.pairs.lth.partition.and.truth.function(M1=Q, K6=K, ngen6=ngen,
 final.clme.per.run.mat1=part.frmm.mat, genelab.m6=genelab.m, true.tildetildeP6=true.tildetildeP, 
 m.gen.clu6=as.numeric(table(true.tildetildeP)),
 genelab.ordered.true.partition=subject.lab.ordered.true.partition,
 achieved.convergence.run.vec1=1) 
)
ari.tpr.furmi.mat <- res.list$prop.correct.ones.zeros.acc.FP.FN.mat.ari[c("ARI", "Acc", "TPR", "TNR"), , drop=FALSE]*100
rownames(ari.tpr.furmi.mat) <- c("%ARI", "%RI", "%TPR", "%TNR")
colnames(ari.tpr.furmi.mat)
ari.tpr.furmi.mat
final.ari.vec["K12n500"] <- ari.tpr.furmi.mat["%ARI", ]
final.ari.vec
#            l1
# %ARI 16.38456
# %RI  86.77756
# %TPR 25.04917
# %TNR 92.25533
# ##
load(gitPath, dataPath, "ARIperIterationInFRMM_78seasonR_K12L2n1000AR1.RData")
binary.mat.true.partition <- binary.matrix.truth.function(K1=K, ngen1=ngen, true.tildetildeP1=true.tildetildeP) 
genelab.ordered.true.partition <- rownames(binary.mat.true.partition) 
summary(achieved.convergence.run.vec)
system.time(
prop.pairs <- prop.correct.pairs.lth.partition.and.truth.function(M1=M.run, K6=K, ngen6=ngen, final.clme.per.run.mat1=final.clme.per.run.mat, genelab.m6=genelab.m,
 true.tildetildeP6=true.tildetildeP, m.gen.clu6=as.numeric(table(true.tildetildeP)),
 genelab.ordered.true.partition6=genelab.ordered.true.partition,
 achieved.convergence.run.vec1=parent.frame()$achieved.convergence.run.vec) 
) 
consensus.mat <- prop.pairs$B.consensus.mat 
set.seed(23)
clme.kmeans <- kmeans(consensus.mat, centers=K, iter.max=100, nstart=50)
clme.kmeans$cluster[1:5] 
subject.lab.ordered.true.partition <- genelab.ordered.true.partition
Q <- 1 
part.frmm.mat <- matrix(clme.kmeans$cluster) 
rownames(part.frmm.mat) <- genelab.m 
colnames(part.frmm.mat) <- c("FuRMi")
system.time(
res.list <- prop.correct.pairs.lth.partition.and.truth.function(M1=Q, K6=K, ngen6=ngen,
 final.clme.per.run.mat1=part.frmm.mat, genelab.m6=genelab.m, true.tildetildeP6=true.tildetildeP, 
 m.gen.clu6=as.numeric(table(true.tildetildeP)),
 genelab.ordered.true.partition=subject.lab.ordered.true.partition,
 achieved.convergence.run.vec1=1) 
) 
ari.tpr.furmi.mat <- res.list$prop.correct.ones.zeros.acc.FP.FN.mat.ari[c("ARI", "Acc", "TPR", "TNR"), , drop=FALSE]*100
rownames(ari.tpr.furmi.mat) <- c("%ARI", "%RI", "%TPR", "%TNR")
colnames(ari.tpr.furmi.mat)
ari.tpr.furmi.mat
final.ari.vec["K12n1000"] <- ari.tpr.furmi.mat["%ARI", ]
final.ari.vec 
save(final.ari.vec, file=paste0(gitPath, dataPath,
"data_Figure5_2_FinalARI_FuRMi_78seasonR.RData"))

87.77#            l187.77 
# %ARI 36.73303
# %RI  90.14454
# %TPR 43.49738
# %TNR 94.33446
# ################################## #
#   K3n500  K3n1000  K12n500 K12n1000 
# 87.77833 93.49990 16.38456 36.73303 
# ################################## #

##############################
## Figure 6
##############################
# DISTRIBUTION OF ARI AGAINST THE NUMBER OF RUNS - FuRMi, 10-100 RUNS; K=12; N=1000, L2, AR(1).

load(paste0(gitPath, dataPath, "data_Figure4_DistributionARIagainstNoRunsFRMM_All10to100Plus5-15_82seasonR_21-5-20.RData") )


## Susana, I've set the figure width in a way that we can include the figures in the .tex without rescaling them; please don't make them wider.
pdf(file=paste0(gitPath, figPath, 'fig6.pdf'), width=figFullWidth, height=figHeight, pointsize=10)
par(mgp=c(3.0, 1, 0), mar=c(4, 4, 1, 1)+0.1) 
X = perc.ari.finalfrmmpartition.mat.total12b
boxplot(X, at=as.integer(dimnames(X)[[2]]), xlab="Number of runs", ylab="Adjusted Rand Index (%)", pars = list(boxwex=4))  # boxplots_ARIvsNoRunsFRMM_Every10To100Plus5-15Runs_82seasonR_21-5-20.pdf
dev.off()



##############################
## Figure 7 (A)
##############################
# MSE versus no. clusters, real seasonal data set.

load(paste0(gitPath, dataPath, "MeanSquaredErrorTotal_K2-18_79AseasonR_15-7-20.RData") )
K <- 18
x <- c(2:K) 
temp <- mse.vec
min.temp <- min(temp)
max.temp <- max(temp)
cextemp2 <- 1.4
lwdtemp2 <- 3

pdf(paste0(gitPath,  figPath, "lineplot_MSE_FRMMrealSeasonSet_K2-", K, "_79AseasonR_15-7-20.pdf"), width=figFullWidth2, height=figHeight3)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) 
plot(x, temp, type='o', pch=pchtemp, cex=0.6, lwd=lwdtemp2,
 ylim=c(min.temp, max.temp), xaxt="none", cex.axis=cextemp2, xlab="K", 
 ylab="Mean Squared Error", col="darkgreen")
axis(1, x) # lineplot_MSE_FRMMrealSeasonSet_K2-18_79AseasonR_15-7-20.pdf
dev.off()

# Make of MeanSquaredErrorTotal_K2-18_79AseasonR_15-7-20.RData

# 7a) Mean squared error (adding up the errors of the FDboost fittings in all the clusters) versus the number
# of clusters, real seasonal data set, 23 time points.

rm(list = ls(all = TRUE))
library(FDboost) 
min.K <- 2 
max.K <- 18
mse.K.vec <- rep(NA, max.K-min.K+1)
names(mse.K.vec) <- paste0("K", min.K:max.K) 
#
source(paste0(gitPath, "R/functions1_3-6-20.R"))
system.time(
for(K in min.K:max.K) 
{
 # ################################################################################### #
 # ################################################################################### #
 # ################################################################################### #
 # MEAN SQUARED ERROR FOR THE FRMM FINAL PARTITION WITH K CLUSTERS: MeanSquaredError_K...RData #
 # Load Objects3..., FinalPartitFRMMandTotalSumResid100runs_K, 
 # for the kth cluster in part.frmm.Ktemp.mat, fit an FDboost model,
 # get the predicted values of each observation (gene) for the model fitted to that cluster:
 # predicted.temp.testing.k, estimate the squared residuals, in squared.residuals.l2.mat;
 # finally, find the mean squared error of the partition, using these residuals. Save
 # it in MeanSquaredError_K.._79seasonR_20-4-20.RData
 # ################################################################################### #
 # ################################################################################### #
 # ################################################################################### #
 cat("# ############ # K = ", K, "# ############ #\n")       
 load(paste0(gitpath, dataPath, "Objects3_K", K, "_100RunsRealSeasonalDataSet_79A_24-6-20.RData") )
 load(gitPath, dataPath, paste0("FinalPartitFRMMandTotalSumResid100runs_K", 
  K, "_79AseasonR_24-6-20.RData") ) 
 cat("# ############ # K = ", K, "# ############ #\n")
 m.train <- as.vector(table(part.frmm.Ktemp.mat) ) 
 rownames(part.frmm.Ktemp.mat)[part.frmm.Ktemp.mat==1][1:5] 
 genelab.clu.j.list <- list() # Gene labels for genes in the k-th cluster:
 for(k in 1:K)
 {
  genelab.clu.j.list[[k]] <- rownames(part.frmm.Ktemp.mat)[part.frmm.Ktemp.mat==k ]
 }
 # mean squared error: observations of the same cluster in each testing set #
 squared.residuals.l2.mat <- matrix(NA, nrow=ngen, ncol=nT)
 rownames(squared.residuals.l2.mat) <- genelab.m 
 colnames(squared.residuals.l2.mat) <- timep 
 seed2 <- 5 # before call to FDboost #
 for(k in 1:K) # n.knots=6 # FOR EACH OF THE CLUSTERS k IN THE PARTITION FOR K=K #
 {
  mtemp <- m.train[k]
  dataset.k.train.x1x2x3.loess.list <- list()
  temp <- y.summer.m.transf.loess[genelab.clu.j.list[[k]], ] # mtemp x 24 time points
  dataset.k.train.x1x2x3.loess.list$Ytemp <- temp
  temp1 <- x1.spring.m.transf.loess[genelab.clu.j.list[[k]], ]
  dataset.k.train.x1x2x3.loess.list$x1temp.centred <- scale(temp1, center=TRUE, scale=FALSE)
  temp2 <- x2.autumn.m.transf.loess[genelab.clu.j.list[[k]], ]
  dataset.k.train.x1x2x3.loess.list$x2temp.centred <- scale(temp2, center=TRUE, scale=FALSE)
  temp3 <- x3.winter.m.transf.loess[genelab.clu.j.list[[k]], ]
  dataset.k.train.x1x2x3.loess.list$x3temp.centred <- scale(temp3, center=TRUE, scale=FALSE)
  dataset.k.train.x1x2x3.loess.list$timep <- timep
  dataset.k.train.x1x2x3.loess.list$s.timep <- timep
  dataset.k.train.x1x2x3.loess.list$df.temp5 <- df.temp5
  dataset.k.train.x1x2x3.loess.list$genelabs.k.clu <- genelab.clu.j.list[[k]]
  set.seed(seed2) # E(Yi|Xi)=beta0(t)+int beta2(s, t)*x2(s)*ds+int beta3(s, t)*x3(s)*ds # usually seed2=5 say
  a1 <- try(fof.functional.temp <- FDboost( Ytemp ~ 1 + bsignal(x=x1temp.centred, s=s.timep, knots=n.knots, df = df.temp5) +
   bsignal(x=x2temp.centred, s=s.timep, knots=n.knots, df = df.temp5) + bsignal(x=x3temp.centred, s=s.timep, knots=n.knots, df = df.temp5),
   timeformula = ~ bbs(timep, knots=n.knots, df = df.temp5), data = dataset.k.train.x1x2x3.loess.list, control = boost_control(mstop = 300)) )
  if(class(a1)[1]=="try-error")
  { # https://stat.ethz.ch/R-manual/R-devel/library/base/html/try.html
   cat("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n")
   cat(" k = ", k, "th cluster or model. FDboost command with\n",  " no.knots = ",
    n.knots, " reports an error. It is \n", a1[1])
   cat("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n")
   n.knots2 <- 7
   b <- try(fof.functional.temp <- FDboost( Ytemp ~ 1 + bsignal(x=x1temp.centred, s=s.timep, knots=n.knots2, df = df.temp5) +
    bsignal(x=x2temp.centred, s=s.timep, knots=n.knots2, df = df.temp5) +
    bsignal(x=x3temp.centred, s=s.timep, knots=n.knots2, df = df.temp5),
    timeformula = ~ bbs(timep, knots=n.knots2, df = df.temp5),
    data = dataset.k.train.x1x2x3.loess.list, control = boost_control(mstop = 300)) )
    if(class(b)[1]=="try-error")
     stop("With 7 knots, there is an error too.\n")
   } # end if
    # ########################################### #
   # # ################################################################################ #
   predicted.temp.testing.k <- predict(fof.functional.temp, newdata=dataset.k.train.x1x2x3.loess.list) # 2094 24
    ntemp <- prod(dim(temp)) 
   if(prod(dim(squared.residuals.l2.mat[genelab.clu.j.list[[k]], ]))!=prod(dim((temp-predicted.temp.testing.k)^2)))
    stop("Dimension of squared.residuals.l2.mat differs from dimension of (Y_k-^Y_k)^2.\n")
    squared.residuals.l2.mat[genelab.clu.j.list[[k]], ] <- (temp-predicted.temp.testing.k)^2 # SQUARE OF EACH RESIDUAL ################### # 
   # ###################################################### #
  } # End for(k in 1:K) #
  explanation <- paste0(c("mse.K: calculated using equal training and testing sets, and these are formed by the
   observations (=genes) in each cluster, and fitting each FDboost model to those. FRMM partition found with K =",
   K, "\n") )
  print(mse.K <- sum(squared.residuals.l2.mat)/(ngen*nT)) 
  mse.K.vec[paste0("K", K)] <- mse.K
  save(mse.K, K, explanation, itemnumber, datetemp, ngen, nT, file=paste0(gitPath, dataPath, "MeanSquaredError_K", K, "_79AseasonR_24-6-20.RData"))
  # ###################################################### #
 } # End of for(K in 2:18) # ############################# # 
)

mse.vec <- mse.K.vec 
mse.vec 
#        K2        K3        K4        K5        K6        K7        K8        K9       K10       K11       K12       K13       K14       K15       K16       K17       K18 
# 0.9467330 0.8322436 0.7631034 0.7313826 0.6860181 0.6663203 0.6452519 0.6183390 0.6014889 0.6029251 0.5987175 0.5783899 0.5678991 0.5658196 0.5497500 0.5572379 0.5530234
save(mse.vec, file=paste0(gitPath, dataPath, "MeanSquaredErrorTotal_K2-", K, "_79AseasonR_15-7-20.RData")) # A_SeasonalSet12TimePoints


##############################
## Figure 7 (B)
##############################
# MSE versus no. clusters, real seasonal data set, 24 time points over 2 days.

rm(list = ls(all = TRUE))
load(paste0(gitPath, dataPath, "MeanSquaredErrorTotal_K2-18_79AseasonR_15-7-20.RData") )

# Plot of the MSEs by number of clusters K.
x <- c(2:K) 
print(temp <- mse.K.vec)
min.temp <- min(temp)
max.temp <- max(temp)
cextemp2 <- 1.4
lwdtemp2 <- 3
figFullWidth2 <- 10
figHeight3 <- 10
pchtemp <- 19
cextemp2 <- 1.4
lwdtemp2 <- 3
pdf(paste0(gitPath, figPath, "lineplot_MSE_FRMMrealSeasonSet_K2-15_79BseasonR_12-11-20.pdf"), width=figFullWidth2, height=figHeight3)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) 
plot(x, temp, type='o', pch=pchtemp, cex=0.6, lwd=lwdtemp2,
     ylim=c(min.temp, max.temp), cex.axis=cextemp2, xaxt="none", xlab="K", ylab="Mean Squared Error", col="darkgreen")
axis(1, x) # lineplot_MSE_FRMMrealSeasonSet_K2-15_79BseasonR_12-11-20.pdf
dev.off() # BUT DIMENSIONS: 'ROUGH ESTIMATION' AS THE DIMENSIONS OF THE PDF DO NOT WORK WELL.

# To produce RData files above:
# 9) Observed mean squared error of the models found by FRMM against the number of clusters. #
# Real seasonal data set, 24 time points "spread" over 2 days per season. #

rm(list = ls(all = TRUE)) # from 79B) #
library(FDboost)
min.K <- 2
max.K <- 15
mse.K.vec <- rep(NA, max.K-min.K+1)
names(mse.K.vec) <- paste0("K", min.K:max.K) 
system.time(
  for(K in min.K:max.K) 
  { # K <- 2 
    # ################################################################################### #
    # MEAN SQUARED ERROR FOR THE FRMM FINAL PARTITION WITH K = K: MeanSquaredError_K...RData #
    # Load Objects3..., FinalPartitFRMMandTotalSumResid100runs_K, for each of the clusters
    # in the FRMM partition (part.frmm.Ktemp.mat), fit an FDboost model,  with Y=y.summer.m.transf.loess, 
    # x1, x2, x3 defined in the cluster, get the predicted values of each observation (gene) 
    # for the model fitted to that cluster: predicted.temp.testing.k; estimate the squared  
    # residuals (squared.residuals.l2.mat), and finally, find the mean squared error of the 
    # partition, using these residuals. MeanSquaredError_K4-5_79seasonR_20-4-20.RData
    # ################################################################################### #
    cat("# ############ # K = ", K, "# ############ #\n")
    load(paste0(gitPath, dataPath, "Objects3_K", K, "_100RunsRealSeasonalDataSet_79_20-4-20.RData") )
    load(paste0(gitPath, dataPath, "FinalPartitFRMMandTotalSumResid100runs_K", K, "_79seasonR_20-4-20.RData") ) 
    cat("# ############ # K = ", K, "# ############ #\n")
    m.train <- as.vector(table(part.frmm.Ktemp.mat) )
    genelab.clu.j.list <- list() # Gene labels for genes in the k-th cluster
    for(k in 1:K)
    {
      genelab.clu.j.list[[k]] <- rownames(part.frmm.Ktemp.mat)[part.frmm.Ktemp.mat==k ] 
    }
    # mean squared error: using observations of the same cluster in each testing set #
    squared.residuals.l2.mat <- matrix(NA, nrow=ngen, ncol=nT)
    rownames(squared.residuals.l2.mat) <- genelab.m 
    colnames(squared.residuals.l2.mat) <- timep # 1:24
    seed2 <- 5 
    for(k in 1:K) # n.knots=6 # 
    {
      mtemp <- m.train[k] # no. observations in this cluster
      dataset.k.train.x1x2x3.loess.list <- list()
      temp <- y.summer.m.transf.loess[genelab.clu.j.list[[k]], ] # mtemp x 24 time points
      dataset.k.train.x1x2x3.loess.list$Ytemp <- temp
      temp1 <- x1.spring.m.transf.loess[genelab.clu.j.list[[k]], ]
      dataset.k.train.x1x2x3.loess.list$x1temp.centred <- scale(temp1, center=TRUE, scale=FALSE)
      temp2 <- x2.autumn.m.transf.loess[genelab.clu.j.list[[k]], ]
      dataset.k.train.x1x2x3.loess.list$x2temp.centred <- scale(temp2, center=TRUE, scale=FALSE)
      temp3 <- x3.winter.m.transf.loess[genelab.clu.j.list[[k]], ]
      dataset.k.train.x1x2x3.loess.list$x3temp.centred <- scale(temp3, center=TRUE, scale=FALSE)
      dataset.k.train.x1x2x3.loess.list$timep <- timep
      dataset.k.train.x1x2x3.loess.list$s.timep <- timep
      dataset.k.train.x1x2x3.loess.list$df.temp5 <- df.temp5
      dataset.k.train.x1x2x3.loess.list$genelabs.k.clu <- genelab.clu.j.list[[k]]
      set.seed(seed2) # E(Yi|Xi)=beta0(t)+int beta2(s, t)*x2(s)*ds+int beta3(s, t)*x3(s)*ds 
      a1 <- try(fof.functional.temp <- FDboost( Ytemp ~ 1 + bsignal(x=x1temp.centred, s=s.timep, knots=n.knots, df = df.temp5) +
        bsignal(x=x2temp.centred, s=s.timep, knots=n.knots, df = df.temp5) + bsignal(x=x3temp.centred, s=s.timep, knots=n.knots, df = df.temp5),
        timeformula = ~ bbs(timep, knots=n.knots, df = df.temp5), data = dataset.k.train.x1x2x3.loess.list, control = boost_control(mstop = 300)) )
      if(class(a1)[1]=="try-error")
      { 
        cat("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n")
        cat(" k = ", k, "th cluster or model. FDboost command with\n",  " no.knots = ",
            n.knots, " reports an error. It is \n", a1[1])
        cat("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n")
        n.knots2 <- 7
        b <- try(fof.functional.temp <- FDboost( Ytemp ~ 1 + bsignal(x=x1temp.centred, s=s.timep, knots=n.knots2, df = df.temp5) +
                                                   bsignal(x=x2temp.centred, s=s.timep, knots=n.knots2, df = df.temp5) +
                                                   bsignal(x=x3temp.centred, s=s.timep, knots=n.knots2, df = df.temp5),
                                                 timeformula = ~ bbs(timep, knots=n.knots2, df = df.temp5),
         data = dataset.k.train.x1x2x3.loess.list, control = boost_control(mstop = 300)) )
        if(class(b)[1]=="try-error")
          stop("Investigate this. With 7 knots, there is an error too.\n")
      } # end if class of (a1)[1]==try-error # ########################################### #
      # # ################################################################################ #
      predicted.temp.testing.k <- predict(fof.functional.temp, newdata=dataset.k.train.x1x2x3.loess.list) # 2094 24
      ntemp <- prod(dim(temp)) 
      if(prod(dim(squared.residuals.l2.mat[genelab.clu.j.list[[k]], ]))!=prod(dim((temp-predicted.temp.testing.k)^2)))
        stop("Dimension of squared.residuals.l2.mat differs from dimension of (Y_k-^Y_k)^2.\n")
      squared.residuals.l2.mat[genelab.clu.j.list[[k]], ] <- (temp-predicted.temp.testing.k)^2 # SQUARE OF EACH RESIDUAL ################### # 
      # ###################################################### #
      # ###################################################### #
    } # End for(K in 2:15) #
    explanation <- paste0(c("mse.K: calculated using equal training and testing sets, and these are formed by the
     observations (=genes) in each cluster, and fitting each FDboost model to those. FRMM partition found with K =", K, "\n") )
    print(mse.K <- sum(squared.residuals.l2.mat)/(ngen*nT)) 
    mse.K.vec[paste0("K", K)] <- mse.K
    save(mse.K, K, explanation, itemnumber, datetemp, ngen, nT, file=paste0(gitPath, dataPath, "MeanSquaredError_K", K, "_79seasonR_20-4-20.RData"))
    # ###################################################### #
  } # End of for(K in 2:15) # ############################# #
)

rm(list = ls(all = TRUE)) # this bit below: added, 12-11-20
path1 <- c(gitPath, dataPath)
min.K <- 2
max.K <- 15
mse.K.vec <- rep(NA, max.K-min.K+1)
names(mse.K.vec) <- paste0("K", min.K:max.K)
for(K in min.K:max.K) # MeanSquaredError_K2_79seasonR_20-4-20.RData # K <- 2
{
  load(paste0(path1, "MeanSquaredError_K", K, "_79seasonR_20-4-20.RData") )
  mse.K.vec[paste0("K", K)] <- mse.K
  cat("K = ", K, "\n")
}
mse.K.vec
#        K2        K3        K4        K5        K6        K7        K8        K9       K10 
# 0.5661322 0.5044359 0.4570869 0.4268043 0.3991912 0.3851744 0.3716133 0.3692003 0.3538955 
#       K11       K12       K13       K14       K15 
# 0.3559160 0.3437409 0.3379414 0.3324962 0.3301616 
explanation <- paste0(c("mse.K: calculated using equal training and testing sets, and these are formed by the
 observations (=genes) in each cluster, and fitting each FDboost model to those. FRMM partition found with K =", K, " 12-11-20 \n") )
itemnumber <- c("79B")
save(mse.K.vec, K, explanation, itemnumber, datetemp, ngen, nT, file=paste0(path1, "MeanSquaredError_K", K, "_79seasonR_20-11-12.RData"))

