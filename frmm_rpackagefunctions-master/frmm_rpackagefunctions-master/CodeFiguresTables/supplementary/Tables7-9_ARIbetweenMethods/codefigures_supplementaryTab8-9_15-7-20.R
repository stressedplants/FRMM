# codefigures_supplementaryTab8-9ARIbetweenMethodsSeasonK10_83BUsingYX-Y_15-7-20.R #
# plus Figure 10.

# ################################################################################# #
# ################################################################################# #
# ########################### TABLE 8 ############################################# #
# Table of ARIs between partitions, real seasonal data set, 24 time points, using Y X
# ################################################################################# #
# ################################################################################# #
# Within data in Gitlab:
#Objects3_K10_100RunsRealSeasonalDataSet_79_20-4-20.RData, FinalPartitFRMMandTotalSumResid100runs_K10_79BseasonR_20-4-20.RData,
#AlternMethods2-5RealSeasonSet_83BseasonR_YX_21-7-20_24TimeP.RData, AlternativeMethodsPartitRealSeasonSet_83AseasonR_15-7-20.RData,
#AlternMethods2-5RealSeasonSet_83BseasonR_Y_25-7-20_24TimeP.RData, AlternMethods2-5RealSeasonSet_83BseasonR_Y_25-7-20_24TimeP.RData,
#SolsBrief_FinalPartition_Alt1WithY_X_MFPCA_RealSeasonalSetK10n5378_21-7-20_83BseasonR_24TimeP.RData,
#SolsBrief_FinalPartition_Alt1WithY_FPCA_RealSeasonalSetK10n5378_25-7-20_83BseasonR_24TimeP.RData,
#ConsensusMat_FinalPartition_Alt1WithY_FPCA4_RealSeasonalSetK10_25-7-20_83BseasonR_TimeP.RData,
#Solutions_FRMM_Alt1-5PartitionsK10RealSeasSet24TimeP_YX_83BseasonR_21-7-20.RData, Solutions_FRMM_Alt1-5PartitionsK10RealSeasSet24TimeP_Y_83BseasonR_25-7-20.RData
# In box.com:
#AlternativeMethodsPartitRealSeasonSet_All_83AseasonR_15-7-20.RData, AlternMethods2-5All_RealSeasonSet_All_83BseasonR_Y_25-7-20_24TimeP.RData,    
# remaining but probably not needed (almost half GB): Alt1_All_withY_K10_FPCA_83BseasonR_25-7-20_24TimeP.RData.
#
rm(list = ls(all = TRUE)) 
K <- 10
load(paste0(gitPath, dataPath, "Objects3_K", K, "_100RunsRealSeasonalDataSet_79_20-4-20.RData") )
source("functions_8-8-19.R")
dim(y.summer.m.transf.loess) # [1] 5378   24
timep.hours <- c(seq(16, 22, by=2), seq(0, 22, by=2), seq(0, 14, by=2))
nickname <- c(paste0("RealSeasonSetK", K, ".24TimeP") )
# ############ HDDC.CATTELL, HDDC.BIC # ##################################################
if(!require(fda)) install.packages("fda"); library(fda)
if(!require(HDclassif)) install.packages("HDclassif"); library(HDclassif)
if(!require(clues)) install.packages("clues"); library(clues) # adjustedRand
K1 <- K # = 10 here
temp <- cbind(y.summer.m.transf.loess, x1.spring.m.transf.loess, x2.autumn.m.transf.loess, x3.winter.m.transf.loess)
set.seed(K+5)
system.time(
 prms <- hddc(temp, K=K1, kmeans.control=list(iter.max=100, nstart=50)) )
system.time(prms1 <- hddc(temp, K=K1, d_select = "BIC", kmeans.control=list(iter.max=100, nstart=50)) )
table(prms$class) 
table(prms1$class) 
table(prms$class, prms1$class)
# ##
part.alt2.default.mat <- matrix(prms$class) 
rownames(part.alt2.default.mat) <- genelab.m
colnames(part.alt2.default.mat) <- nickname
part.alt3.bic.mat <- part.alt2.default.mat
part.alt3.bic.mat[, 1] <-  prms1$class
# Alternative 4-5: # ####################################################### #
if(!require(funHDDC)) install.packages("funHDDC"); library(funHDDC) 
part.alt4.default.mat <- matrix(NA, nrow=ngen, ncol=1)
rownames(part.alt4.default.mat) <- rownames(temp) 
colnames(part.alt4.default.mat) <- nickname
part.alt5.bic.mat <- part.alt4.default.mat
temp1 <- y.summer.m.transf.loess
rng <- c(1, 24) 
norder <- 4
nT <- 24 
bspl.bas.no <- 12 
wbasis <- create.bspline.basis(rangeval=rng, nbasis=bspl.bas.no, norder) 
tempfd <- smooth.basis(timep, t(temp1), wbasis)$fd 
names(tempfd$fdnames) <- c("time", "gene labels", "gene expression") 
temp2 <- x1.spring.m.transf.loess
temp3 <- x2.autumn.m.transf.loess
temp4 <- x3.winter.m.transf.loess
temp2fd <- smooth.basis(timep, t(temp2), wbasis)$fd
temp3fd <- smooth.basis(timep, t(temp3), wbasis)$fd
temp4fd <- smooth.basis(timep, t(temp4), wbasis)$fd
K1 <- K # = 10
# ############################################################# #
# # Plot of the estimated smooth of Y=gene expression at summer:
# #plot(wbasis) 
# par(las=1, cex.axis=2, cex.lab=2) # increase left margin:
# par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) # par(mar = c(12, 5, 4, 2)+ 0.1) # 5.5, 6.35, 1, 2)+0.1
# plot(tempfd, xlim=c(1, 24), xaxt="none", ylim=c(min.y, max.y), xlab="t", ylab="Y") 
# plot(timep, temp[i, ], type='o', pch=19, cex=0.5, xaxt="none", xlab="t", ylim=c(min.y, max.y), ylab="Y")
# axis(1, timep, labels=timep.hours, cex.axis=1.2) #
# ############################################################# #
# ##
usingY <- 0 # 1 Y, 0 YX usingY <- 1 # for Y
set.seed(6)
if(usingY==1) # All models diverged here below for any usingY: #
 prms4.temp <- funHDDC(data=tempfd, K=K1, d_select="Cattell", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) else 
if(usingY==0) prms4.temp <- funHDDC(data=list(tempfd, temp2fd, temp3fd, temp4fd), K=K1, d_select="Cattell", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) 
system.time( 
if(length(warnings())>0) 
{
 for(K1 in K:1) 
 { #
  cat("K1 =", K1, "\n")
  for(i1 in 1:5)
  { 
   cat("i =", i1, "\n")
   set.seed(2+(i1-1)*50)
   if(usingY==1) prms4.temp <- funHDDC(data=tempfd, K=K1, d_select="Cattell", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) else
   if(usingY==0) prms4.temp <- funHDDC(data=list(tempfd, temp2fd, temp3fd, temp4fd), K=K1, d_select="Cattell", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50))
   if(length(table(prms4.temp$class) )==0) next else   
   {
    cat(paste0("\tI found K1=", K1, " clusters, funHDDC, Y X, Alt 4.\n") )
    cat("# ############################################################# #")
    break
   }  
  } # end for i = 1:5 #
  if(length(table(prms4.temp$class) )==0) next else break
 } # end for K1 in K:1, ALT. 4 # ####################### #
} # end if length(warnings())>0
) # model with 2 clusters.
warnings() # YX: 26, all models diverged # Y: K = 3 clusters, 35 warnings. All are "All models diverged".
assign("last.warning", NULL, envir = baseenv())
part.alt4.default.mat[, 1] <- prms4.temp$class 
prms4 <- prms4.temp
#######################################################
############################# #
# K1 = 5 FOR Y, X
# i1 = 2 # model K threshold complexity BIC
# 1 AKJBKQKDK 5 0.2 1,439 -837,901.65
# SELECTED: model AKJBKQKDK with 5 clusters.
# Selection Criterion: BIC.
# I found K1=5 clusters, funHDDC, Y X, Alt 4.#
#######################################################
############################# #
# K1 = 3 FOR USING Y ONLY
# i1 = 1 
# model K threshold complexity BIC
# 1 AKJBKQKDK 3 0.2 198 -1,500,495.61
# 
# SELECTED: model AKJBKQKDK with 3 clusters.
# Selection Criterion: BIC.
# I found K1=3 clusters, funHDDC, Y X, Alt 4.
# End of Alt. 4 Method #
#######################################################
## #
# End of FunHDDC.Cattell Method # ########################################################## #
# FunHDDC.BIC # ################################################################## #
K1 <- K # set it again
set.seed(8)
if(usingY==1) 
 prms5.temp <- funHDDC(data=tempfd, K=K1, d_select = "BIC", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) else # "All models  diverged" 
if(usingY==0) 
 prms5.temp <- funHDDC(data=list(tempfd, temp2fd, temp3fd, temp4fd), K=K1, d_select = "BIC", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) # "All models diverged" 
# All models diverged above. Then: #
if(length(warnings())>0) 
{
 for(K1 in K:1) 
 { #
  cat("K1 =", K1, "\n")
  for(i in 1:5)
  { # check 5 times whether Alt. 4, 5 found a partition with K1 clusters; if not, check for K1-1 #
   cat(" i =", i, "\n")
   set.seed(5+(i-1)*20)  
   if(usingY==1) prms5.temp <- funHDDC(data=tempfd, K=K1, d_select="BIC", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) else
   prms5.temp <- funHDDC(data=list(tempfd, temp2fd, temp3fd, temp4fd), K=K1, d_select="BIC", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) #Y X
   if(length(table(prms5.temp$class) )==0) next else 
   {
    cat(paste0("\tI found K1=", K1, " clusters, funHDDC, Y X, Alt 5.\n") )
    cat("# ############################################################ #")
    break
   } 
  } # end for i = 1:5 #
  if(length(table(prms5.temp$class) )==0) next else break
 } # end for K1 in K:1 USING Y X, funHDDC.Cattell # ####################### #
} # end if length(warnings())>0
# ############################################################### #
#
#######################################################
####### #
warnings() # YX: 10: Quick-TRANSfer stage steps exceeded maximum (= 268900). Rest: All models diverged.
#######################################################
############################# #
# K1 = 6, YX. It found a model with 6 clusters. 514 1682 1300 751 759 372
# model K threshold complexity BIC
# 1 AKJBKQKDK 6 bic 6,175 -1,386,204,826.22
# SELECTED: model AKJBKQKDK with 6 clusters.
# Selection Criterion: BIC.
#######################################################
############################# #
# K = 9, Y: model K threshold complexity BIC
# 1 AKJBKQKDK 9 bic 692 -100,206,479.63
# SELECTED: model AKJBKQKDK with 9 clusters.
# Selection Criterion: BIC.
# I found K1=9 clusters, funHDDC, Y X, Alt 5. Only warnings: all models diverged
assign("last.warning", NULL, envir = baseenv()) # clear the previous warnings; then there are 0 again
table(prms5.temp$class)
prms5 <- prms5.temp
part.alt5.bic.mat[, 1] <- prms5$class
table(prms4$class, prms5$class)
# 1 2 3 4 5 6
# 1 72 752 1172 972 233 262# 2 0 0 0 109 0 0
# 3 10 44 0 0 94 0
# 4 579 126 0 8 5 0
# 5 89 9 0 840 0 2
# 1 2 3 4 5 6 7 8 9 # Y
# 1 95 295 156 548 0 0 257 448 664
# 2 0 37 0 0 0 0 0 0 0
# 3 0 0 647 892 1134 202 3 0 0

Ks.alt.4.5 <- c(length(table(prms4.temp$class)), length(table(prms5.temp$class))) # 5 6
save(part.alt2.default.mat, part.alt3.bic.mat, part.alt4.default.mat, part.alt5.bic.mat,
 file=paste0(gitPath, dataPath, "AlternativeMethodsPartitRealSeasonSet_83AseasonR_15-7-20.RData"))
save.image(
 file=paste0(gitPath, dataPath, "AlternativeMethodsPartitRealSeasonSet_All_83AseasonR_15-7-20.RData"))
save(part.alt2.default.mat, part.alt3.bic.mat, part.alt4.default.mat,
 part.alt5.bic.mat, Ks.alt.4.5, 
 file=paste0(gitPath, dataPath, "AlternMethods2-5RealSeasonSet_83BseasonR_Y_25-7-20_24TimeP.RData"))
save.image( 
 file=paste0(gitPath, dataPath, "AlternMethods2-5All_RealSeasonSet_All_83BseasonR_Y_25-7-20_24TimeP.RData"))
save(part.alt2.default.mat, part.alt3.bic.mat, part.alt4.default.mat, part.alt5.bic.mat, Ks.alt.4.5,
file=paste0(gitPath, dataPath, "AlternMethods2-5RealSeasonSet_83BseasonR_Y_25-7-20_24TimeP.RData") )
save.image( 
file=paste0(gitPath, dataPath, "AlternMethods2-5All_RealSeasonSet_All_83BseasonR_Y_25-7-20_24TimeP.RData"))
#
# ################################ #
# End of Alt. 2-5 here. (85B)
# ################################ #
# FPCA.oracle: # ############################################################# #
# ALT 1 using Y X.
# ##
rm(list = ls(all = TRUE)) 
K <- 10
load(paste0(gitPath, dataPath, "Objects3_K", K, "_100RunsRealSeasonalDataSet_79_20-4-20.RData"))
source(gitPath, functionsPath, "functions_8-8-19.R")
if(!require(fda)) install.packages("fda"); library(fda)
if(!require(HDclassif)) install.packages("HDclassif"); library(HDclassif)
if(!require(clues)) install.packages("clues"); library(clues)
library(MFPCA)
# FPC per variable
temp <- y.summer.m.transf.loess 
rng <- c(1, nT) 
norder <- 4 
bspl.bas.no <- 12 
wbasis <- create.bspline.basis(rangeval=rng, nbasis=bspl.bas.no, norder) 
tempfd <- smooth.basis(timep, t(temp), wbasis)$fd 
names(tempfd$fdnames) <- c("time", "gene labels", "gene expression") 
fpca.temp <- pca.fd(tempfd, nharm=bspl.bas.no)
a <- cbind(fpca.temp$values, cumsum(fpca.temp$values)/sum(fpca.temp$values), fpca.temp$varprop)
colnames(a) <- c("Eigenvalue", "% variance", "Prop.var")
rownames(a) <- 1:nrow(a)
cat("Variance of each fPC (=eigenvalue) and cumulative variance:\n")
print(a)
#      Eigenvalue % variance     Prop.var
# 1  1.115777e+01  0.5115512 5.115512e-01
# 2  3.887754e+00  0.6897934 1.782422e-01
# 3  2.814294e+00  0.8188206 1.290272e-01
# 4  2.286158e+00  0.9236343 1.048137e-01
# 5  1.437625e+00  0.9895452 6.591092e-02 # pick 4 FPCs for Y as they explain >= 90% variability.
temp1 <- x1.spring.m.transf.loess
tempfd1 <- smooth.basis(timep, t(temp1), wbasis)$fd 
fpca.temp <- pca.fd(tempfd1, nharm=bspl.bas.no)
a <- cbind(fpca.temp$values, cumsum(fpca.temp$values)/sum(fpca.temp$values), fpca.temp$varprop)
colnames(a) <- c("Eigenvalue", "Cum%variance", "PropVar")
rownames(a) <- 1:nrow(a)
print(a)
#      Eigenvalue Cum%variance      PropVar
# 1  6.495193e+00    0.6468585 6.468585e-01
# 2  1.428653e+00    0.7891385 1.422800e-01
# 3  9.818949e-01    0.8869257 9.778725e-02
# 4  6.132502e-01    0.9479995 6.107380e-02 # pick 4 FPCs for x1 as they explain >=90%
temp1 <- x2.autumn.m.transf.loess
tempfd1 <- smooth.basis(timep, t(temp1), wbasis)$fd 
fpca.temp <- pca.fd(tempfd1, nharm=bspl.bas.no)
a <- cbind(fpca.temp$values, cumsum(fpca.temp$values)/sum(fpca.temp$values), fpca.temp$varprop)
colnames(a) <- c("Eigenvalue", "Cum%variance", "PropVar")
rownames(a) <- 1:nrow(a)
print(a)
#      Eigenvalue Cum%variance      PropVar
# 1  6.192066e+00    0.5668594 5.668594e-01
# 2  1.962338e+00    0.7465037 1.796443e-01
# 3  1.122647e+00    0.8492775 1.027739e-01
# 4  9.951932e-01    0.9403836 9.110603e-02 # pick 4 FPCs for x2
temp1 <- x3.winter.m.transf.loess
tempfd1 <- smooth.basis(timep, t(temp1), wbasis)$fd 
fpca.temp <- pca.fd(tempfd1, nharm=bspl.bas.no)
a <- cbind(fpca.temp$values, cumsum(fpca.temp$values)/sum(fpca.temp$values), fpca.temp$varprop)
colnames(a) <- c("Eigenvalue", "Cum%variance", "PropVar")
rownames(a) <- 1:nrow(a)
print(a)
#      Eigenvalue Cum%variance      PropVar
# 1  7.655210e+00    0.7713119 7.713119e-01
# 2  1.160176e+00    0.8882072 1.168952e-01
# 3  4.975285e-01    0.9383364 5.012922e-02 # pick 3 FPCs for x3
# Summary: pick 4, 4, 4, 3 individual FPCs for Y, x1, x2, x3 respectively. # ##
YfunData <- funData(argvals = 1:nT, X = y.summer.m.transf.loess) 
x1funData <- funData(argvals = 1:nT, X = x1.spring.m.transf.loess)
x2funData <- funData(argvals = 1:nT, X = x2.autumn.m.transf.loess)
x3funData <- funData(argvals = 1:nT, X = x3.winter.m.transf.loess)
Y.X.funData <- multiFunData(YfunData, x1funData, x2funData, x3funData) 
Y.funData <- multiFunData(YfunData) 
selected.individual.fpcs <- c(4, 4, 4, 3) 
a <- selected.individual.fpcs
uniExpansions.temp <- list(list(type = "uFPCA", npc = a[1]), list(type = "uFPCA", npc = a[2]),  
 list(type = "uFPCA", npc = a[3]), list(type = "uFPCA", npc = a[4])) 
number.mfpc <- 15 
set.seed(10) 
MFPCA.Y.X <- MFPCA(Y.X.funData, M = number.mfpc, uniExpansions = uniExpansions.temp) 
round(summary(MFPCA.Y.X), 5); 
no.selected.mfpc <- 8 # 8 explain >= 90%; see below # ################################################################ #
# dim(MFPCA.Y.X$scores); # 5378 15
# MFPCA.Y.X$scores[1:10, ]; 
# screeplot(MFPCA.Y.X); #scoreplot(MFPCA.Y.X); 
#15 multivariate functional principal components estimated with 4 elements, each.
#* * * * * * * * * *                     
#                                      PC 1     PC 2     PC 3     PC 4     PC 5     PC 6     PC 7     PC 8     PC 9    PC 10    PC 11    PC 12    PC 13    PC 14    PC 15
# Eigenvalue                       14.51196 11.80930  5.53528  4.44234  3.12415  2.28610  1.52388  1.20836  1.18300  0.97727  0.71579  0.56906  0.39328  0.30365  0.01519
# Proportion of variance explained  0.29861  0.24300  0.11390  0.09141  0.06428  0.04704  0.03136  0.02486  0.02434  0.02011  0.01473  0.01171  0.00809  0.00625  0.00031
# Cumulative proportion             0.29861  0.54160  0.65550  0.74691  0.81120  0.85824  0.88959  0.91446  0.93880  0.95891  0.97364  0.98535  0.99344  0.99969  1.00000
# 15 MFPCs, 8 explain > 90% of the total variability.
# ##################################################################################################################################### #
# ##################################################################################################################################### #
# Next: K-means 100 times with MFPCA.Y.X$scores[, 1:no.selected.mfpc] i.e. 8 1st columns
usingY <- 1 
if(usingY==1) 
{
 pc.dataset.temp <- fpca.temp$scores[, c(1:4)] 
 colnames(pc.dataset.temp) <- paste0("FPC", 1:4) 
} else if(usingY==0)
{
 pc.dataset.temp <- MFPCA.Y.X$scores[, c(1:no.selected.mfpc)] 
 colnames(pc.dataset.temp) <- paste0("MFPC", 1:no.selected.mfpc) 
}
rownames(pc.dataset.temp) <- genelab.m 
Q1 <- 100 
kmeans.partitions.mat <- matrix(NA, nrow=ngen, ncol=Q1) 
rownames(kmeans.partitions.mat) <-  genelab.m
colnames(kmeans.partitions.mat) <- paste0("q", 1:Q1)
set.seed(11+K)
for(q in 1:Q1) 
{
  # ################################################################################################## #
  # Given MFPCA.Y.X$scores[, c(1:no.selected.mfpc)], K,
  # the n genes=5378, perform K-means Q1 times. Finally compute the consensus matrix, K-means, and 
  # this is the partition found by Alt. 1 method.
  # ################################################################################################### #
  if(!(q%%10) ) cat("q =", q, "\n")
  kmeans.temp <- kmeans(pc.dataset.temp, centers=K, iter.max=100, nstart=50)
  kmeans.partitions.mat[, q] <- kmeans.temp$cluster
} # End for q=   # 4 warnings: 4: Quick-TRANSfer stage steps exceeded maximum (= 268900) # using Y only: 1 of these warnings
ys# ########################### # # I think that these are those Daphne said not to be worried.
table(kmeans.partitions.mat[, 1], kmeans.partitions.mat[, 3]) # same
table(kmeans.partitions.mat[, 1], kmeans.partitions.mat[, 2]) # same
table(kmeans.partitions.mat[, 1], kmeans.partitions.mat[, 100]) 
kmeans.partitions.mat[1:5, 1:10] # USING Y ONLY: NOT EXACTLY SAME PARTITIONS
# ############################################################################ #
# 2nd part: 
system.time(
B.consens <- consensus.matrix.fun(M.run3=M.run, final.clme.per.run.mat3=kmeans.partitions.mat, 
 ngen3=ngen, genelab.m3=genelab.m, K3=K)
)  
B.consens$B.consensus.counts[1:5, 1:5]
B.consens$B.consensus.mat[1:5, 1:5]
#           Ahg311263 Ahg311286 Ahg311317 Ahg311324 Ahg311338
# Ahg311263       100        84         2        84         3 # not all partitions are the same #
# Ahg311286        84       100        18       100         3
# Ahg311317         2        18       100        18         0
# Ahg311324        84       100        18       100         3
# Ahg311338         3         3         0         3       100
#           Ahg311263 Ahg311286 Ahg311317 Ahg311324 Ahg311338
# Ahg311263      1.00      0.84      0.02      0.84      0.03 # ditto. 
# Ahg311286      0.84      1.00      0.18      1.00      0.03
# Ahg311317      0.02      0.18      1.00      0.18      0.00
# Ahg311324      0.84      1.00      0.18      1.00      0.03
# Ahg311338      0.03      0.03      0.00      0.03      1.00
# ################################################################################################### #
print(a <- sum(as.vector(B.consens$B.consensus.mat)==0) / 2 ) # 12712646 # USING Y ONLY: 12652924
print(b <- (ngen^2-ngen) / 2 ) # 14458753
a/b*100 
summary(as.vector(B.consens$B.consensus.mat)) 
length(table(B.consens$B.consensus.mat)) 
consensus.mat.temp <- B.consens$B.consensus.mat
set.seed(6)
kmeans.from.Qkmeans.temp <- kmeans(consensus.mat.temp, centers=K, iter.max=100, nstart=50)
if(usingY==1) part.alt1.Y.4fpc <- matrix(kmeans.from.Qkmeans.temp$cluster) 
part.alt1.temp <- part.alt1.Y.4fpc; rownames(part.alt1.temp) <- genelab.m
part.alt1.Y.X.8mfpc.4.4.4.3 <- matrix(kmeans.from.Qkmeans.temp$cluster)
rownames(part.alt1.Y.X.8mfpc.4.4.4.3) <- genelab.m
optimal.FPCs.MFPCs.combination <- c("i8m4.4.4.3") # 8 MFPCs; rest: no. of individual FPCs #
explanation <- c("Alt 1, partition with 4, 4, 4, 3, individual FPCs (which explain each >= 90% of total variability) and
 8 MFPCs (which also explain >= 90% of total variability), 83BseasonR,
 part.alt1.Y.X.8mfpc.4.4.4.3: partition of Alt. 1, optimal.FPCs.MFPCs.combination: optimal
 combination e.g. i8m4.4.4.3 = 4, 4, 4, 3 individual FPCs resp. per variable and 8 MFPCs.")
save(explanation, part.alt1.Y.X.8mfpc.4.4.4.3, optimal.FPCs.MFPCs.combination,  
 y.summer.m.transf.loess, x1.spring.m.transf.loess, x2.autumn.m.transf.loess, x3.winter.m.transf.loess,
 tempfd, ngen, nT, genelab.m, bspl.bas.no, K,
 file=paste0(gitPath, dataPath, "SolsBrief_", "FinalPartition_Alt1WithY_X_MFPCA_RealSeasonalSetK10n5378_21-7-20_83BseasonR_24TimeP.RData"))
# consensus.mat.temp: half GB.
# ###################### USING Y ONLY: ########################################################### #
part.alt1.Y.4fpc <- part.alt1.temp
explanation <- c("Alt 1, USING Y ONLY, partition with 4 individual FPCs (which explain each >= 90% of total variability) of Y, 83BseasonR,
 part.alt1.Y.4fpc: partition of Alt. 1, 25-7-20.")
save(explanation, part.alt1.Y.4fpc, y.summer.m.transf.loess, x1.spring.m.transf.loess, x2.autumn.m.transf.loess, x3.winter.m.transf.loess,
 tempfd, ngen, nT, genelab.m, bspl.bas.no, K,
 file=paste0(gitPath, dataPath, "SolsBrief_", "FinalPartition_Alt1WithY_FPCA_RealSeasonalSetK10n5378_25-7-20_83BseasonR_24TimeP.RData"))
save(consensus.mat.temp, 
 file=paste0(gitPath, dataPath, "ConsensusMat_", "FinalPartition_Alt1WithY_FPCA4_RealSeasonalSetK10_25-7-20_83BseasonR_TimeP.RData"))
save.image(file=paste0(gitPath, dataPath, "Alt1_All_withY_K10_FPCA_83BseasonR_25-7-20_24TimeP.RData")) 
# ###################################################################### #
# ##
rm(list = ls(all = TRUE))
load(paste0(gitPath, dataPath, "FinalPartitFRMMandTotalSumResid100runs_K10_79BseasonR_20-4-20.RData") ) 
load(paste0(gitPath, dataPath, "AlternMethods2-5RealSeasonSet_83BseasonR_YX_21-7-20_24TimeP.RData") ) 
load(paste0(gitPath, dataPath, "SolsBrief_FinalPartition_Alt1WithY_X_MFPCA_RealSeasonalSetK10n5378_21-7-20_83BseasonR_24TimeP.RData"))   
nickname <- paste0("K", K, "realSeaSet", nT, "TimePointsYX.83BseasonR") # nT=24, K = 10 # K10realSeaSet24TimePointsYX.83BseasonR
ls() 
all.equal(rownames(part.frmm.Ktemp.mat), rownames(part.alt1.Y.X.8mfpc.4.4.4.3), rownames(part.alt2.default.mat),
 rownames(part.alt3.bic.mat), rownames(part.alt4.default.mat), rownames(part.alt5.bic.mat), genelab.m) 
all.equal(rownames(part.alt4.default.mat), genelab.m) # TRUE
save.image(file=paste0(paste0(gitPath, dataPath, "Solutions_FRMM_Alt1-5PartitionsK10RealSeasSet24TimeP_YX_83BseasonR_21-7-20.RData")))

rm(list = ls(all = TRUE)) # MATRIX OF ADJUSTED RAND INDICES BETWEEN THESE PARTITIONS #
load(paste0(gitPath, dataPath, "Solutions_FRMM_Alt1-5PartitionsK10RealSeasSet24TimeP_YX_83BseasonR_21-7-20.RData") )
ls()
newnamestemp <- c("FuRMi", "FPCA", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
partit.mat <- cbind(part.frmm.Ktemp.mat, part.alt1.Y.X.8mfpc.4.4.4.3, part.alt2.default.mat,
 part.alt3.bic.mat, part.alt4.default.mat, part.alt5.bic.mat)
colnames(partit.mat) <- newnamestemp
partit.mat[1:5, ]
#           FuRMi FPCA HDDC.Cattell HDDC.BIC FunHDDC.Cattell FunHDDC.BIC
# Ahg311263     8    3           10        7               4           1
# Ahg311286     8    3            3        7               5           1
# Ahg311317    10    7            8        7               5           4
library(clues)
adjustedRand(partit.mat[, 1], partit.mat[, 2])["HA"]
ari.frmm.all.alt.mat <- matrix(NA, nrow=6, ncol=6)
rownames(ari.frmm.all.alt.mat) <- colnames(partit.mat)
colnames(ari.frmm.all.alt.mat) <- colnames(partit.mat)
diag(ari.frmm.all.alt.mat) <- 1
for(l2 in 1:5) 
{
 for(l3 in (l2+1):6)
  ari.frmm.all.alt.mat[l2, l3] <- adjustedRand(partit.mat[, l2], partit.mat[, l3])["HA"]
} 
ari.frmm.all.alt.mat 
round(ari.frmm.all.alt.mat, 3)
library(xtable) 
options(xtable.floating = FALSE) 
options(xtable.timestamp = "")
a <- xtable(round(ari.frmm.all.alt.mat, 3))
digits(a) <- xdigits(a)
a
# #################################################################################### #
#                 FuRMi       FPCA HDDC.Cattell   HDDC.BIC FunHDDC.Cattell  FunHDDC.BIC
# FuRMi               1 0.07149641  0.009528053 0.06025190      0.02320135  0.041879458
# FPCA               NA 1.00000000  0.018860070 0.46271835      0.12514861  0.258700609
# HDDC.Cattell       NA         NA  1.000000000 0.02873931      0.02028330 -0.002094342
# HDDC.BIC           NA         NA           NA 1.00000000      0.09620990  0.180304425
# FunHDDC.Cattell    NA         NA           NA         NA      1.00000000  0.141801086
# FunHDDC.BIC        NA         NA           NA         NA              NA  1.000000000
#                 FuRMi  FPCA HDDC.Cattell HDDC.BIC FunHDDC.Cattell FunHDDC.BIC
# FuRMi               1 0.071        0.010    0.060           0.023       0.042
# FPCA               NA 1.000        0.019    0.463           0.125       0.259
# HDDC.Cattell       NA    NA        1.000    0.029           0.020      -0.002
# HDDC.BIC           NA    NA           NA    1.000           0.096       0.180
# FunHDDC.Cattell    NA    NA           NA       NA           1.000       0.142
# FunHDDC.BIC        NA    NA           NA       NA              NA       1.000
# #################################################################################### #
# usingY=1:
rm(list = ls(all = TRUE)) 
usingY <- 1 
load(paste0(gitPath, dataPath, "FinalPartitFRMMandTotalSumResid100runs_K10_79BseasonR_20-4-20.RData")) 
load(paste0(gitPath, dataPath, "AlternMethods2-5RealSeasonSet_83BseasonR_Y_25-7-20_24TimeP.RData"))
load(paste0(gitPath, dataPath, "SolsBrief_FinalPartition_Alt1WithY_FPCA_RealSeasonalSetK10n5378_25-7-20_83BseasonR_24TimeP.RData"))
nickname <- paste0("K", K, "realSeaSet", nT, "TimePointsY.83BseasonR") 
ls() 
save.image(file=paste0(gitPath, dataPath, "Solutions_FRMM_Alt1-5PartitionsK10RealSeasSet24TimeP_Y_83BseasonR_25-7-20.RData")) 

rm(list = ls(all = TRUE)) # MATRIX OF ADJUSTED RAND INDICES BETWEEN THESE PARTITIONS #
load(paste0(gitPath, dataPath, "Solutions_FRMM_Alt1-5PartitionsK10RealSeasSet24TimeP_Y_83BseasonR_25-7-20.RData"))
ls()
newnamestemp <- c("FuRMi", "FPCA", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
partit.mat <- cbind(part.frmm.Ktemp.mat, part.alt1.Y.4fpc, part.alt2.default.mat,
 part.alt3.bic.mat, part.alt4.default.mat, part.alt5.bic.mat)
colnames(partit.mat) <- newnamestemp
partit.mat[1:5, ]
#           FuRMi FPCA HDDC.Cattell HDDC.BIC FunHDDC.Cattell FunHDDC.BIC
# Ahg311263     8    2            7        2               2           2
# Ahg311286     8    2           10        2               1           2
# Ahg311317    10    1            7        3               3           4
library(clues)
adjustedRand(partit.mat[, 1], partit.mat[, 2])["HA"] 
ari.frmm.all.alt.mat <- matrix(NA, nrow=6, ncol=6)
rownames(ari.frmm.all.alt.mat) <- colnames(partit.mat)
colnames(ari.frmm.all.alt.mat) <- colnames(partit.mat)
diag(ari.frmm.all.alt.mat) <- 1
system.time(
for(l2 in 1:5)
{
 for(l3 in (l2+1):6)
  ari.frmm.all.alt.mat[l2, l3] <- adjustedRand(partit.mat[, l2], partit.mat[, l3])["HA"]
} 
)
ari.frmm.all.alt.mat 
round(ari.frmm.all.alt.mat, 3)
library(xtable) 
options(xtable.floating = FALSE) 
options(xtable.timestamp = "") 
a <- xtable(round(ari.frmm.all.alt.mat, 3))
digits(a) <- xdigits(a)
a
# ########################################################################################## #
#                 FuRMi      FPCA HDDC.Cattell   HDDC.BIC FunHDDC.Cattell  FunHDDC.BIC
# FuRMi               1 0.1987324  0.004225112 0.12460659     0.071765319 0.1579446037
# FPCA               NA 1.0000000  0.009500570 0.21116818     0.145647546 0.2512765363
# HDDC.Cattell       NA        NA  1.000000000 0.03281363     0.007765357 0.0009364944
# HDDC.BIC           NA        NA           NA 1.00000000     0.100239380 0.1193002157
# FunHDDC.Cattell    NA        NA           NA         NA     1.000000000 0.1738609420
# FunHDDC.BIC        NA        NA           NA         NA              NA 1.0000000000
#                 FuRMi  FPCA HDDC.Cattell HDDC.BIC FunHDDC.Cattell FunHDDC.BIC
# FuRMi               1 0.199        0.004    0.125           0.072       0.158
# FPCA               NA 1.000        0.010    0.211           0.146       0.251
# HDDC.Cattell       NA    NA        1.000    0.033           0.008       0.001
# HDDC.BIC           NA    NA           NA    1.000           0.100       0.119
# FunHDDC.Cattell    NA    NA           NA       NA           1.000       0.174

# ################################################################################# #
# ################################################################################# #
# ########################### FIGURE 10 ########################################### #
# ################################################################################# #
# ################################################################################# #
rm(list = ls(all = TRUE))
library(FDboost) 
min.K <- 2
max.K <- 15
mse.K.vec <- rep(NA, max.K-min.K+1)
names(mse.K.vec) <- paste0("K", min.K:max.K) # to save plot with betas:
gitPath <- c('~/git/frmm_rpackagefunctions/') # SHAHIN
#gitPath = c("https://gitlab.com/sconde778/frmm_rpackagefunctions/") ## Susana # CodeFiguresPaper/GitlabRepository/frmm_rpackagefunctions/
#gitPath = ## Daphne
## relative paths
dataPath <- 'data/'
figPath <- 'figures/'
source(paste0(gitPath, "R/functions1_3-6-20.R"))

  for(K in min.K:max.K)
  { 
    # ################################################################################### #
    # ################################################################################### #
    # ################################################################################### #
    # MEAN SQUARED ERROR FOR THE FRMM FINAL PARTITION WITH K = K: MeanSquaredError_K...RData #
    # Load Objects3..., FinalPartitFRMMandTotalSumResid100runs_K, for K = k,
    # for each of the clusters in the FRECL partition (part.frmm.Ktemp.mat), fit an FDboost model,
    # using Y=y.summer.m.transf.loess, x1, x2, x3 defined in the cluster,
    # get the predicted values of each observation (gene) for the model fitted to that cluster:
    # predicted.temp.testing.k, estimate the squared residuals, in squared.residuals.l2.mat,
    # and finally, find the mean squared error of the partition, using these residuals. Save
    # it in MeanSquaredError_K..._79seasonR_20-4-20.RData
    # ################################################################################### #
    # ################################################################################### #
    # ################################################################################### #
    load(paste0(gitPath, dataPath, "Objects3_K", K, "_100RunsRealSeasonalDataSet_79_20-4-20.RData") )
    if(K==10) load(paste0(gitPath, dataPath, "FinalPartitFRMMandTotalSumResid100runs_K",
      K, "_79BseasonR_20-4-20.RData") ) else
    load(paste0(gitPath, dataPath, "FinalPartitFRMMandTotalSumResid100runs_K",
     K, "_79seasonR_20-4-20.RData") ) # y.summer.m.transf.loess, x1.spring.m.transf.loess, x2.autumn.m.transf.loess, (5378 x 24), x3.winter.m.transf.loess, timep, M.run
    # part.frmm.Ktemp.mat, K, genelab.m, ngen, sum.l2.resid.convergclme.per.run.vec #part.frmm.Ktemp.mat <- part.frmm.K4.mat
    cat("# ############ # K = ", K, "# ############ #\n")
    table(part.frmm.Ktemp.mat) 
    m.train <- as.vector(table(part.frmm.Ktemp.mat) )
    part.frmm.Ktemp.mat[1:5, ]
    rownames(part.frmm.Ktemp.mat)[part.frmm.Ktemp.mat==1][1:5] # ""Ahg311317" ...
    genelab.clu.j.list <- list() # Gene labels for genes in the k-th cluster:
    names(timep) <- colnames(y.summer.m.transf.loess)
    for(k in 1:K)
    {
     genelab.clu.j.list[[k]] <- rownames(part.frmm.Ktemp.mat)[part.frmm.Ktemp.mat==k ] # "Ahg311317" ... 2094
    }
    # mean squared error #
    squared.residuals.l2.mat <- matrix(NA, nrow=ngen, ncol=nT)
    rownames(squared.residuals.l2.mat) <- genelab.m # "Ahg311263" ...
    colnames(squared.residuals.l2.mat) <- timep 
    seed2 <- 5 
    for(k in 1:K) # n.knots=6 # FOR k-th CLUSTER #
    {
      mtemp <- m.train[k] 
      dataset.k.train.x1x2x3.loess.list <- list()
      temp <- y.summer.m.transf.loess[genelab.clu.j.list[[k]], ] # mtemp x 24
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
      set.seed(seed2) # E(Yi|Xi)=beta0(t)+int beta2(s, t)*x2(s)*ds+int beta3(s, t)*x3(s)*ds # seed2=5 
      a1 <- try(fof.functional.temp <- FDboost( Ytemp ~ 1 + bsignal(x=x1temp.centred, s=as.vector(s.timep), knots=n.knots, df = df.temp5) +
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
          stop("Investigate this. With 7 knots, there is an error too.\n")
      } # end if class of (a1)[1]==try-error # ########################################### #
      # # ################################################################################ #
      predicted.temp.testing.k <- predict(fof.functional.temp, newdata=dataset.k.train.x1x2x3.loess.list) 
      ntemp <- prod(dim(temp)) 
      if(prod(dim(squared.residuals.l2.mat[genelab.clu.j.list[[k]], ]))!=prod(dim((temp-predicted.temp.testing.k)^2)))
        stop("Dimension of squared.residuals.l2.mat differs from dimension of (Y_k-^Y_k)^2.\n")
      squared.residuals.l2.mat[genelab.clu.j.list[[k]], ] <- (temp-predicted.temp.testing.k)^2 # SQUARE OF EACH RESIDUAL ################### # 2094 x 24
      # ###################################################### #
      # # Y = temp = y.summer.m.transf.loess[genelab.clu.j.list[[k]], ]
      # # and the predicted values Yhat = predicted.temp.testing.k, both n_k x nT where n_K=no. observations in the k-th cluster.
      # # For the i-th observation (=gene) in the k-th cluster, temp[i, ] = Y_i, x1[i, ]=x1_i and so on.
      # ###################################################### #
    } # End for(K in 2:15) #
    explanation <- paste0(c("mse.K: calculated using equal training and testing sets, and these are formed by the
  observations (=genes) in each cluster, and fitting each FDboost model to those. FRMM partition found with K =",
                            K, "\n") )
    print(mse.K <- sum(squared.residuals.l2.mat)/(ngen*nT)) # K5: 0.4268043 K4: 0.4570869 # 0.5044359 # 0.5661322
    mse.K.vec[paste0("K", K)] <- mse.K
    path.saving <- paste0(gitPath, dataPath)
    save(mse.K, K, explanation, itemnumber, datetemp, ngen, nT, file=paste0(path.saving,
                                                                            "MeanSquaredError_K", K, "_79seasonR_20-4-20.RData"))
    # ###################################################### #
  } # End of for(K in 2:15) # ############################# #
)

# Plot of the MSEs by number of clusters K.
x <- c(2:K) 
print(temp <- mse.K.vec)
#        K2        K3        K4        K5        K6        K7        K8        K9       K10       K11       K12       K13       K14       K15 
# 0.5661322 0.5044359 0.4570869 0.4268043 0.3991912 0.3851744 0.3716133 0.3692003 0.3538955 0.3559160 0.3437409 0.3379414 0.3324962 0.3301616
min.temp <- min(temp)
max.temp <- max(temp)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) 
plot(x, temp, type='o', pch=19, cex=0.6, lwd=3,
     ylim=c(min.temp, max.temp), xaxt="none", xlab="", ylab="", col="darkgreen")
axis(1, x, cex.axis=2) # lineplot_MSE_FRMMrealSeasonSet_K2-15_79seasonR_15-7-20.pdf
# ####################################################### #
# ####################################################### #



# FunHDDC.BIC        NA    NA           NA       NA              NA       1.000
# ############################################################################################ #








