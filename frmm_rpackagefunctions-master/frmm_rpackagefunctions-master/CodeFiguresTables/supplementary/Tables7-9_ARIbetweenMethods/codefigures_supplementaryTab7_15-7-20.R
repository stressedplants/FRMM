# codefigures_supplementaryTab7ARIbetweenMethodsSeasonK10_83AUsingY_15-7-20.R #

# ######################################################################### #
# ######################################################################### #
# ########################### TABLE 7 #################################### #
# Table of ARIs between partitions, real seasonal data set, 23 time points, using Y
# ######################################################################### #
# ######################################################################### #
# 
# 83A).# remaining: from "Objects3_K10_100RunsRealSeasonalDataSet_79A_24-6-20.RData" back. #
rm(list = ls(all = TRUE)) 
load("~/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/Objects3_K10_100RunsRealSeasonalDataSet_79A_24-6-20.RData")
source("functions_8-8-19.R") #dim(y.summer.m.transf.loess) # [1] 5378   23#y.summer.m.transf.loess[c(1:2, ngen), c(1:2, nT)]
timep.hours2 <- c(seq(16, 22, by=1), seq(0, 15, by=1)) # 16 17 18 ... 14 15; length=23
nickname <- c("RealSeasonSetK10.23TimeP")
# ###############################################################
# ##
# ALTERNATIVE METHOD 2, 3: HDDC. #
if(!require(fda)) install.packages("fda"); library(fda)
if(!require(HDclassif)) install.packages("HDclassif"); library(HDclassif)
if(!require(clues)) install.packages("clues"); library(clues) 
temp <- y.summer.m.transf.loess 
K1 <- K # 10 #temp <- cbind(y.summer.m.transf.loess, x1.spring.m.transf.loess, x2.autumn.m.transf.loess, x3.winter.m.transf.loess)
set.seed(K) 
system.time( # warning in prms: Quick-TRANSfer stage steps exceeded maximum (= 268900) 
 prms <- hddc(temp, K=K1, kmeans.control=list(iter.max=100, nstart=50)) ) 
system.time(prms1 <- hddc(temp, K=K1, d_select = "BIC", kmeans.control=list(iter.max=100, nstart=50)) )
table(prms$class) # Y: 360  519  412  440   11  367  300  327 2356  286 #
table(prms1$class) #Y: 318  426  345  581 1207  604  590  603  302  402
table(prms$class, prms1$class)
part.alt2.default.mat <- matrix(prms$class) 
rownames(part.alt2.default.mat) <- genelab.m
colnames(part.alt2.default.mat) <- nickname
part.alt3.bic.mat <- part.alt2.default.mat
part.alt3.bic.mat[, 1] <-  prms1$class
# Alternative 4-5: FunHDDC # ####################################################### #
if(!require(funHDDC)) install.packages("funHDDC"); library(funHDDC) # Y variable is in temp at the moment
part.alt4.default.mat <- matrix(NA, nrow=ngen, ncol=1)
rownames(part.alt4.default.mat) <- rownames(temp) # == genelab.m, checked
colnames(part.alt4.default.mat) <- nickname
part.alt5.bic.mat <- part.alt4.default.mat
# Estimate a smooth of Y
# no. of basis functions = order + no. of interior knots. # nbasis = nbreaks + norder - 2
# 12 order 4 B-splines with equally spaced knots (i.e. piecewise cubic polynomials) over [1, 23]. 
# n. basis functions = 4 + 8 = 12. We will thus have 12
temp1 <- y.summer.m.transf.loess
rng <- c(1, 23) 
norder <- 4
nT <- 23 
bspl.bas.no <- 12 
wbasis <- create.bspline.basis(rangeval=rng, nbasis=bspl.bas.no, norder) 
tempfd <- smooth.basis(timep, t(temp1), wbasis)$fd 
names(tempfd$fdnames) <- c("time", "gene labels", "gene expression") 
K1 <- K # = 10
# ############################################################# #
# # Plot of the estimated smooth of Y
# #plot(wbasis) 
# par(las=1, cex.axis=2, cex.lab=2) # increase left margin:
# par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) # default: 5, 4, 4, 2 + 0.1 # par(mar = c(12, 5, 4, 2)+ 0.1) # 5.5, 6.35, 1, 2)+0.1
# plot(tempfd, xlim=c(1, 24), xaxt="none", ylim=c(min.y, max.y), xlab="t", ylab="Y") 
# ############################################################# #
set.seed(1)
prms4b.temp <- funHDDC(data=tempfd, K=K1, d_select="Cattell", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50))  
system.time(
if(length(warnings())>0) # All models diverged
{
    for(K1 in K:1) # CHECK WHEN USING Y, ALT. 4 #
    { #
      cat("K1 =", K1, "\n")
      if(i==1) imax <- 4 else 
        imax <- 5
      for(i in 1:imax)
      { # check 5 times whether Alt. 4, 5 found a partition with K1 clusters; if not, check for K1-1 #
        cat("i =", i, "\n")
        prms4.temp <- funHDDC(data=tempfd, K=K1, d_select="Cattell", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) 
        if(length(table(prms4.temp$class) )==0) next else   # or if(!is.null(prms4.temp$class))
        {
          cat(paste0("\tI found K1=", K1, " clusters, funHDDC, Y X, Alt 4.\n") )
          cat("# ############################################################# #")
          break
        }  
      } # end for i = 1:5 #
      if(length(table(prms4.temp$class) )==0) next else break
    } # end for K1 in K:1 USING Y, ALT. 4 # ####################### #
} # end if length(warnings())>0
) 
# ############################################################### #
# model with 3 clusters.
# ##
warnings() # 31, all models diverged, and 25: Quick-TRANSfer stage steps exceeded maximum (= 268900)
# ##################################################################################### #
# K1 = 3
# i1 = 1 	Length of warnings() =  0
# i1 = 2 	      model K threshold complexity           BIC
#         1 AKJBKQKDK 3       0.2        191 -1,199,986.08
# SELECTED: model  AKJBKQKDK  with  3  clusters.
# Selection Criterion: BIC. # Length of warnings() =  0
# Warning message:  In funHDDC(data = tempfd, K = K1, mc.cores = 5) : All models diverged.
# ##################################################################################### #
table(prms4.temp$class) # 2234 1871 1273
prms4 <- prms4.temp
part.alt4.default.mat[, 1] <- prms4$class
# End of Alt. 4 Method # ########################################################## #
# Alt.5 Method # ################################################################## #
K1 <- K 
set.seed(2)
prms5.temp <- funHDDC(data=tempfd, K=K1, d_select = "BIC", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) # "All models diverged" 
system.time( # All models diverged above. Then: #
if(length(warnings())>0) # All models diverged
  {
    for(K1 in K:1) 
    { #
      cat("K1 =", K1, "\n")
      for(i in 1:5)
      { # check 5 times whether Alt. 4, 5 found a partition with K1 clusters; if not, check for K1-1 #
        cat("i =", i, "\n")
        set.seed(1+(i-1)*30) # assign("last.warning", NULL, envir = baseenv()) # clear the previous warning. Then length(warnings)=0 #https://stackoverflow.com/questions/3903157/how-can-i-check-whether-a-function-call-results-in-a-warning
        prms5.temp <- funHDDC(data=list(tempfd, temp2fd, temp3fd, temp4fd), K=K1, d_select="BIC", mc.cores=5, kmeans.control=list(iter.max=100, nstart=50)) # Y X
        if(length(table(prms5.temp$class) )==0) next else   # or if(!is.null(prms4.temp$class))
        {
          cat(paste0("\tI found K1=", K1, " clusters, funHDDC, Y, Alt 5.\n") )
          cat("# ############################################################# #")
          break
        }  
      } # end for i = 1:5 #
      if(length(table(prms5.temp$class) )==0) next else break
    } # end for K1 in K:1 USING Y, FunHDDC.BIC # ####################### #
  } # end if length(warnings())>0
) 
# ############################################################### #
# It found a model with 4 clusters, for i1=1
# ##
warnings() # 18: Quick-TRANSfer stage steps exceeded maximum (= 268900); rest: all models diverged.
# ##################################################################################### #
# K1 = 8. It found a model with 8 clusters instead of 10. 761  679  664  951  235  386  517 1185 # before setting seed: 572 663 1095  371  434  926  673  644
#       model K threshold complexity                BIC
# 1 AKJBKQKDK 4       bic      4,423  -4,409,675,021.08
# SELECTED: model  AKJBKQKDK  with  4  clusters.
# Selection Criterion: BIC.
# I found K1=4 clusters, funHDDC, Y X, Alt 5.
# ##################################################################################### #
assign("last.warning", NULL, envir = baseenv()) # clear the previous warnings; then there are 0 again
table(prms5.temp$class) #1277 1431 1862  808
prms5 <- prms5.temp
part.alt5.bic.mat[, 1] <- prms5$class
table(prms4$class, prms5$class)
#      1    2    3    4
# 1    0    0  643    1
# 2 1277 1431 1219  807
Ks.alt.4.5 <- c(length(table(prms4.temp$class)), length(table(prms5.temp$class))) # 2 4 here
save(part.alt2.default.mat, part.alt3.bic.mat, part.alt4.default.mat, part.alt5.bic.mat,
 file="/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/AlternativeMethodsPartitRealSeasonSet_83AseasonR_15-7-20.RData")

rm(list = ls(all = TRUE)) 
load("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/Objects3_K10_100RunsRealSeasonalDataSet_79A_24-6-20.RData")
source("functions_8-8-19.R")
load("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/Objects3_K10_100RunsRealSeasonalDataSet_79A_24-6-20.RData")
# FPCA.oracle: # ############################################################# #
# 12 order 4 B-splines with 10 equally spaced knots.
if(!require(fda)) install.packages("fda"); library(fda)
if(!require(HDclassif)) install.packages("HDclassif"); library(HDclassif)
if(!require(clues)) install.packages("clues"); library(clues) 
temp <- y.summer.m.transf.loess
rng <- c(1, 23)
norder <- 4
nT <- 23 
bspl.bas.no <- 12 
wbasis <- create.bspline.basis(rangeval=rng, nbasis=bspl.bas.no, norder) 
tempfd <- smooth.basis(timep, t(temp), wbasis)$fd 
names(tempfd$fdnames) <- c("time", "gene labels", "gene expression") # FPCs:
fpca.temp <- pca.fd(tempfd, nharm=bspl.bas.no)
a <- cbind(fpca.temp$values, cumsum(fpca.temp$values)/sum(fpca.temp$values))
colnames(a) <- c("Eigenvalue", "% variance")
rownames(a) <- 1:nrow(a)
a <- cbind(a, 1/a[, 2], fpca.temp$varprop, fpca.temp$values/sum(fpca.temp$values))
colnames(a)[3:5] <- c("Inverse", "PropVar1", "PropVar2")
cat("Variance of each fPC (=eigenvalue) and cumulative variance:\n")
print(a)
#      Eigenvalue % variance  Inverse     PropVar1  =  PropVar2
# 1  1.712448e+01  0.4766091 2.098156 4.766091e-01 4.766091e-01 
# 2  7.944225e+00  0.6977129 1.433254 2.211039e-01 2.211039e-01
# 3  5.693065e+00  0.8561625 1.168003 1.584495e-01 1.584495e-01
# 4  3.185066e+00  0.9448093 1.058415 8.864683e-02 8.864683e-02
# 5  1.554160e+00  0.9880647 1.012079 4.325543e-02 4.325543e-02
# 6  3.529468e-01  0.9978880 1.002117 9.823225e-03 9.823225e-03
# 
plot(1:12, a[, 3], type='o', pch=19, cex=0.6)
nharm <- 6 # after looking at the graph of the log10, see p. 109 in Ramsay's book
x <- matrix(1, 12-nharm, 2)
x[, 2] <- (nharm+1):12
x # design matrix with intercept (1st column) plus eigenvalue number (=x)
y1 <- log10(fpca.temp$values[(nharm+1):12])
c <- lsfit(x, y1, int=FALSE)$coef
neig <- 12
plot(1:neig, log10(fpca.temp$values[1:neig]), type='o', pch=19, cex=0.6,
     xlab = "Eigenvalue Number", ylab="Log10 Eigenvalue")
lines(1:neig, c[1] + c[2]*(1:neig), lty=2, col=2)
# select 4 FPCs at the moment.
# if selecting instead 1st 5 FPCs:
nharm2 <- 5 
x2 <- matrix(1, 12-nharm2, 2)
x2[, 2] <- (nharm2+1):12
x2 
y2 <- log10(fpca.temp$values[(nharm2+1):12])
c2 <- lsfit(x2, y2, int=FALSE)$coef
lines(1:neig, c2[1] + c2[2]*(1:neig), lty=2, col=4)
nharm3 <- 4 
x3 <- matrix(1, 12-nharm3, 2)
x3[, 2] <- (nharm3+1):12; x3 
y3 <- log10(fpca.temp$values[(nharm3+1):12])
c3 <- lsfit(x3, y3, int=FALSE)$coef
neig <- 12
plot(1:neig, log10(fpca.temp$values[1:neig]), type='o', pch=19, cex=0.6,
 xlab = "Eigenvalue Number", ylab="Log10 Eigenvalue")
lines(1:neig, c3[1] + c3[2]*(1:neig), lty=2, col=3)
# ################################################################################################################# #
# pick 4 FPCs.
harmfd <- fpca.temp$harmonics # class=fd = functional PCs 
harmvals <- eval.fd(timep, harmfd) 
npc.temp <- 4
pc.dataset.temp <- fpca.temp$scores[, c(1:npc.temp)]
rownames(pc.dataset.temp) <- rownames(temp) 
colnames(pc.dataset.temp) <- paste0("fPC", 1:npc.temp)
Q1 <- 100 
# K-means Q (=Q1=100) times:
kmeans.partitions.mat <- matrix(NA, nrow=ngen, ncol=Q1) 
rownames(kmeans.partitions.mat) <-  genelab.m
colnames(kmeans.partitions.mat) <- paste0("q", 1:Q1)
set.seed(555) .oracle
system.time(
  for(q in 1:Q1) 
  {
    # ################################################################################################## #
    # Given the set of npc.temp functional PC scores calculated from tempfd, K,
    # the n genes=5378, perform K-means Q1=100 times in the set of X fPCs. 
    # ################################################################################################### #
    
    cat("q = ", q, "\n")
    kmeans.pc.temp <- kmeans(pc.dataset.temp, centers=K, iter.max=100, nstart=50)
    kmeans.partitions.mat[, q] <- kmeans.pc.temp$cluster
  }
) 
# 5 warnings: 5: Quick-TRANSfer stage steps exceeded maximum (= 268900)
# ############################################################################ #
# 2nd part: consensus matrix, K-means in it, and this is the final partition found by
# Alt. 1.
system.time(
 B.consens <- consensus.matrix.fun(M.run3=M.run, final.clme.per.run.mat3=kmeans.partitions.mat, ngen3=ngen, genelab.m3=genelab.m, K3=K) 
) 
B.consens$B.consensus.counts[1:5, 1:5]
#           Ahg311263 Ahg311286 Ahg311317 Ahg311324 Ahg311338
# Ahg311263       100       100         3         0         0
# Ahg311286       100       100         3         0         0
# Ahg311317         3         3       100         0        97
# Ahg311324         0         0         0       100         0
# Ahg311338         0         0        97         0       100
# ################################################################################################### #
table(kmeans.partitions.mat[, 1], kmeans.partitions.mat[, 2]) # not exactly same partition
table(kmeans.partitions.mat[, 1], kmeans.partitions.mat[, 4]) # ditto.
#       1    2    3    4    5    6    7    8    9   10
# 1   274    0    1    0    0    0    0    0    0    0
# 2     0    0    0    0 1151    0    0    0    9    0
# 3     0    0    0    0    0  483    0    0    0    0
# 4     0  542    0    0    8    0    0    0    7    0
# 5     0    0    0    0    0    0    0    0  547    0
# 6     0    0  482    0    0    0    0    0    0    6
# 7     0    8    0    0    1    0    0    1    0  502
# 8     0    0    0    0    0    0    0  432    0    0
# 9     0    0    0    0    1    0  773    0    0    0
# 10    0    0    0  150    0    0    0    0    0    0
# ###################################################################### #
consensus.mat <- B.consens$B.consensus.mat 
set.seed(5)
system.time( 
  clme.kmeans <- kmeans(consensus.mat, centers=K, iter.max=100, nstart=50) 
) 
clme.kmeans$cluster[1:5] 
table(clme.kmeans$cluster) # K10 432  275 1162  774  488  548  150  554  483  512
table(kmeans.partitions.mat[, 1], clme.kmeans$cluster) 
part.alt1.4fpc.mat <- kmeans.partitions.mat[, 1, drop=FALSE] 
rownames(part.alt1.4fpc.mat) <- genelab.m 
colnames(part.alt1.4fpc.mat) <- paste0(c("K"), K)
# RData file with all the final partitions found by FuRMi, Alt. 1-5
path.saving <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/")
save(part.alt1.4fpc.mat, nT, K, genelab.m, df.temp5, bspl.bas.no, n.knots, norder, kmeans.partitions.mat, B.consens, timep, 
 y.summer.m.transf.loess, x1.spring.m.transf.loess, x2.autumn.m.transf.loess, x3.winter.m.transf.loess, timepoint2,
 file=paste0(path.saving, "Solutions_Alt1_Only_ConsensusMatRealSeasSetK10_83AseasonR.RData"))
load("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/AlternativeMethodsPartitRealSeasonSet_83AseasonR_15-7-20.RData") 
Ks.alt.4.5 <- c(3, 8) # ##################################################################### #
save(part.alt1.4fpc.mat, nT, K, genelab.m, df.temp5, bspl.bas.no, n.knots, norder, kmeans.partitions.mat, timep, 
 y.summer.m.transf.loess, x1.spring.m.transf.loess, x2.autumn.m.transf.loess, x3.winter.m.transf.loess, timepoint2,
 part.alt2.default.mat, part.alt3.bic.mat, part.alt4.default.mat, part.alt5.bic.mat, Ks.alt.4.5,
 file=paste0(path.saving, "Solutionsb_Alt1-5_PartitionsRealSeasSetK10_83AseasonR.RData"))
# ############################################################################################################################ #

rm(list = ls(all = TRUE)) 
load("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/FinalPartitFRMMandTotalSumResid100runs_K10_79AseasonR_24-6-20.RData")
load("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/Solutionsb_Alt1-5_PartitionsRealSeasSetK10_83AseasonR.RData")
nickname <- c("K10realSeaSet23TimePoints83AseasonR")
all.equal(rownames(part.frmm.Ktemp.mat), rownames(part.alt1.4fpc.mat), rownames(part.alt2.default.mat),
 rownames(part.alt3.bic.mat), rownames(part.alt4.default.mat), rownames(part.alt5.bic.mat), genelab.m) # TRUE. Good, as expected.
all.equal(rownames(part.alt4.default.mat), genelab.m) # TRUE
path.saving <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/")
save.image(file=paste0(path.saving, "Solutions_FRMM_Alt1-5PartitionsK10RealSeasSet23TimeP_83AseasonR.RData"))
# ##

rm(list = ls(all = TRUE)) # MATRIX OF ADJUSTED RAND INDICES BETWEEN THESE PARTITIONS #
load("/home/susanacon.oraclede/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/supplementary/Table7_ARIbetweenMethodsSeasoSet_K10UsingY/Solutions_FRMM_Alt1-5PartitionsK10RealSeasSet23TimeP_83AseasonR.RData")
ls()
partit.mat <- cbind(part.frmm.Ktemp.mat, part.alt1.4fpc.mat, part.alt2.default.mat,
 part.alt3.bic.mat, part.alt4.default.mat, part.alt5.bic.mat)
colnames(partit.mat) <- c("FRMM", paste0("Alt", 1:5))
partit.mat[1:5, ]
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
newnamestemp <- c("FuRMi", "FPCA", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
colnames(ari.frmm.all.alt.mat) <- newnamestemp
rownames(ari.frmm.all.alt.mat) <- newnamestemp
# round(ari.frmm.all.alt.mat, 3)
# ######################################## #
# greatest: between Alt1-Alt5 (0.295) and between Alt4-Alt5 (0.288) # K=10 clusters # updated with seeds set, 15-7-20
#                 FuRMi        FPCA HDDC.Cattell HDDC.BIC FunHDDC.Cattell FunHDDC.BIC
# FuRMi               1       0.187        0.009    0.031           0.124       0.175
# FPCA               NA       1.000       -0.005    0.047           0.236       0.295
# HDDC.Cattell       NA          NA        1.000    0.004           0.010       0.010
# HDDC.BIC           NA          NA           NA    1.000           0.022       0.011
# FunHDDC.Cattell    NA          NA           NA       NA           1.000       0.288
# FunHDDC.BIC        NA          NA           NA       NA              NA       1.000
# ######################################### #
library(xtable) 
options(xtable.floating = FALSE)
options(xtable.timestamp = "") 
a <- xtable(round(ari.frmm.all.alt.mat, 3))
digits(a) <- xdigits(a)
a



