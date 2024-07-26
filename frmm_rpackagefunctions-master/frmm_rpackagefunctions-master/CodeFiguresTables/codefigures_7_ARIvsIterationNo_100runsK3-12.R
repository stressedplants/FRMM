# codefigures_7_ARIvsIterationNo_100runsK3-12.R
# ######################################################################### #
# ######################################################################### #
# ########################### FIGURE 7 #################################### #
# ######################################################################### #
# ######################################################################### #

# ARI AGAINST ITERATION NUMBER FOR EACH RUN IN FuRMi, 100 RUNS; K=3, 12; N=500, 1000.

rm(list = ls(all = TRUE))
path1 <- c("https://gitlab.com/sconde778/frmm_rpackagefunctions/data/") # ? path1 <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/ClustersMixtureMod_SeasonalGeneExpression/Sols_HistogramsL2_25-9-19/n500/)
setwd(path1)  #path3 <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/")
load("data_Figure7_ARIvsIteration_K3-12n500-1000.RData") # ari.mat.list.K3n500, ari.mat.list.K3n1000, ari.mat.list.K12n500, ari.mat.list.K12n1000
pathtemp <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/")

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

l <- 1
ytemp <- ari.mat.list[[l]]["%ARI", c(2:mtemp)]
pdf(paste0(pathtemp,  "Fig3aARI_versus_IterationNo_K3n500L2AR1.pdf"), width=10, height=9.5)
par(las=1, cex.axis=2, cex.lab=2) # increase left margin:
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) # default: 5, 4, 4, 2 + 0.1 # par(mar = c(12, 5, 4, 2)+ 0.1)
plot(xtemp, ytemp, type='o', xlim=c(1, mtemp), xlab="", ylab="Adjusted Rand Index (*100)", pch=19, cex=cextemp,
     lwd=lwdtemp, col=colvec[l], ylim=c(min.y, max.y))
for(l in 2:M.run)
{
  points(xtemp, ari.mat.list[[l]]["%ARI", c(2:mtemp)], pch=19, col=colvec[l], cex=cextemp, lwd=lwdtemp, type='o')
} # StudyARIperIteration_78seasonR/Graph_PercARIperIteration_FRMM_K3L2n500AR1_78seasonR_11-5-20.pdf
dev.off()

ari.mat.list <- ari.mat.list.K3n1000
l <- 1
ytemp <- ari.mat.list[[l]]["%ARI", c(2:mtemp)]
pdf(paste0(pathtemp,  "Fig3bARI_versus_IterationNo_K3n1000L2AR1.pdf"), width=10, height=9.5)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) # default: 5, 4, 4, 2 + 0.1 # par(mar = c(12, 5, 4, 2)+ 0.1)
plot(xtemp, ytemp, type='o', xlim=c(1, mtemp), xlab="", ylab="Adjusted Rand Index (*100)", pch=19, cex=cextemp,
     lwd=lwdtemp, col=colvec[l], ylim=c(min.y, max.y))
for(l in 2:M.run)
  points(xtemp, ari.mat.list[[l]]["%ARI", c(2:mtemp)], pch=19, col=colvec[l], cex=cextemp, lwd=lwdtemp, type='o')
# StudyARIperIteration_78seasonR/Graph_PercARIperIteration_FRMM_K3L2n1000AR1_78seasonR_11-5-20.pdf
dev.off()

ari.mat.list <- ari.mat.list.K12n500
mtemp <- 70
min.y <- 0
max.y <- 42
xtemp <- c(2:mtemp)

l <- 1
ytemp <- ari.mat.list[[l]]["%ARI", c(2:mtemp)]
pdf(paste0(pathtemp,  "Fig3cARI_versus_IterationNo_K12n500L2AR1.pdf"), width=10, height=9.5)
par(las=1, cex.axis=2, cex.lab=2) # increase left margin:
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) # default: 5, 4, 4, 2 + 0.1 # par(mar = c(12, 5, 4, 2)+ 0.1)# 5.5, 6.35, 1, 2)+0.1
plot(xtemp, ytemp, type='o', xlim=c(1, mtemp), xlab="", ylab="Adjusted Rand Index (*100)", pch=19, cex=cextemp,
     lwd=lwdtemp, col=colvec[l], ylim=c(min.y, max.y))
for(l in 2:M.run)
  points(xtemp, ari.mat.list[[l]]["%ARI", c(2:mtemp)], pch=19, col=colvec[l], cex=cextemp, lwd=lwdtemp, type='o')
dev.off()

ari.mat.list <- ari.mat.list.K12n1000
l <- 1
ytemp <- ari.mat.list[[l]]["%ARI", c(2:mtemp)]
pdf(paste0(pathtemp,  "Fig3dARI_versus_IterationNo_K12n1000L2AR1.pdf"), width=10, height=9.5)
par(las=1, cex.axis=2, cex.lab=2) # increase left margin:
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) # default: 5, 4, 4, 2 + 0.1 # par(mar = c(12, 5, 4, 2)+ 0.1)# 5.5, 6.35, 1, 2)+0.1
plot(xtemp, ytemp, type='o', xlim=c(1, mtemp), xlab="", ylab="Adjusted Rand Index (*100)", pch=19, cex=cextemp,
     lwd=lwdtemp, col=colvec[l], ylim=c(min.y, max.y))
#indextemp <- achieved.convergence.run.vec[l]-1
#text(indextemp, ari.mat.list[[l]]["%ARI", achieved.convergence.run.vec[l]-1], labels=l, col=colvec[l], cex=2)
for(l in 2:M.run)
  points(xtemp, ari.mat.list[[l]]["%ARI", c(2:mtemp)], pch=19, col=colvec[l], cex=cextemp, lwd=lwdtemp, type='o')
# StudyARIperIteration_78seasonR/Graph_PercARIperIteration_FRMM_K12L2n1000AR1_78seasonR_16-4-20.pdf
dev.off()

# END OF 3) ARI AGAINST ITERATION NUMBER FOR EACH RUN IN FuRMi, 100 RUNS; K=3, 12; N=500, 1000.




# ######################################################################### #
# ######################################################################### #
# ########################### FIGURE 4 #################################### #
# ######################################################################### #
# ######################################################################### #

# 4) DISTRIBUTION OF ARI AGAINST THE NUMBER OF RUNS - FuRMi, 10-100 RUNS; K=12; N=1000, L2, AR(1).

rm(list = ls(all = TRUE))
pathtemp <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/")
load(paste0(pathtemp, "data_Figure4_DistributionARIagainstNoRunsFRMM_All10to100Plus5-15_82seasonR_21-5-20.RData") )


par(las=1, cex.axis=2, cex.lab=2) # increase left margin:
par(mgp=c(4.5, 1, 0), mar=c(6, 6, 1, 2)+0.1) # default: 5, 4, 4, 2 + 0.1 # par(mar = c(12, 5, 4, 2)+ 0.1) # 5.5, 6.35, 1, 2)+0.1
boxplot(perc.ari.finalfrmmpartition.mat.total12b, xlab="Number of runs",
        ylab="Adjusted Rand Index (*100)")  # boxplots_ARIvsNoRunsFRMM_Every10To100Plus5-15Runs_82seasonR_21-5-20.pdf

# END OF 4) DISTRIBUTION OF ARI AGAINST THE NUMBER OF RUNS - FuRMi, 10-100 RUNS; K=12; N=1000, L2, AR(1).


# ######################################################################### #
# ######################################################################### #
# ########################### FIGURE 5 #################################### #
# ######################################################################### #
# ######################################################################### #


# 5) Mean squared error (adding up the errors of the FDboost fittings in all the clusters) versus the number
# of clusters, real seasonal data set, 23 time points.

rm(list = ls(all = TRUE))
library(FDboost) 
min.K <- 2 
max.K <- 18
mse.K.vec <- rep(NA, max.K-min.K+1)
names(mse.K.vec) <- paste0("K", min.K:max.K) 
path2 <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/ClustersMixtureMod_SeasonalGeneExpression/SolsSimul1Replic_Jan20/Study_79seasonR_FRMM_RealSeasonalDataSet17-4-20/A_SeasonalSet12TimePoints/Plots_BetasFinalFRMM_Model/")
source("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/ClustersMixtureMod_SeasonalGeneExpression/R_functions/functions_8-8-19.R")
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
    # predicted.temp.testing.k, estimate the squared residuals, in squared.residuals.l2.mat,
    # and finally, find the mean squared error of the partition, using these residuals. Save
    # it in MeanSquaredError_K.._79seasonR_20-4-20.RData
    # ################################################################################### #
    # ################################################################################### #
    # ################################################################################### #
    cat("# ############ # K = ", K, "# ############ #\n")       
    load(paste0("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/Figure5MSEvsKseasonal/Objects3_K", K, "_100RunsRealSeasonalDataSet_79A_24-6-20.RData") )
    load(paste0("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/CodeFiguresPaper/data/Figure5MSEvsKseasonal/FinalPartitFRMMandTotalSumResid100runs_K", 
                K, "_79AseasonR_24-6-20.RData") ) 
    cat("# ############ # K = ", K, "# ############ #\n")
    m.train <- as.vector(table(part.frmm.Ktemp.mat) ) 
    rownames(part.frmm.Ktemp.mat)[part.frmm.Ktemp.mat==1][1:5] 
    genelab.clu.j.list <- list() # Gene labels for genes in the k-th cluster:
    for(k in 1:K)
    {
      genelab.clu.j.list[[k]] <- rownames(part.frmm.Ktemp.mat)[part.frmm.Ktemp.mat==k ]
    }
    # mean squared error: using observations of the same cluster in each testing set #
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
          stop("Investigate this. With 7 knots, there is an error too.\n")
      } # end if class of (a1)[1]==try-error # ########################################### #
      # # ################################################################################ #
      predicted.temp.testing.k <- predict(fof.functional.temp, newdata=dataset.k.train.x1x2x3.loess.list) # 2094 24
      ntemp <- prod(dim(temp)) 
      if(prod(dim(squared.residuals.l2.mat[genelab.clu.j.list[[k]], ]))!=prod(dim((temp-predicted.temp.testing.k)^2)))
        stop("Dimension of squared.residuals.l2.mat differs from dimension of (Y_k-^Y_k)^2.\n")
      squared.residuals.l2.mat[genelab.clu.j.list[[k]], ] <- (temp-predicted.temp.testing.k)^2 # SQUARE OF EACH RESIDUAL ################### # 
      # ###################################################### #
    } # End for(K in 2:15) #
    explanation <- paste0(c("mse.K: calculated using equal training and testing sets, and these are formed by the
    observations (=genes) in each cluster, and fitting each FDboost model to those. FRMM partition found with K =",
                            K, "\n") )
    print(mse.K <- sum(squared.residuals.l2.mat)/(ngen*nT)) 
    mse.K.vec[paste0("K", K)] <- mse.K
    path.saving <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/ClustersMixtureMod_SeasonalGeneExpression/SolsSimul1Replic_Jan20/Study_79seasonR_FRMM_RealSeasonalDataSet17-4-20/A_SeasonalSet12TimePoints/")
    save(mse.K, K, explanation, itemnumber, datetemp, ngen, nT, file=paste0(path.saving, "MeanSquaredError_K", K, "_79AseasonR_24-6-20.RData"))
    # ###################################################### #
  } # End of for(K in 2:18) # ############################# # 
)

mse.vec <- mse.K.vec # c(mse.vec, mse.K.vec)
mse.vec # after the update where we set the seed before the kmeans call:
#        K2        K3        K4        K5        K6        K7        K8        K9       K10       K11       K12       K13       K14       K15       K16       K17       K18 
# 0.9467330 0.8322436 0.7631034 0.7313826 0.6860181 0.6663203 0.6452519 0.6183390 0.6014889 0.6029251 0.5987175 0.5783899 0.5678991 0.5658196 0.5497500 0.5572379 0.5530234
save(mse.vec, file=paste0(path.saving, "MeanSquaredErrorTotal_K2-", K, "_79AseasonR_15-7-20.RData")) # A_SeasonalSet12TimePoints

round(mse.vec, 5)
# Plot of the MSEs by number of clusters K.
x <- c(2:K) 
temp <- mse.vec
min.temp <- min(temp)
max.temp <- max(temp)
path.plot <- paste0(path.saving, "Plots_BetasFinalFRMM_Model/")

# ####################################################### #
pdf(paste0(path.plot,  "lineplot_MSE_FRMMrealSeasonSet_K2-", K, "_79AseasonR_15-7-20.pdf"), width=10, height=10)
par(las=1, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(4., 5., 1, 2)+0.1) 
plot(x, temp, type='o', pch=19, cex=0.6, lwd=3,
     ylim=c(min.temp, max.temp), xaxt="none", xlab="", ylab="", col="darkgreen")
axis(1, x, cex.axis=1.4) # lineplot_MSE_FRMMrealSeasonSet_K2-12_79AseasonR_24-6-20.pdf
dev.off()
# ####################################################### #
# ####################################################### #
