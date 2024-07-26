# codefigures_1_BoxplotsARIetc_50SimulK3L2iid_55seasonR_15-8-20.R #

# ######################################################################### #
# ######################################################################### #
# ########################### FIGURE 1 #################################### #
# Boxplots of the distribution of ARI*100, RI, etc., 50 simulations, K=3 clu.,
# L2 distance, i.i.d., all methods
# ######################################################################### #
# ######################################################################### #
# 
# remaining: from Res_50simulations_48seasonR_K3n500L2iids2_1_FuRMi_Alt1_15-8-20.RData, Res_Alt2-3, Res_Alt4-5 (48, 55) Res_50simulations_50seasonR_K3n1000L2iids2_1_FuRMi_Alt1_15-8-20.RData, etc. in 55, 56
# n=500 #
rm(list = ls(all = TRUE)) # ari.tpr.our.algorithm 4 x 50, ari.tpr.altern.fpc, K, ngen, genelab.m; ari.tpr.alt2.hddc.default, ari.tpr.alt2.hddc.bic
path1 <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/Project_Packrat/R_Files/")
setwd(path1)
# setwd("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/ClustersMixtureMod_SeasonalGeneExpression/Sols_HistogramsL2_25-9-19/iid_All/ReplicSumm/")
load("Res_50simulations_48seasonR_K3n500L2iids2_1_FuRMi_Alt1_15-8-20.RData")
load("Res_Alt2-3_K3L2n500iid_15-8-20.RData")
load("Res_Alt4-5_K3L2n500iid_15-8-20.RData")
#dim(ari.tpr.our.algorithm) # ari.tpr.our.algorithm[, c(1:3)]
mat <- cbind(ari.tpr.our.algorithm["%ARI", ], ari.tpr.altern.fpc["%ARI", ], ari.tpr.alt2.hddc.default["%ARI", ], 
 ari.tpr.alt2.hddc.bic["%ARI", ], ari.tpr.alt4.default.mat["%ARI", ], ari.tpr.alt5.bic.mat["%ARI", ])
colnames(mat) <- c("FuRMi", "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
path.plot <- path1
x <- 1:6
pdf(paste0(path.plot,  "boxplot_ARIby100_50simul_K", K, "L2distn", ngen, "_50seasonR_15-8-20.pdf"), width=12, height=9)
par(las=2, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(12, 5., 1, 2)+0.1) 
boxplot(mat, xaxt="none")
axis(1, at=x, labels=colnames(mat), cex.axis=1.6) 
dev.off()
labtemp <- c("%RI")
mat <- cbind(ari.tpr.our.algorithm[labtemp, ], ari.tpr.altern.fpc[labtemp, ], ari.tpr.alt2.hddc.default[labtemp, ], 
 ari.tpr.alt2.hddc.bic[labtemp, ], ari.tpr.alt4.default.mat[labtemp, ], ari.tpr.alt5.bic.mat[labtemp, ])
colnames(mat) <- c("FuRMi", "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
pdf(paste0(path.plot,  "boxplot_RIby100_50simul_K", K, "L2distn", ngen, "_50seasonR_15-8-20.pdf"), width=12, height=9)
par(las=2, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(12, 5., 1, 2)+0.1) 
boxplot(mat, xaxt="none")
axis(1, at=x, labels=colnames(mat), cex.axis=1.6) 
dev.off()
labtemp <- c("%TPR")
mat <- cbind(ari.tpr.our.algorithm[labtemp, ], ari.tpr.altern.fpc[labtemp, ], ari.tpr.alt2.hddc.default[labtemp, ], 
 ari.tpr.alt2.hddc.bic[labtemp, ], ari.tpr.alt4.default.mat[labtemp, ], ari.tpr.alt5.bic.mat[labtemp, ])
colnames(mat) <- c("FuRMi", "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
pdf(paste0(path.plot,  "boxplot_TPRby100_50simul_K", K, "L2distn", ngen, "_50seasonR_15-8-20.pdf"), width=12, height=9)
par(las=2, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(12, 5., 1, 2)+0.1) 
boxplot(mat, xaxt="none")
axis(1, at=x, labels=colnames(mat), cex.axis=1.6) 
dev.off()
labtemp <- c("%TNR")
mat <- cbind(ari.tpr.our.algorithm[labtemp, ], ari.tpr.altern.fpc[labtemp, ], ari.tpr.alt2.hddc.default[labtemp, ], 
 ari.tpr.alt2.hddc.bic[labtemp, ], ari.tpr.alt4.default.mat[labtemp, ], ari.tpr.alt5.bic.mat[labtemp, ])
colnames(mat) <- c("FuRMi", "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
pdf(paste0(path.plot,  "boxplot_TNRby100_50simul_K", K, "L2distn", ngen, "_50seasonR_15-8-20.pdf"), width=12, height=9)
par(las=2, cex.axis=2, cex.lab=2)
par(mgp=c(3, 1, 0), mar=c(12, 5., 1, 2)+0.1) 
boxplot(mat, xaxt="none")
axis(1, at=x, labels=colnames(mat), cex.axis=1.6) 
dev.off()

# n = 1000 # (55))
rm(list = ls(all = TRUE)) # ari.tpr.our.algorithm 4 x 50, ari.tpr.altern.fpc, K, ngen, genelab.m; ari.tpr.alt2.hddc.default, ari.tpr.alt2.hddc.bic
path1 <- c("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/Project_Packrat/R_Files/")
setwd(path1)
# setwd("/home/susanaconde/MyFilesDell/susana_Warwick/PapersWarwick/MyPapersWarwick/ClustersMixtureMod_SeasonalGeneExpression/Sols_HistogramsL2_25-9-19/iid_All/ReplicSumm/")
load("Res_50simulations_50seasonR_K3n1000L2iids2_1_FuRMi_Alt1_15-8-20.RData")
load("Res_50simulations_50seasonR_K3n1000L2iids2_1_Alt2-3_15-8-20.RData")
load("Res_Alt4-5_K3L2n1000iid_15-8-20.RData") # ari.tpr.alt4.default.mat, ari.tpr.alt5.bic.mat
# dim(ari.tpr.our.algorithm) 4 50 # ari.tpr.our.algorithm[, c(1:3)]
mat <- cbind(ari.tpr.our.algorithm["%ARI", ], ari.tpr.altern.fpc["%ARI", ], ari.tpr.alt2.hddc.default["%ARI", ], 
 ari.tpr.alt2.hddc.bic["%ARI", ], ari.tpr.alt4.default.mat["%ARI", ], ari.tpr.alt5.bic.mat["%ARI", ])
colnames(mat) <- c("FuRMi", "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
path.plot <- path1
x <- 1:6
pdf(paste0(path.plot,  "boxplot_ARIby100_50simul_K", K, "L2distn", ngen, "_50seasonR_15-8-20.pdf"), width=12, height=9)
par(las=2, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(12, 5., 1, 2)+0.1) 
boxplot(mat, xaxt="none")
axis(1, at=x, labels=colnames(mat), cex.axis=1.6)
dev.off()

labtemp <- c("%RI")
mat <- cbind(ari.tpr.our.algorithm[labtemp, ], ari.tpr.altern.fpc[labtemp, ], ari.tpr.alt2.hddc.default[labtemp, ], 
 ari.tpr.alt2.hddc.bic[labtemp, ], ari.tpr.alt4.default.mat[labtemp, ], ari.tpr.alt5.bic.mat[labtemp, ])
colnames(mat) <- c("FuRMi", "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
pdf(paste0(path.plot,  "boxplot_RIby100_50simul_K", K, "L2distn", ngen, "_50seasonR_15-8-20.pdf"), width=12, height=9)
par(las=2, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(12, 5., 1, 2)+0.1) 
boxplot(mat, xaxt="none")
axis(1, at=x, labels=colnames(mat), cex.axis=1.6) 
dev.off()
labtemp <- c("%TPR")
mat <- cbind(ari.tpr.our.algorithm[labtemp, ], ari.tpr.altern.fpc[labtemp, ], ari.tpr.alt2.hddc.default[labtemp, ], 
 ari.tpr.alt2.hddc.bic[labtemp, ], ari.tpr.alt4.default.mat[labtemp, ], ari.tpr.alt5.bic.mat[labtemp, ])
colnames(mat) <- c("FuRMi", "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
pdf(paste0(path.plot,  "boxplot_TPRby100_50simul_K", K, "L2distn", ngen, "_50seasonR_15-8-20.pdf"), width=12, height=9)
par(las=2, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(12, 5., 1, 2)+0.1) 
boxplot(mat, xaxt="none")
axis(1, at=x, labels=colnames(mat), cex.axis=1.6)
dev.off()
labtemp <- c("%TNR")
mat <- cbind(ari.tpr.our.algorithm[labtemp, ], ari.tpr.altern.fpc[labtemp, ], ari.tpr.alt2.hddc.default[labtemp, ], 
 ari.tpr.alt2.hddc.bic[labtemp, ], ari.tpr.alt4.default.mat[labtemp, ], ari.tpr.alt5.bic.mat[labtemp, ])
colnames(mat) <- c("FuRMi", "FPCA.oracle", "HDDC.Cattell", "HDDC.BIC", "FunHDDC.Cattell", "FunHDDC.BIC")
pdf(paste0(path.plot,  "boxplot_TNRby100_50simul_K", K, "L2distn", ngen, "_50seasonR_15-8-20.pdf"), width=12, height=9)
par(las=2, cex.axis=2, cex.lab=2) 
par(mgp=c(3, 1, 0), mar=c(12, 5., 1, 2)+0.1) 
boxplot(mat, xaxt="none")
axis(1, at=x, labels=colnames(mat), cex.axis=1.6) 
dev.off()






