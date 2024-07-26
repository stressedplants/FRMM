# functions1_test_3-6-20.R
# binary.matrix.truth.function (used in iterative.clustering.function, later),
# binary.vectors.lth.partition.function (used in proportions..., later), consensus.matrix.fun,
# iterative.clustering.function, prop.correct.pairs.lth.partition.and.truth.function requires R clues package (adjustedRand)

# rm(list = ls())

binary.matrix.truth.function <- function(K1=K, ngen1=ngen,
 true.tildetildeP1=true.tildetildeP)
{
  # ##################################################################### #
  # Given a true partition true.tildetildeP, find a binary matrix Atrue with
  # observations in rows and columns ordered by cluster such
  # that those in the same cluster appear together; 1 row and column
  # genes are clustered together 0 otherwise
  # ##################################################################### #
  if(length(true.tildetildeP1)!=ngen1)
    stop("Length of true.tildetildeP1 must be ngen1.\n")
  if(length(table(true.tildetildeP1))!=K1)
    stop("The partition true.tildetildeP1 should have exactly K1 clusters.\n")
  # ##################################################################### #
  # true.tildetildeP[1:5] # Ahg311263 Ahg311286 Ahg311317, 2, 1, 1, etc.
  binary.mat.true.partition <- matrix(0, nrow=ngen1, ncol=ngen1)
  rownames(binary.mat.true.partition) <- 1:ngen1
  colnames(binary.mat.true.partition) <- 1:ngen1
  m.gen.clu.temp <- table(true.tildetildeP1)
  m.gen.clu.temp.cum <- cumsum(m.gen.clu.temp) # 1793 3586 5378
  # ################################################################# #
  # added 3-6-20: if no names in true.tildetildeP1, we use the indices
  if(is.null(names(true.tildetildeP1) )) names(true.tildetildeP1) <- c(1:ngen1)
  # only null in example when doing the R package. Otherwise it always had names e.g. gene labels.
  # ################################################################# #
  k <- 1
  mtemp <- m.gen.clu.temp[k] # 1793
  rownames(binary.mat.true.partition)[1:mtemp] <-
    names(which(true.tildetildeP1==k)) # Ahg311286 ...
  #colnames(binary.mat.true.partition)[1:mtemp] <- names(which(true.tildetildeP==k))
  binary.mat.true.partition[1:mtemp, 1:mtemp] <- matrix(1, nrow=mtemp, ncol=mtemp)
  for(k in 2:K1) #  k <- 2
  {
    mtemp <- table(true.tildetildeP1)[k]
    rownames(binary.mat.true.partition)[(m.gen.clu.temp.cum[k-1]+1):(m.gen.clu.temp.cum[k])] <-
      names(which(true.tildetildeP1==k))
    binary.mat.true.partition[(m.gen.clu.temp.cum[k-1]+1):(m.gen.clu.temp.cum[k]),
                              (m.gen.clu.temp.cum[k-1]+1):(m.gen.clu.temp.cum[k])] <-
      matrix(1, nrow=mtemp, ncol=mtemp)
  }
  colnames(binary.mat.true.partition) <- rownames(binary.mat.true.partition)
  return(binary.mat.true.partition)
}


binary.vectors.lth.partition.function <- function(K2=K, ngen2=ngen,
  hatP.l2)
{
  # ##################################################################### #
  # Given a partition ^P(l)=hatP.l2, find a binary vector v_k[i] with
  # 1 ith observation (=gene) belongs to the kth cluster, 0 otherwise i.e. a "dummy" variable
  # for the kth cluster. Store the K2 vectors in v.binary.k.mat.temp by rows and
  # return this matrix, which is K2 x ngen2.
  # hatP.l2=final.clme.per.run.mat[, l]
  # ##################################################################### #
  if(length(hatP.l2)!=ngen2)
    stop("Length of hatP.l2 must be ngen2.\n")
  if(length(table(hatP.l2))!=K2) # modified, 7-2-20
    cat("binary.vectors...() message: the partition hatP.l2 does not have exactly K2 clusters. If K2> no. of clu. in hatP.l2 the returning matrix
in binary.vectors.lth.partition.function() will have some rows with all 0s, which indicate that 0 observations
(in the columns) belong to those clusters, i.e. they are empty.\n")
  # Example, 65), with K=12 clu. Alt. 4 and Alt. 5, 7-2-20, output matrix of this function:
  #    Ahg486775 Ahg940488 Ahg924569 ... (n=500 genes)
  # 1          0         1         1
  # 2          1         0         0
  # 3          0         0         0 # from this row on, all zeros
  # 4          0         0         0
  # 5          0         0         0
  # 6          0         0         0
  # 7          0         0         0
  # 8          0         0         0
  # 9          0         0         0
  # 10         0         0         0
  # 11         0         0         0
  # 12         0         0         0
  # ##################################################################### #
  # ##################################################################### #
  genelab.m.temp <- names(hatP.l2)
  v.binary.k.mat.temp <- matrix(0, nrow=K2, ncol=ngen2) # 3 x 5378
  rownames(v.binary.k.mat.temp) <- 1:K2
  colnames(v.binary.k.mat.temp) <- genelab.m.temp
  for(i in 1:ngen2) # i <- 1
  {
    # create v_k[i] and store it in v.binary.k.mat.temp[k, i]
    ktemp <- hatP.l2[i]
    v.binary.k.mat.temp[ktemp, i] <- 1
  } # end for ith gene
  # recall that usually, the colnames of v.binary.k.mat.temp do not have
  # to be in the same order as the one for "the true partition"
  return(v.binary.k.mat.temp)
}






consensus.matrix.fun <- function(M.run3=1, final.clme.per.run.mat3=final.clme.per.run.mat,
  ngen3=ngen, genelab.m3=genelab.m, K3=K)
{
  # ######################################################################### #
  # Given the M.run3 partitions from runs of the iterative algorithm (= no. of columns of final.clme.per.run.mat3)
  # find binary matrices A(l)=1 if row and column genes clustered together and 0 otherwise,
  # then get the consensus matrix B = sum(A(l)/M.run. Return a list with
  # $B.consensus.counts, $B.consensus.mat.
  # Input: FRMM M.run3 partitions, final.clme.per.run.mat3: matrix with dimensions ngen3 x M.run3 with the
  # partitions by columns, ngen3=no. of observations, genelab.m3 = character vector with the observation labels
  # in the order as they appear in the rows in final.clme.per.run.mat3, K3 = number of clusters in the partitions.
  # ############################################################################ #
  # Output: list with B.consensus.mat, ngen3 x ngen3, plus A(l)
  # Return a list with 3 components: each of the A(l) in $A.l.binary.hatP[[l]], plus the
  # consensus matrix with the counts and the consensus matrix with the proportions.
  # $B.consensus.counts and $ B.consensus.mat respectively.
  # Call to binary.vectors.lth.partition.function.
  # ************************************************************************** #
  # if M.run = 1 (the default) then it returns the binary matrix A(l) corresponding to
  # the input partition - rows and columns are ordered 'in some way'
  # ########################################################################## #
  start.time <- Sys.time()
  if(!all.equal(dim(final.clme.per.run.mat3), c(ngen3, M.run3)))
    stop("final.clme.per.run.mat3 should have ngen3 rows and M.run3 columns.\n")
  if(length(genelab.m3)!=ngen3 && !(is.null(genelab.m3)))
    stop("genelab.m3 has to have length=ngen3.\n")
  # All the columns in final.clme.per.run.mat3 should have K3 different values. #
  # ########################################################################### #
  # ########################################################################### #
  result <- list()
  result$A.l.binary.hatP <- list() # return A(l) in a sublist
  B.consensus.mat <- matrix(0, nrow=ngen3, ncol=ngen3)
  rownames(B.consensus.mat) <- genelab.m3
  colnames(B.consensus.mat) <- genelab.m3
  v.binary.k.mat.list <- list() # [[l]]: binary vectors by rows indicating belonging to kth cluster
  V.mat.K.list <- list() # [[k]]: the K3 matrices here = v_k %*% t*v_k = V_k(l)
  for(l in 1:M.run3)
  {
    # ########################################################################## #
    # Compute A(l) with a_ij = 1 row and column genes are clustered together, 0 otherwise.
    # This is Reduce('+', sum.V.mat.K.list) below. Add up these matrices by l and the result is
    # B.consensus.mat = consensus matrix (b_ij = no. of times row and column
    # gene clustered together within the M.run3 partitions).
    # a) Get the K3 binary vectors v_k [i] = 1 ith gene is in kth clu., 0 otherwise,
    # using binary.vectors.lth.partition.function
    # b) calculate the matrices V_k(l) = v_k %*% t(v_k) = V.mat.K.list[[l]][[k]]
    # c) add V.mat.K.list[[l]][[k]] up over k and this is A(l) = A.l.binary.mat.hatP.list[[l]] # https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
    # Call to binary.vectors.lth.partition.function
    # Using K3, ngen3, final.clme.per.run.mat, m.gen.clu3, cum.m.gen, m, no.correct.ones
    # ########################################################################### #
    v.binary.k.mat.list[[l]] <- binary.vectors.lth.partition.function(K2=K3,
                                                                      ngen2=ngen3, hatP.l2=final.clme.per.run.mat3[, l])
    # we will have K3 ngen3 x ngen3 matrices for each l
    sum.V.mat.K.list <- list() # [[k]] the K3 matrices.
    for(k in 1:K3)
    {
      v.k.temp <- v.binary.k.mat.list[[l]][k, ] # == v_k for this l
      sum.V.mat.K.list[[k]] <- v.k.temp %*% t(v.k.temp) # colnames are kept but rownames are not
      rownames(sum.V.mat.K.list[[k]]) <- colnames(sum.V.mat.K.list[[k]]) # dim(sum.V.mat.K.list[[k]])=ngen3 x ngen3
    }
    A.l.binary.hatP.temp <- Reduce('+', sum.V.mat.K.list)
    B.consensus.mat <- B.consensus.mat + A.l.binary.hatP.temp # 'lth part' of the consensus matrix
    result$A.l.binary.hatP[[l]] <- A.l.binary.hatP.temp
  } #end for the lth run
  result$B.consensus.counts <- B.consensus.mat
  result$B.consensus.mat <- B.consensus.mat/M.run3
  end.time <- Sys.time()
  result$timing.min <- (end.time-start.time) # /60 # in minutes
  return(result)
}








iterative.clustering.function <- function(seed1=seed3, ngen5=ngen, genelab.m5=genelab.m, K5=K, m.gen.clu5=m.gen.clu, df.temp5=3,
 m.it5=m.it, Ytemp1=y.summer.m.transf.loess, x1temp1=x1.spring.m.transf.loess, x2temp1=x2.autumn.m.transf.loess,
 x3temp1=x3.winter.m.transf.loess, timep5=timep, l=1, nT=24, norm="l1", bspl.bas.no1=12, bspl.bas.ord1=4,
 n.knots=parent.frame()$n.knots, seed2=parent.frame()$seed2, ari.per.iteration=0, true.partit.for.ari=NA)
{
  # ############################################################################## #
  # Functional Regression Mixture Modelling - Run l of the iterative clustering algorithm.
  # Starting with a random partition of K5 clusters, fit an FDboost model with Ytemp1 versus x1temp1, x2temp2, x3temp3 by cluster. Get
  # the predicted values for all observations and all models. Update the partition that minimises the residual per observation.
  # Repeat the procedure i.e. fit an FDboost model by cluster for the updated clusters, get the predicted values, etc. either
  # until convergence (=find the same partition) or until reaching the maximum no. of iterations (m.it5)
  # Input: seed1: seed to draw a starting random partition P_0, ngen5= no. of observations (= genes),
  # genelab.m5: character vector of length=ngen5 with the observations (genes) labels, K5=number of clusters,
  # m.gen.clu5= vector of K5 natural numbers that add up to ngen5 (= the sample sizes of the clusters of the starting random partition),
  # ############################################################################## #
  # df.temp5: the degrees of freedom used in the FDboost fit. It is essential that they are read indeed in the main environment though.
  # ############################################################################## #
  # m.it5: the (maximum) number of iterations, Ytemp1 = a matrix with dimensions ngen5 x nT indicating a functional response variable
  # x1temp1, x2temp1, x3temp1 = ngen5 x nT matrices with the values of functional explanatory variables
  # timep5: vector with length = nT indicating the time labels
  # l = number of run in the whole FRMM algorithm procedure.
  # nT: number of time points; norm: type of norm used in FRMM. It can be l1, l2 or l2cont.
  # bspl.bas.no1, bspl.bas.ord1: when norm="l2cont", the number and order of B-splines used in the smooth of the data respectively. Otherwise not used
  # n.knots: number of knots of the FDboost calls; seed2: seed set before the FDboost call
  # ari.per.iteration: 1 computes the %ARI between the resulting partition in the l-th iteration and true.partit.for.ari, 0 otherwise;
  # if 1, then it saves it in $ari.iteration.mat (4 x m.it5 matrix with the %ARI, %RI, %TPR, %TNR)
  # true.partit.for.ari: if ari.per.iteration=0 it should be an NA. Otherwise, it is a vector of length=
  # ngen5 with the cluster membership observed variable of the partition we want to compare to in order to get
  # the ARI for the partition in each iteration.
  # #################################################################################
  # creating P.0[i], P.0.lab[[k]], cum.m.gen[k], m.test.temp, distances.obs.fitted.mat,
  # clumem.m.mat, m.train.temp, m.train.mat, functional variables: Ytemp1, x1temp1, x2temp1, x3temp1,
  # dataset.k.test, mtemp, dataset.k.train.x1x2x3.loess.list, select.these.k
  # Calls to binary.matrix.truth.function, prop.correct.pairs.lth.partition.and.truth.function
  # Output: list res.P.list with final convergent partition P = final.clme.per.run.mat[, l]
  # and total sum of the norms of the residuals from fitting the FDboost models with P
  # Algorithm set with m.it iterations. If (and only if) it did not reach convergence with this no. of iterations,
  # then the output of res.P.list$achieved.convergence.lth.run will be NA. # ###
  # ######## #
  if(length(genelab.m5)!=ngen5)
    stop("Length of the genelab.m5 vector differs from ngen5. \n")
  if((sum(dim(Ytemp1)==c(ngen5, nT))!=2) | (sum(dim(x1temp1)!=c(ngen5, nT))==2) |
     (sum(dim(x2temp1)==c(ngen5, nT))!=2) | (sum(dim(x3temp1)==c(ngen5, nT))!=2))
    stop("Dimensions of all the functional variables have to be ngen5 x nT. \n")
  if(length(timep5)!= nT)
    stop("Length of timep5, with the timing observations, has to be nT. Check\n")
  if(sum(rownames(Ytemp1)== genelab.m5)!=ngen5)
    stop("Gene labels in Ytemp1 do not coincide with genelab.m5. \n")
  if(!(norm %in% c("l1", "l2", "l2cont")) ) stop("Norm should be l1, l2 or l2cont.")
  if(ari.per.iteration!=0 && ari.per.iteration!=1) # ADDED 11-4-20 # ### #
    stop("ari.per.iteration has to be either 1 (compute %ARI per iteration, in the run), or 0 (do not).\n")
  if((ari.per.iteration==1) && length(true.partit.for.ari) != ngen5) # ADDED 11-4-20; 1st bit reupdated 13-5-20 # ### #
    stop("true.partit.for.ari has to be a partition in order to calculate the %ARI when compared with the resulting FRMM
         partition in each iteration.\n")
  if((ari.per.iteration==1) && length(table(true.partit.for.ari))!=K5)
    stop("true.partit.for.ari has to be a partition with K5 clusters.\n")
  if((!is.vector(m.gen.clu5)) | (length(m.gen.clu5)!=K5) | (sum(m.gen.clu5)!=ngen5) )
    stop("m.gen.clu5 has to be a vector of natural numbers adding up to ngen5 and of length = K5. It indicates the sample sizes of the starting random partition.\n")
  stopifnot( all(m.gen.clu5 == floor(m.gen.clu5)) ) # all entries in m.gen.clu5 have to be integer
  # ##################################################################################### #
  cum.m.gen <- cumsum(m.gen.clu5)
  m.train.mat.list <- list()
  res.P.list <- list() # output: final convergent partition and total sum of L.-norms of residuals
  timing.min <- NA # timing in minutes
  achieved.convergence.lth.run <- NA # iteration number when it reaches convergence
  if(ari.per.iteration==1)
  { # save the %ARI, %RI, etc. of the i-th iteration then #
    ari.iteration.mat <- matrix(NA, nrow=4, ncol=m.it5)
    rownames(ari.iteration.mat) <- c("%ARI", "%RI", "%TPR", "%TNR")
    colnames(ari.iteration.mat) <- paste0("it", 1:m.it5)
  } # end if
  start.time <- Sys.time()
  # ################################################################## #
  # ################################################################## #
  # ################################################################## #
  set.seed(seed1)
  temp <- sample(1:ngen5)
  P.0 <- rep(NA, ngen5)
  names(P.0) <- genelab.m5 # "Ahg311263" ...
  P.0.lab <- list()
  k <- 1
  P.0.lab[[k]] <- genelab.m5[temp[1:m.gen.clu5[k]]]
  P.0.lab[[k]][1:3] # "Ahg934533" "Ahg478938" ...
  P.0[P.0.lab[[k]] ] <- k
  for(k in 2:K5) # k <- 2
  {
    # for l-th run,create initial partition (m.it5=1) #
    P.0.lab[[k]] <- genelab.m5[temp[(cum.m.gen[k-1]+1):cum.m.gen[k]]]
    P.0[P.0.lab[[k]] ] <- k # assign genes in kth cluster (updating, lth run)
    # P.0[1:5] # Ahg311263, 2, Ahg311286, 1, ...
  }
  m.test.temp <- ngen5
  distances.obs.fitted.mat <- matrix(NA, nrow=ngen5, ncol=K5) # L.-norms of residuals #
  rownames(distances.obs.fitted.mat) <- genelab.m5
  colnames(distances.obs.fitted.mat) <- paste0("cl", 1:K5)
  clumem.m.mat <- matrix(NA, nrow=ngen5, ncol=m.it5)
  rownames(clumem.m.mat) <- genelab.m5
  colnames(clumem.m.mat) <- paste0("it", 1:m.it5)
  m.train.temp <- m.gen.clu5
  m.train.mat <- matrix(NA, nrow=m.it5, ncol=K5)
  colnames(m.train.mat) <- paste0("cl", 1:K5)
  rownames(m.train.mat) <- paste0("it", 1:m.it5)
  m.train.mat[1, ] <- m.train.temp # no. genes in clu. in jth iteration #
  genelab.clu.j.list <- P.0.lab # [[k]]: gene labels in kth clu., jth iteration #
  dataset.k.test <- list()
  dataset.k.test$Ytemp <- Ytemp1
  dataset.k.test$timep <- timep5
  dataset.k.test$s.timep <- timep5
  dataset.k.test$genelabs.k.clu <- rownames(Ytemp1) # "Ahg311263" "Ahg311286" ...
  clumem.m.mat[, 1] <- P.0
  # ############################################################### #
  for(j in 1:m.it5) # Update the resulting partition in jth iteration # j <- 1 #
  {
    cat("# ############################################################# #\n")
    cat("l=", l, " j = ", j, "Sample sizes in this iteration: ", m.train.temp, "\n")
    for(k in 1:K5) # k <- 1 #
    { # centre explan. vars, fit FDboost model, jth iteration, kth clu., subtract mean
      # in testing set, predicted values using predicted.temp.testing[itest, k], L.-norm
      # of residuals by gene, # using m.train.temp, genelab.clu.j.list[[k]]: gene labels kth clu.
      mtemp <- m.train.temp[k]
      # List with the variables in the kth cluster in the jth iteration (start)
      dataset.k.train.x1x2x3.loess.list <- list()
      temp <- Ytemp1[genelab.clu.j.list[[k]], ] # mtemp x 24 time points
      dataset.k.train.x1x2x3.loess.list$Ytemp <- temp  # e.g. 1793 x 24 time points# response for genes in kth clu.
      temp1 <- x1temp1[genelab.clu.j.list[[k]], ] # 1793 x 24 x1.spring.m.transf.loess, etc.
      dataset.k.train.x1x2x3.loess.list$x1temp.centred <- scale(temp1, center=TRUE, scale=FALSE)
      temp2 <- x2temp1[genelab.clu.j.list[[k]], ] # x2.autumn.m.transf.loess
      #dataset.k.train.x1x2x3.loess.list$x2temp <- temp2 # training x2
      dataset.k.train.x1x2x3.loess.list$x2temp.centred <- scale(temp2, center=TRUE, scale=FALSE)
      temp3 <- x3temp1[genelab.clu.j.list[[k]], ] # x3.winter.m.transf.loess
      #dataset.k.train.x1x2x3.loess.list$x3temp <- temp3 # training x3
      dataset.k.train.x1x2x3.loess.list$x3temp.centred <- scale(temp3, center=TRUE, scale=FALSE)
      dataset.k.train.x1x2x3.loess.list$timep <- timep5
      dataset.k.train.x1x2x3.loess.list$s.timep <- timep5
      dataset.k.train.x1x2x3.loess.list$df.temp5 <- df.temp5
      dataset.k.train.x1x2x3.loess.list$genelabs.k.clu <- genelab.clu.j.list[[k]] # 1793
      # Fit FDboost model with the 3 functional main effects
      # integrals over all the range of time points: from 1 to 24
      set.seed(seed2) # E(Yi|Xi)=beta0(t)+int beta2(s, t)*x2(s)*ds+int beta3(s, t)*x3(s)*ds # usually seed2=5 say
      a1 <- try(fof.functional.temp <- FDboost::FDboost( Ytemp ~ 1 +
         bsignal(x=x1temp.centred, s=s.timep, knots=n.knots, df = df.temp5) +
         bsignal(x=x2temp.centred, s=s.timep, knots=n.knots, df = df.temp5) +
         bsignal(x=x3temp.centred, s=s.timep, knots=n.knots, df = df.temp5),
       timeformula = ~ bbs(timep, knots=n.knots, df = df.temp5),
       data = dataset.k.train.x1x2x3.loess.list, control = boost_control(mstop = 300)) )
      if(class(a1)[1]=="try-error")
      { # https://stat.ethz.ch/R-manual/R-devel/library/base/html/try.html
        cat("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n")
        cat("l = ", l, ", j = ", j, "th iteration, k = ", k, "th cluster or model. FDboost command with\n",
            " no.knots = ", n.knots, " reports an error. It is \n", a1[1])
        cat("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n")
        n.knots2 <- 7
        b <- try(fof.functional.temp <- FDboost::FDboost( Ytemp ~ 1 +
         bsignal(x=x1temp.centred, s=s.timep, knots=n.knots2, df = df.temp5) +
          bsignal(x=x2temp.centred, s=s.timep, knots=n.knots2, df = df.temp5) +
          bsignal(x=x3temp.centred, s=s.timep, knots=n.knots2, df = df.temp5),
        timeformula = ~ bbs(timep, knots=n.knots2, df = df.temp5),
        data = dataset.k.train.x1x2x3.loess.list,
         control = boost_control(mstop = 300)) )
        if(class(b)[1]=="try-error")
          stop("Investigate this. With 7 knots, there is an error too.\n")
      } # end if class of (a1)[1]==try-error #
      a <- matrix(rep(colMeans(temp1), m.test.temp), nrow=m.test.temp,
                  byrow=TRUE) # 5378 x 24
      dataset.k.test$x1temp.centred <- x1temp1-a # centred x1 with mean from training observations #
      a <- matrix(rep(colMeans(temp2), m.test.temp), nrow=m.test.temp,
                  byrow=TRUE) # by time point and variable #
      dataset.k.test$x2temp.centred <- x2temp1-a
      a <- matrix(rep(colMeans(temp3), m.test.temp), nrow=m.test.temp,
                  byrow=TRUE)
      dataset.k.test$x3temp.centred <- x3temp1-a
      # ############################################################# #
      predicted.temp.testing.k <- predict(fof.functional.temp, newdata=dataset.k.test) # ngen5 x nT with the fitted Y values in all genes #
      # L1, L2 or L2cont distance between the vector of observed (over time points) and
      # fitted by gene, kth clu., jth it. #
      if(norm=="l1")
      {
        dist.l1.mat.temp <- abs(Ytemp1-predicted.temp.testing.k)
        distances.obs.fitted.mat[, k] <- rowSums(dist.l1.mat.temp) # over time points #
      } else
        if(norm=="l2") # although named 'l1', it is using indeed the L2 distance #
        {
          dist.l1.mat.temp <- (Ytemp1-predicted.temp.testing.k)^2
          distances.obs.fitted.mat[, k] <- sqrt(rowSums(dist.l1.mat.temp)) # over time points #
        } else # L2 continuous i.e. doing a smooth of Y and another smooth of ^Y #
          if(norm=="l2cont") # bspl.bas.no1=12, bspl.bas.ord1=4
          {
            # Calculate a smooth of Y, another smooth of ^Y, using bspl.bas.no1, bsp.bas.ord1, subtract them #
            # and get the L2-continuous norm of the resulting 'fd' smooth object or "vector"
            # using 'norm.fd' in the R fda.usc package
            hourbasis24b <- create.bspline.basis(rangeval=c(0, nT), nbasis=bspl.bas.no1, norder=bspl.bas.ord1)
            tempfd1 <- smooth.basis(timep5, t(Ytemp1), hourbasis24b)$fd # smooth of Y #
            names(tempfd1$fdnames) <- c("time", "gene labels", "gene expression")
            tempfd2 <- smooth.basis(timep5, t(predicted.temp.testing.k), hourbasis24b)$fd # smooth of ^Y #
            names(tempfd2$fdnames) <- c("time", "gene labels", "gene expression")
            resid.temp.fd <- tempfd1 - tempfd2
            a <- fda.usc::norm.fd(resid.temp.fd)
            if(length(diag(a))!=ngen5) stop("Dimensions of the L2-continuous distance matrix are not ngen5 x ngen5. Check\n")
            if(any(is.na(diag(a))) ) stop("norm.fd command in fda.usc package produced some L2-continuous missing norm. Check\n")
            distances.obs.fitted.mat[, k] <- diag(a) # https://www.rdocumentation.org/packages/fda.usc/versions/1.5.0/topics/norm.fdata
            # ################################################################################################### #
            # this is an ngen x ngen matrix with all the L2-cont. distances. There are quite a lot of NaN
            # but we are only interested in the elements in the diagonal, and none of those has an NaN.
            # ################################################################################################### #
          } else
            stop("Norm can only be l1 (default), l2, or l2cont.\n")
      # Output of k loop: distances.obs.fitted.mat[, k] # ########### #
    } # end of kth model in jth iteration #
    # ################################# #
    # If this iteration produces a partition that indicates convergence then we want to keep the L. norms
    # in the fittings of the models for all clusters
    # Calculate P_{j+1} using distances.obs.fitted.mat: clumem.m.mat[, j+1]
    # update genelab.clu.j.list[[k]]
    # if j=m.it5 then we do not have a future partition, m.train.temp, gene lab in training set,
    # so we finish here. We do not know whether j=m.it5 exactly reached convergence.
    if(j==m.it5)
    {
      cat("Attention: it did not reach convergence up to and including the m.it5-1 iteration.")
      break # this breaks produces to come out from the jth iteration for loop. Consequently, the value
      # of achieved.convergence.lth.run will be NA. Vice versa, when achieved...=NA, it implies
      # necessarily that it did not reach convergence. #
    } # command below: if there are 'ties', it picks one of the elements of the tie at random #
    clumem.m.mat[, j+1] <- max.col(-distances.obs.fitted.mat) # col. where min. lies by row #
    # partition in j+1 is already computed.
    # call to prop.correct.pairs.lth.partition.and.truth.function #
    if(ari.per.iteration==1) # if requested to compute the %ARI, %RI, etc. per iteration. I.e. compare
    { # clumem.m.mat[, j+1] with true.partit.for.ari, and keep it in ari.iteration.mat[, j+1].
      # Need subject.lab.ordered.true.partition. CALLS TO BINARY.MATRIX.TRUTH.FUNCTION,
      # PROP.CORRECT.PAIRS.LTH.PARTITION.AND.TRUTH.FUNCTION
      binary.mat.true.partitiontemp <- binary.matrix.truth.function(K1=K5, ngen1=ngen5, true.tildetildeP1=true.partit.for.ari)
      subject.lab.ordered.true.partitiontemp <- rownames(binary.mat.true.partitiontemp)
      restemp <- prop.correct.pairs.lth.partition.and.truth.function(M1=1, K6=K5, ngen6=ngen5, final.clme.per.run.mat1=clumem.m.mat[, j+1, drop=FALSE],
       genelab.m6=genelab.m5, true.tildetildeP6=true.partit.for.ari, m.gen.clu6=as.numeric(table(true.partit.for.ari)),
       genelab.ordered.true.partition6=subject.lab.ordered.true.partitiontemp, achieved.convergence.run.vec1=1)
      ari.iteration.mat[, j+1] <- restemp$prop.correct.ones.zeros.acc.FP.FN.mat.ari[1:4, , drop=FALSE]*100
    } # End if ari.per.iteration==1 # ################################################################################# #
    # clumem.m.mat[, j+1]==clumem.m.mat[, j] <=> convergence at jth iteration so all future
    # partitions will be equal to the jth one. #
    if(sum(abs(clumem.m.mat[, j+1]-clumem.m.mat[, j]))==0)
    { # copy clumem.m.mat[, j+1] in the future and finish
      achieved.convergence.lth.run <- j
      if(j==m.it5-1) break
      clumem.m.mat[, (j+2):m.it5] <- matrix(clumem.m.mat[, j+1], nrow=ngen5,
                                            ncol=m.it5-j-2+1, byrow=FALSE)
      cat("l = ", l, " End of ", j, "th iteration. Convergence achieved - - - - - - - - - - - - - - \n")
      break   # break of the j for loop #
    } # end if convergence found in jth iteration; save later final convergent partition #
    # If no convergence was achieved, update m.train.temp, m.train.mat[j+1, ],  # updated partition genelab.clu.j.list[[k]]
    if(length(table(clumem.m.mat[, j+1]))!=K5) stop("Attention, maybe (at least) one cluster is empty.")
    if(sum(as.numeric(names(table(clumem.m.mat[, j+1])) ) != c(1:K5)) )
      stop("Recheck the partition in this iteration. Are there K5 nonempty clusters?")
    m.train.temp <- as.vector(table(clumem.m.mat[, j+1]))
    m.train.mat[j+1, ] <- m.train.temp
    cat("l=", l, ", j = ", j, "Sample sizes of partition in next iteration: ",
        m.train.temp, "\n")
    # Gene labels of clusters in (j+1)th iteration: #
    genelab.clu.j.list <- list()
    for(k in 1:K5) # k <- 1
    {
      select.these.k <- clumem.m.mat[, j+1]==k
      genelab.clu.j.list[[k]] <- names(select.these.k)[select.these.k=="TRUE"]
    }
    # then get the cluster membership variable for the next iteration until the last iteration #
    cat("l=", l, ", end of ", j, "th iteration - - - - - - - - - - - - - - \n")
    next # go to the next iteration until convergence #
    # Output of this loop: clumem.m.mat #
  } # end of jth iteration (lth run) # ###################################################### #
  # Computer will arrive here only if it found that the iteration is convergent (expectably).
  # Save final convergent partition, total sum of L.- distances of residuals from the fittings
  # with this partition (done in jth loop before breaking) as well as m.train.mat[j, k]=no. genes in kth clu. in final convergent partition.
  m.train.mat.list[[l]] <- m.train.mat[1:j, ]
  res.P.list$final.clme.l <- clumem.m.mat[, m.it5, drop=FALSE] # sometimes matrix has a single column #
  res.P.list$m.i.convergent.clme.l <- m.train.mat[1:j, ]
  res.P.list$total.norm.resid.l <- sum(distances.obs.fitted.mat)
  res.P.list$achieved.convergence.lth.run <- achieved.convergence.lth.run
  # List with gene labels in final convergent partition: #
  final.clme.labels.list <- list()
  temp <- clumem.m.mat[, m.it5]
  for(k in 1:K5)
  {
    final.clme.labels.list[[k]] <- names(temp)[temp==k]
  }
  res.P.list$final.clme.l.lab.list <- final.clme.labels.list
  # ################## IF WE WANT THE FUNCTION TO COMPUTE THE %ARI PER ITERATION, IN EACH RUN: # ############ #
  if(ari.per.iteration==1)
  { # IF REQUESTED TO COMPUTE THE %ARI, %RI, ETC. PER ITERATION, REPORT IT IN THE OBJECT BELOW: # ##### #
   res.P.list$ari.tpr.per.iteration.mat <- ari.iteration.mat # 4 x m.it5 matrix # each column has the %ARI, %RI, %TPR, %TNR in the iteration #
  } # End if ari.per.iteration==1 # ############################################################################################### #
  end.time <- Sys.time()
  res.P.list$timing.min <- (end.time-start.time) # in minutes #
  return(res.P.list)
  # End of iterative.clustering.function #
}






prop.correct.pairs.lth.partition.and.truth.function <- function(M1=1, K6=K, ngen6=ngen,
final.clme.per.run.mat1=final.clme.per.run.mat, genelab.m6=NULL, true.tildetildeP6=true.tildetildeP,
 m.gen.clu6=as.vector(table(true.tildetildeP)), genelab.ordered.true.partition6=NULL, # genelab.ordered.true.partition
 achieved.convergence.run.vec1=1) # last argument: added 10-1-20
{ # Adjusted Rand Index, added 29-11-19, $ari.vec # ############ #
  # genelab.ordered...=NULL: names(truetildetildeP6) should be NULL too, as well as probably rownames(final.clme...)
  # but even if they were not NULL,the function might work correctly. 5-6-20
  # modified so that it deals correctly when possibly some run is nonconvergent (indicated by an NA in achieved...):
  # when studying the L2 continuous norm, we saw some nonconvergent runs. These are not detected in
  # final.clme.per.run.mat, as when happening, somehow there is a value therein. Thus, this
  # function needs to have an extra argument: achieved.convergence.run.vec1; when =NA, it indicates, without
  # any possible confusion, that the run was not convergent. Added this argument in the function, 10-1-20
  # It is a vector of length=M.run. The total of NAs in it indicates the no. of nonconvergent runs # ##### #
  # Note 2: this code also works when the partition in final.clme.per.run.mat1[, l] has less clusters than K6
  # (checked: 7-2-20, 65), Alt. 4 or 5, with K=12 and 2, 5 clu. resp.). I did not "check" the case when it has more
  # clusters than K6 but it might happen that it works too in that case.
  # ############################################################# #
  # ############################################################# #
  # m.gen.clu6 == cardinalities of the clusters in true.tildetildeP6, # ############################# #
  # Given (essentially) a "true" partition Pt=true.tildetildeP6 and the lth estimated partition hatP.l[, l]=final.clme.per.run.mat1[, l]
  # In the space of the combinations of observations (=genes) taken 2 at a time (so cardinality = n choose 2),
  # given the events: A.t = the pair clustered together in Pt; A.l=the pair clustered together in P.l,
  # CALCULATE THE ADJUSTED RAND INDEX (HUBERT AND ARABIE, 1985), USING THE R CLUES PACKAGE - adjustedRandIndex(P.truth, ^P.l);
  # PROPORTIONS OF CORRECT CLASSIFIED PAIRS (=CORRECT ONES + CORRECT ZEROS = ACCURACY = RAND INDEX)
  # = (TRUE 'POSITIVES' + TRUE 'NEGATIVES') / (N CHOOSE 2);  THE TRUE POSITIVE RATE = No. (Pt=1, P.l=1)/No. (Pt=1) AND
  # THE TRUE NEGATIVE RATE
  # TRUE NEGATIVES = No. (Pt=0, P.l=0)/No. (Pt=0), FALSE POSITIVES = No. (Pt=0, P.l=1)/No. (Pt=0),
  # FALSE NEGATIVES = No. (Pt=1, P.l=0)/No. (Pt=1)
  # OF AN ESTIMATED PARTITION ^P_l WHEN COMPARED WITH Pt
  # Given the lth estimated partition final.clme.per.run.mat1[, l] and a true partition true.tildetildeP6,
  # find A(l), a_ij(l)=1 if ith and jth observations (genes) clustered together in
  # final.clme.per.run.mat1[, l], 0 otherwise. This is Reduce('+', sum.V.mat.K.list) below.
  # Note that rows and columns in A(l) are not ordered as in the truth=binary.mat.true.partition.
  # Add up these matrices by l, divide over M1 and the result is
  # B.consensus.mat = consensus matrix (b_ij = no. of times row and column gene clustered together).
  # ACCURACY (=RI): compare A.sor(l) (after ORDERING its rows and columns) with binary.mat.true.partition
  # a) Get the K6 binary vectors v_k [i] = 1 ith gene is in kth clu., 0 otherwise,
  # using binary.vectors.lth.partition.function
  # b) calculate the matrices V_k(l) = v_k %*% t(v_k) = V.mat.K.list[[l]][[k]]
  # c) sum V.mat.K.list[[l]][[k]] over k and this is A(l) = A.l.binary.mat.hatP.list[[l]]
  # d) Reorder rows and columns in A.l.binary.mat.hatP.list[[l]] as in
  # binary.mat.true.partition: A.l.binary.mat.hatP.list.sor[[l]]
  # e) Calculate no.correct.ones[l]=true 'positives' (we define a positive as a pair that clustered
  # together in hatP.l[, l]), no.correct.zeros[l] = true 'negatives' to get later the accuracy
  # of ^P(l) when compared to Pt with the pairs of genes that are clustered together or not.
  # Call to binary.vectors.lth.partition.function # ############################## #
  # Using K6, ngen6, final.clme.per.run.mat1, m.gen.clu6, cum.m.gen, m, no.correct.ones
  # Input: mainly the true partition ~~P, the estimated partition hatP.l, ...
  # ---------------------------------------------------------------------------------------------- #
  # Output: list with 10 components (from 7-6-20, 10-1-20 on): $prop.correct.ones.zeros.acc.FP.FN.mat.ari (6xM1 matrix), $B.consensus.mat[i, j] (only includes results for
  # the convergent partitions), $M.effective=M1 - possible no. of nonconverg. partitions (i.e. latter possibly !=0), $no.nonconvergent.partitions (vector length=M1),
  # plus: $no.Acc.TP.TN.FP.FN.mat, $ari.vec[l], $true.tildetildeP, $hatP.l, (prop.correct.ones.zeros.acc.FP.FN.mat left out now)
  # $PercAccuracy, $timing.vec. The no. NAs in ari.vec, columns of prop., etc. should be exactly = $no.nonconvergent.partitions
  # ################################################################## #
  # genelab.m6 should be a label vector exactly as the names in true.tildetildeP, in the same order; or both NULL (added 5-6-20)
  # final.clme.per.run.mat1 should have row names, and these should be exactly like genelab.m6.
  # the cardinalities of the clusters in true.tildetildeP6 should be exactly m.gen.clu6. ********* #
  # ################################################################## #
  if(!is.matrix(final.clme.per.run.mat1) | !(all.equal(dim(final.clme.per.run.mat1),
                                                       c(ngen6, M1))) )
    stop("final.clme.per.run.mat1 should be a matrix and it should have dimensions = ngen6 x M1.\n")
  if(length(m.gen.clu6)!=K6)
    stop("Length of m.gen.clu6 should be K6 = the no. of clusters.\n")
  # added 5-6-20: possibility that genelab.m6 is null, as well as names(true.tildetildeP6) and rownames(final...mat1)
  if(!is.null(genelab.m6))
  {
    if(!all.equal(rownames(final.clme.per.run.mat1), genelab.m6))
      stop("The row names of final.clme.per.run.mat1, with the genes, should coincide with genelab.m6.\n")
  }
  if(sum(abs(as.numeric(table(true.tildetildeP6)) -m.gen.clu6) )) # added 9-12-19 #
    stop("The no. of elements in each cluster in true.tildetildeP6 should be equal to m.gen.clu6 and in the same order.\n")
  if(!is.vector(achieved.convergence.run.vec1) | length(achieved.convergence.run.vec1)!=M1) # 10-1-20
    stop("achieved.convergence.run.vec1 should be a vector of length=M1 (and if NA it means nonconvergent run). Added 10-1-20\n")
  # modified below, 5-6-20: for an R package, one should write e.g. clues::adjustedRand(... etc.). NEVER CALL TO A LIBRARY
  # https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html
  # if(!require(clues)) install.packages("clues"); library(clues)
  # ############################################################################# #
  res <- list()
  res$true.tildetildeP <- true.tildetildeP6
  res$hatP.l <- final.clme.per.run.mat1
  m <- m.gen.clu6
  # ################################################################################################ #
  # Accuracy == no. correct 1s + no. correct 0s; TP=correct 1s; TN=correct 0s; FN=1 in the truth and
  # 0 in the lth partition; FP=0 in the truth and 1 in the lth partition. Adjusted RI=(RI-E(RI))/(1-E(RI))
  # see Hubert and Arabie's paper, R clues package with adjustedRand
  # ################################################################################################ #
  no.correct.ones.zeros.acc.FP.FN.mat <- matrix(0, nrow=5, ncol=M1)
  rownames(no.correct.ones.zeros.acc.FP.FN.mat) <- c("Acc", "TP", "TN", "FP", "FN")
  colnames(no.correct.ones.zeros.acc.FP.FN.mat) <- paste0("q", 1:M1)
  prop.correct.ones.zeros.acc.FP.FN.mat <- matrix(0, nrow=5, ncol=M1)
  rownames(prop.correct.ones.zeros.acc.FP.FN.mat) <- c("Acc", "TPR", "TNR", "FPR", "FNR")
  colnames(prop.correct.ones.zeros.acc.FP.FN.mat) <- paste0("q", 1:M1)
  #no.correct.ones <- rep(0, M1) # = true positives
  #names(no.correct.ones) <- paste0("q", 1:M1)
  #no.correct.zeros <- no.correct.ones
  cum.m.gen <- cumsum(m.gen.clu6)
  v.binary.k.mat.list <- list() # [[l]]: binary vectors by rows indicating the belonging to the kth clu.
  # (the order of the components of these vectors will typically differ from that in genelab.ordered.true.partition6)
  B.consensus.mat <- matrix(0, nrow=ngen6, ncol=ngen6) # DO NOT KEEP ALL THE MATRICES BY L; MUCH MORE EFFICIENT IN THIS OTHER WAY
  if(!is.null(genelab.m6)) # modified, 5-6-20
  {
    rownames(B.consensus.mat) <- genelab.m6 # sum over all the l's
    colnames(B.consensus.mat) <- genelab.m6
  } else
  {
    rownames(B.consensus.mat) <- 1:ngen6 # sum over all the l's
    colnames(B.consensus.mat) <- 1:ngen6
  } # end of modification 5-6-20
  timing.vec <- rep(NA, M1)
  names(timing.vec) <- paste0("l", 1:M1)
  ari.vec <- rep(NA, M1) # using clues package. Added 29-11-19. NA will mean nonconvergent partition
  names(ari.vec) <-  paste0("l", 1:M1)
  for(l in 1:M1) # l <- 1
  { # https://davetang.org/muse/2017/09/21/adjusted-rand-index/ # ######################################## #
    # ##################################################################################################### #
    # This for loop computes the adjusted Rand index for the lth partition, the dummy binary vectors
    # v indicating (individual) cluster membership, by cluster, the matrix product of the square of those, ...,
    # no.correct.ones.zeros.acc.FP.FN.mat["TP", l], TN, FP, FN, B.consensus.mat although here its entries are
    # counts, not proportions, total.0s.in.truth; i.e. it fills inter alia ari.vec[l],
    # no.correct.ones.zeros.acc.FP.FN.mat[ , l] total.0s.in.truth #
    # if the lth partition is not convergent then ari.vec[l] will keep an NA (if and only if),
    # no.correct.ones.zeros.acc.FP.FN.mat[, l] will be NAs, and B.consensus.mat will not add anything for this run #
    # ###################################################################################################### #
    if(!is.na(achieved.convergence.run.vec1[l])) # added 10-1-20. I.e. if the partition was a convergent one
      ari.vec[l] <- clues::adjustedRand(true.tildetildeP6, final.clme.per.run.mat1[, l])["HA"] # numeric vector
    # https://www.r-bloggers.com/5-ways-to-measure-running-time-of-r-code/
    start.time <- Sys.time() # matrix with K clu. x ngen observations:
    v.binary.k.mat.list[[l]] <- binary.vectors.lth.partition.function(K2=K6,
                                                                      ngen2=ngen6, hatP.l2=final.clme.per.run.mat1[, l])
    sum.V.mat.K.list <- list() # [[k]] the K6 matrices. DO NOT DO IT [[l]]: R won't have enough memory
    #
    for(k in 1:K6) # k <- 1
    {
      v.k.temp <- v.binary.k.mat.list[[l]][k, ] # == v_k for this l
      sum.V.mat.K.list[[k]] <- v.k.temp %*% t(v.k.temp) # colnames are kept but rownames are not, ngen x ngen
      if(!is.null(colnames(sum.V.mat.K.list[[k]]))) # ADDED 5-6-20 #
        rownames(sum.V.mat.K.list[[k]]) <- colnames(sum.V.mat.K.list[[k]]) else
        {
          colnames(sum.V.mat.K.list[[k]]) <- 1:ngen6
          rownames(sum.V.mat.K.list[[k]]) <- 1:ngen6
        } # End of added 5-6-20 #
      # dim(sum.V.mat.K.list[[k]]) # e.g. 5378 5378 as expected # Ahg311263 Ahg311286 ...
    } # end for the kth cluster
    # c) add up V.mat.K.list[[l]][[k]]
    # https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
    # binary matrix where 1 = if row and column genes clustered together in the lth partition:
    # If the lth partition was not convergent then do not update B.consensus.mat. Otherwise, add the 1s =
    # pairs of observations that clustered together, indicated by A.l.binary.hatP.temp
    A.l.binary.hatP.temp <- Reduce('+', sum.V.mat.K.list)
    if(!is.na(achieved.convergence.run.vec1[l])) B.consensus.mat <- B.consensus.mat + A.l.binary.hatP.temp # lth part of the consensus matrix
    # ####################################################### #
    # ####################################################### #
    # ACCURACY=PROPORTION OF CORRECT 1S AND 0S OF THE LTH ESTIMATED PARTITION #
    # = (no. true positives + no. true negatives)/(n choose 2) ############# #
    # ####################################################### #
    # ####################################################### #
    # Order A.l.binary.hatP.temp as the binary matrix with the truth
    # 5-6-20: I AM NOT SURE WHETHER IT COULD EVER HAPPEN THAT THE LABELS IN TRUE.TILDETILDEP6 AND THOSE IN
    # FINAL.CLME...MAT1 WERE NOT IN THE SAME ORDER. I AM GOING TO CONSIDER THAT IF
    # GENELAB.ORDERED.TRUE.PARTITION6 IS NULL I ASSUME THAT ALL THE PARTITIONS (TRUTH AND THE OTHER/S) ARE
    # IN THE SAME ORDER AND THE LABELS IN A.L.BINARY.HATp.TEMP == 1:ngen6 == ORDER OF LABELS IN TRUE.TILDETILDEP6
    if(!is.null(genelab.ordered.true.partition6)) # added 5-6-20
    {
      atemp <- A.l.binary.hatP.temp[, genelab.ordered.true.partition6]
      A.l.binary.hatP.temp.sor <- atemp[genelab.ordered.true.partition6, ]
    } else # added 5-6-20
    {
      A.l.binary.hatP.temp.sor <- A.l.binary.hatP.temp # i.e. assume that it is the same order as in true.tildetildeP
    } # end else from if above
    # Calculate the accuracy i.e. compare A.l.binary.hatP.temp.sor with true.tildetildeP6
    # Number of overlapping 1s and 0s between A.l.binary.hatP.temp.sor
    # and binary.mat.true.partition = no. pairs of genes found in the same cluster
    # and they were indeed in the same cluster in true.tildetildeP6, and
    # vice versa with the 0s.
    k <- 1 # l <- 1
    total.0s.in.truth <- m[k]*(ngen6-m[k]) # See Report2.pdf formula; evident
    if(any(is.na(achieved.convergence.run.vec1[l]))) # if the lth partition is nonconvergent indeed
    { # ADDED 10-1-20 # ############################################################ #
      # Fill in total.0s.in.truth and next case i.e. do not add anything in B.consensus.mat, and NAs in Acc, etc.
      # no.correct.ones.zeros.acc.FP.FN.mat[, l]. The matrix prop... below will also have an NA.
      # total.1s.in.truth, prop. are calculated outside this lth for loop. # ######################################### #
      for(k in 2:K6) {
        total.0s.in.truth <- total.0s.in.truth + m[k]*(ngen6-cum.m.gen[k]) } # Report2.pdf
      no.correct.ones.zeros.acc.FP.FN.mat[, l] <- NA # Acc, TP, etc., all NA
      # prop.                # #################################### #
      next # i.e. compute the total no. of 0s in the truth, put an NA in acc, TP in no.correct... and go on
    }
    k <- 1 # again, continue, for an lth convergent partition now #
    submatrix.temp <- A.l.binary.hatP.temp.sor[1:m[k], 1:m[k]]
    no.ones.temp <- (sum(submatrix.temp, na.rm=TRUE) - m[k])/2 # in this submatrix
    no.correct.ones.zeros.acc.FP.FN.mat["TP", l] <-
      no.correct.ones.zeros.acc.FP.FN.mat["TP", l] + no.ones.temp
    #no.correct.ones[l] <- no.correct.ones[l] + (sum(submatrix.temp, na.rm=TRUE) - m[k])/2
    no.correct.ones.zeros.acc.FP.FN.mat["FN", l] <-
      no.correct.ones.zeros.acc.FP.FN.mat["FN", l] + (m[k]^2-m[k])/2 - no.ones.temp
    # true negatives (0 in truth and 0 in lth partition) and FP (0 in truth and 1 in lth part.)
    submatrix.temp1 <- A.l.binary.hatP.temp.sor[1:m[k], (m[k]+1):ngen6]
    mtemp <- m[k]*(ngen6-m[k]) # 6427905
    no.correct.ones.zeros.acc.FP.FN.mat["TN", l] <-
      no.correct.ones.zeros.acc.FP.FN.mat["TN", l] + (mtemp - sum(submatrix.temp1))
    no.correct.ones.zeros.acc.FP.FN.mat["FP", l] <-
      no.correct.ones.zeros.acc.FP.FN.mat["FP", l] + sum(submatrix.temp1)
    # accuracy = prop. correct 1s and 0s = correctly classified pairs. Later
    #no.correct.zeros[l] <- no.correct.zeros[l] + (mtemp - sum(submatrix.temp1))
    for(k in 2:K6) # for the lth partition, which is convergent: #
    {
      # Compute binary square and rectangular submatrices to find overlapping 1s in the kth clu. in the truth (~~P) with
      # the kth clu. in the estimated partition (^P_l). Then add up all the overlapping
      # 1s and store their sum in no.correct.ones[l]. Similarly with no.correct.zeros[l] # ##
      total.0s.in.truth <- total.0s.in.truth + m[k]*(ngen6-cum.m.gen[k]) # Report2.pdf
      submatrix.temp <- A.l.binary.hatP.temp.sor[(cum.m.gen[k-1]+1):cum.m.gen[k],
                                                 (cum.m.gen[k-1]+1):cum.m.gen[k]]
      no.ones.temp <- (sum(submatrix.temp, na.rm=TRUE) - m[k])/2
      no.correct.ones.zeros.acc.FP.FN.mat["TP", l] <-
        no.correct.ones.zeros.acc.FP.FN.mat["TP", l] + no.ones.temp
      no.correct.ones.zeros.acc.FP.FN.mat["FN", l] <- no.correct.ones.zeros.acc.FP.FN.mat["FN", l] + (m[k]^2-m[k])/2 - no.ones.temp
      # no.correct.ones[l] <- no.correct.ones[l] + (sum(submatrix.temp, na.rm=TRUE)-m[k])/2
      # ## now with overlapping 0s:
      if(k==K6) break # rectangles with 0s in the truth:
      submatrix.temp2 <- A.l.binary.hatP.temp.sor[
        (cum.m.gen[k-1]+1):cum.m.gen[k], (cum.m.gen[k]+1):ngen6]
      mtemp <- (cum.m.gen[k] - cum.m.gen[k-1])*(ngen6-cum.m.gen[k]) # 3213056
      no.correct.ones.zeros.acc.FP.FN.mat["TN", l] <-
        no.correct.ones.zeros.acc.FP.FN.mat["TN", l] + (mtemp - sum(submatrix.temp2))
      no.correct.ones.zeros.acc.FP.FN.mat["FP", l] <-
        no.correct.ones.zeros.acc.FP.FN.mat["FP", l] + sum(submatrix.temp2)
      #no.correct.zeros[l] <- no.correct.zeros[l] + (mtemp - sum(submatrix.temp2))
      # example: if l2 and l4 were divergent, we would have (NAs in those columns and some numbers in the other columns)
      #         q1 q2    q3 q4    q5  q100
      # Acc  72791 NA     0 NA     0     0
      # TP   31996 NA 16359 NA 16237 15756
      # TN  113586 NA 56559 NA 57689 57122
      # FP   53080 NA 26774 NA 25644 26211
      # FN   50838 NA 25058 NA 25180 25661
    }
    end.time <- Sys.time()
    timing.vec[l] <- end.time-start.time
    # summary: if the lth partition was not convergent, then we did not add anything in consensus.mat and
    # moreover, no.correct.ones.zeros.acc.FP.FN.mat[, l] = NA for all its entries (Acc, TP, etc.).
    # total.1s.in.truth and prop.correct.ones.zeros.acc.FP.FN.mat are calculated below. Ari.vec[l] is NA too#
  } # end for the lth estimated partition #
  # ############################################## #
  # IF LTH PARTITION IS NOT CONVERGENT WE WANT THE PROP...[, l] TO BE NAs. As no.correct...[, l]=NAs
  # for that case, the commands below produce an NA in prop... too (checked).
  # The entry in the matrix will be NA so the result of the assignments below will be NA#
  # we do not need no.correct.ones.zeros either
  total.1s.in.truth <- sum(choose(m, 2)) # 4817792
  # added, 12-9-19, after meeting with Daphne: #
  prop.correct.ones.zeros.acc.FP.FN.mat["TPR", ] <-
    no.correct.ones.zeros.acc.FP.FN.mat["TP", ]/total.1s.in.truth
  prop.correct.ones.zeros.acc.FP.FN.mat["TNR", ] <-
    no.correct.ones.zeros.acc.FP.FN.mat["TN", ]/total.0s.in.truth
  prop.correct.ones.zeros.acc.FP.FN.mat["FNR", ] <-
    no.correct.ones.zeros.acc.FP.FN.mat["FN", ]/total.1s.in.truth
  prop.correct.ones.zeros.acc.FP.FN.mat["FPR", ] <-
    no.correct.ones.zeros.acc.FP.FN.mat["FP", ]/total.0s.in.truth
  no.correct.ones.zeros <- no.correct.ones.zeros.acc.FP.FN.mat[c("TP", "TN"), , drop=FALSE]
  prop.correct.ones.zeros <- prop.correct.ones.zeros.acc.FP.FN.mat[c("TPR", "TNR"), , drop=FALSE]
  #no.correct.ones.zeros <- rbind(no.correct.ones, no.correct.zeros) # = correctly classified pairs = TP+TN
  if(M1==1) # i.e. a single estimated partition so no.correct.ones or zeros is a vector with 1 element
  { # proportions.correct.classified.pairs == prop.correct.ones.zeros.acc.FP.FN.mat["Acc", ]
    prop.correct.ones.zeros.acc.FP.FN.mat["Acc", ] <-
      (no.correct.ones.zeros["TP", ] + no.correct.ones.zeros["TN", ])/
      choose(ngen6, 2)
    no.correct.ones.zeros.acc.FP.FN.mat["Acc", ] <- no.correct.ones.zeros["TP", ] +
      no.correct.ones.zeros["TN", ]
  } else
  {
    prop.correct.ones.zeros.acc.FP.FN.mat["Acc", ] <- colSums(no.correct.ones.zeros)/choose(ngen6, 2)
    no.correct.ones.zeros.acc.FP.FN.mat["Acc", ] <- colSums(no.correct.ones.zeros)
  } # added 29-11-19 #
  prop.correct.ones.zeros.acc.FP.FN.mat.ari <- rbind(ari.vec, prop.correct.ones.zeros.acc.FP.FN.mat)
  rownames(prop.correct.ones.zeros.acc.FP.FN.mat.ari)[1] <- c("ARI") # NA means nonconvergent partition (e.g. with L2cont or whatever)
  #res$proportions.correct.classif.pairs <-  prop.correct.ones.zeros.acc.FP.FN.mat["Acc", ]
  no.nonconvergent.partitions.temp <- sum(is.na(achieved.convergence.run.vec1)) # added, 10-1-20 #
  if(sum(is.na(prop.correct.ones.zeros.acc.FP.FN.mat["TPR", ]))!= no.nonconvergent.partitions.temp) # just a check
    stop("The no. of NA in prop...['TPR',] is not equal to the sum of NAs in achieved.convergence...\n")
  res$ari.vec <- ari.vec # added 29-11-19. R clues package
  res$prop.correct.ones.zeros.acc.FP.FN.mat.ari <- prop.correct.ones.zeros.acc.FP.FN.mat.ari # 29-11-19 #res$proportions.Acc.TP.TN.FP.FN.mat <- prop.correct.ones.zeros.acc.FP.FN.mat # included in object above
  res$B.consensus.mat <- B.consensus.mat # of no.nonconverg. != 0 partitions
  res$M.effective <- M1 - no.nonconvergent.partitions.temp # added, 10-1-20 #
  res$no.nonconvergent.partitions <- no.nonconvergent.partitions.temp # added, 10-1-20 #
  res$no.Acc.TP.TN.FP.FN.mat <- no.correct.ones.zeros.acc.FP.FN.mat
  #res$no.true.positives.M1 <- no.correct.ones #res$no.true.negatives.M1 <- no.correct.zeros
  res$PercAccuracy <- prop.correct.ones.zeros.acc.FP.FN.mat["Acc", ]*100
  res$timing.vec <- timing.vec
  return(res)
}































