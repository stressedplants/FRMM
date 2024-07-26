\name{prop.correct.pairs.lth.partition.and.truth.function}
\alias{prop.correct.pairs.lth.partition.and.truth.function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adjusted Rand Index and Consensus Matrix between Two Partitions
%%  ~~function to do ... ~~
}
\description{ Given two partitions, get (essentially) the adjusted Rand index as well as the consensus matrix, plus the Rand index, true positive rate and true negative rate. Given an ``index'' partition and another set of partitions, get the adjusted Rand index, TPR, TNR using the consensus matrix calculated from the partitions in the columns in final.clme.per.run.mat1.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
prop.correct.pairs.lth.partition.and.truth.function(M1=1, K6=K, ngen6=ngen,
 final.clme.per.run.mat1=final.clme.per.run.mat, genelab.m6=NULL, true.tildetildeP6=true.tildetildeP,
 m.gen.clu6=as.vector(table(true.tildetildeP)), genelab.ordered.true.partition6=NULL, # genelab.ordered.true.partition
 achieved.convergence.run.vec1=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{M1}{The number of columns of final.clme.per.run.mat1 = number of partitions we will compare with true.tildetildeP6}
 \item{K6}{The number of clusters of true.tildetildeP6. }
 \item{ngen6}{The number of observations. It has to coincide with the number of rows of final.clme.per.run.mat1, or the length(true.tildetildeP6), of the length(genelab.m6) or length(genelab.ordered.true.partition6) if they are nonnull}
 \item{final.clme.per.run.mat1}{A matrix of dimensions ngen6 x M1. If the rownames are NULL then in theory the names of true.tildetildeP6 or genelab.m6 should also be NULL (and we assume that the row order in the matrix coincides with the order of the elements in true.tildetildeP6). Its entries should be integers from 1 to K6, indicating the cluster membership by observation for each column partition.}
 \item{genelab.m6}{If nonnull, a character vector that has to coincide with names(true.tildetildeP6). If null, then we assume that the order of the rows in final.clme... is the same as the order of the elements in true.tildetildeP6.}
 \item{true.tildetildeP6}{An integer vector indicating the cluster membership of the partition against which we want to calculate the adjusted Rand Index, Rand Index, True Positive Rate and True Negative Rate when compared with each column in final.clme.per.run.mat1. Its entries have to be numbers from 1 to K6.}
\item{achieved.convergence.run.vec1}{A vector of length = ngen6. If NA, then the corresponding column in final.clme.per.run.mat1 is kept apart and it does not participate in the calculation of the consensus matrix. }
%%     ~~Describe \code{x} here~~
}
\details{ It needs the R clues package (as it uses adjustedRand).
The dimensions of final.clme.per.run.mat1 have to be ngen6 x M1.
If genelab.ordered.true.partition6 = NULL then the names(true.tildetildeP6) should be NULL too (and probably rownames(final.clme...)) but even if those are not NULL, the function might work correctly too (IF THE ORDER OF THE ELEMENTS IN TRUE.TILDE... IS THE SAME AS THAT IN ROWNAMES(FINAL.CLME...)).
If genelab.m6 is not null then it should be exactly names(true.tildetildeP6) and in the same order, = rownames(final.clme...).
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with 10 components:

prop.correct.ones.zeros.acc.FP.FN.mat.ari:  6xM1 matrix with the ARI*100, percentage of RI, \%TPR, \%TNR, \%FPR, \%FNR by rows,

B.consensus.mat[i, j]: consensus matrix calculated with the partitions in the columns of final.clme.per.run.mat1 (it only includes results for the convergent partitions, which are the nonNA values in achieved.convergence.run.vec1),

M.effective=M1 - no. of nonconverg. partitions (i.e. latter possibly !=0),

no.nonconvergent.partitions (vector length=M1), plus:

no.Acc.TP.TN.FP.FN.mat,

ari.vec[l],

true.tildetildeP,

hatP.l,

PercAccuracy,

timing.vec.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Susana Conde
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{iterative.clustering.function}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# 2 partitions of 15 observations with 3 clusters each. Get the consensus matrix plus the ARI, RI, TPR, TNR:
part1 <- c(3, 1, 1, 1, 2, 1, 2, 3, 3, 1, 1, 1, 3, 3, 2)
part2 <-  c(2, 2, 2, 1, 3, 1, 2, 3, 2, 3, 1, 1, 1, 3, 2)
a <- prop.correct.pairs.lth.partition.and.truth.function(M1=1, K6=3, ngen6=15,
 final.clme.per.run.mat1=matrix(part1), genelab.m6=NULL, true.tildetildeP6=part2,
 m.gen.clu6=as.vector(table(part2)), genelab.ordered.true.partition6=NULL,
 achieved.convergence.run.vec1=1)
a$prop.correct.ones.zeros.acc.FP.FN.mat.ari[1:4, ]
a$B.consensus.mat
#               l1
# ARI -0.001695993
# Acc  0.533333333
# TPR  0.258064516
# TNR  0.648648649
#
#    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  1 0 0 0 0 0 0 1 1  0  0  0  1  1  0
# 2  0 1 1 1 0 1 0 0 0  1  1  1  0  0  0
# 3  0 1 1 1 0 1 0 0 0  1  1  1  0  0  0
# 4  0 1 1 1 0 1 0 0 0  1  1  1  0  0  0
# 5  0 0 0 0 1 0 1 0 0  0  0  0  0  0  1
# 6  0 1 1 1 0 1 0 0 0  1  1  1  0  0  0
# 7  0 0 0 0 1 0 1 0 0  0  0  0  0  0  1
# 8  1 0 0 0 0 0 0 1 1  0  0  0  1  1  0
# 9  1 0 0 0 0 0 0 1 1  0  0  0  1  1  0
# 10 0 1 1 1 0 1 0 0 0  1  1  1  0  0  0
# 11 0 1 1 1 0 1 0 0 0  1  1  1  0  0  0
# 12 0 1 1 1 0 1 0 0 0  1  1  1  0  0  0
# 13 1 0 0 0 0 0 0 1 1  0  0  0  1  1  0
# 14 1 0 0 0 0 0 0 1 1  0  0  0  1  1  0
# 15 0 0 0 0 1 0 1 0 0  0  0  0  0  0  1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
