%% iterative_clustering_function.Rd
\name{iterative.clustering.function}
\alias{iterative.clustering.function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run of the Iterative Clustering Algorithm
%%  ~~function to do ... ~~
}
\description{From a random partition of K5 clusters, fit an FDboost model with Ytemp1 versus x1temp1, x2temp2, x3temp3 by cluster. Get the predicted values for all observations and all models. Update the partition that minimises the residual per observation. Repeat the procedure i.e. fit an FDboost model by cluster for the updated clusters, get the predicted values, etc. either until convergence (=find the same partition) or until reaching the maximum no. of iterations (m.it5).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
iterative.clustering.function(seed1=seed3, ngen5=ngen, genelab.m5=genelab.m, K5=K, m.gen.clu5=m.gen.clu, df.temp5=3,
 m.it5=m.it, Ytemp1=y.summer.m.transf.loess, x1temp1=x1.spring.m.transf.loess, x2temp1=x2.autumn.m.transf.loess,
 x3temp1=x3.winter.m.transf.loess, timep5=timep, l=1, nT=24, norm="l1", bspl.bas.no1=12, bspl.bas.ord1=4,
 n.knots=parent.frame()$n.knots, seed2=parent.frame()$seed2, ari.per.iteration=0, true.partit.for.ari=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed1}{number indicating the seed to draw a starting random partition P_0
%%     ~~Describe \code{x} here~~
}
\item{ngen5}{number of observations (= genes)}
\item{genelab.m5}{character vector of length=ngen5 with the observations (genes) labels}
\item{K5}{number of clusters}
\item{m.gen.clu5}{vector of K5 natural numbers that add up to ngen5 (= the sample sizes of the clusters of the starting random partition)}
\item{df.temp5}{the degrees of freedom used in the FDboost fit. It is essential that they are read indeed in the main environment prior to running the function}
\item{m.it5}{the (maximum) number of iterations}
\item{Ytemp1}{an ngen5 x nT matrix indicating a functional response variable}
\item{x1temp1, x2temp1, x3temp1}{ngen5 x nT matrices with the values of functional explanatory variables}
\item{timep5}{vector with length = nT indicating the time labels}
\item{l}{number of run in the whole FRMM algorithm procedure}
\item{nT}{number of time points}
\item{norm}{type of norm used in FRMM. It can be l1, l2 or l2cont}
\item{bspl.bas.no1, bspl.bas.ord1}{when norm="l2cont", the number and order of B-splines used in the smooth of the data respectively}
\item{n.knots}{number of knots of the FDboost calls}
\item{seed2}{seed set before the FDboost call}
\item{ari.per.iteration}{1 computes the \%ARI between the resulting partition in the l-th iteration and true.partit.for.ari, 0 otherwise; if 1, then it saves it in $ari.iteration.mat (4 x m.it5 matrix with the \%ARI, \%RI, \%TPR, \%TNR)}
\item{true.partit.for.ari}{if ari.per.iteration=0 it should be an NA. Otherwise, it is a vector of length= ngen5 with the cluster membership observed variable of the partition we want to compare to in order to get the ARI for the partition in each iteration}
}
\details{df.temp5 needs to be read in within the current environment prior to calling to this function.

    It needs the R FDboost, fda, fda.usc packages.

    The function uses the FDboost function and within it, one of the arguments is a call to boost_control, which is another function in FDboost. If trying to type in FDboost::boost_control(... etc.), R reports an error when trying to run the function: ``Error : 'boost_control' is not an exported object from 'namespace:FDboost' ''. This is overcome by removing FDboost:: in the definition of the function and calling to the FDboost library before executing it.

    Output: list with the final convergent partition (when existing) and the total sum of the norms of the residuals from fitting the FDboost models with P
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# FRMM run, n = 5378 gene expressions collected at spring, summer, autumn and winter during 24 times from Nagano et al (2019) #
# Set the sample sizes of an initial random partition in m.gen.clu and perform an FRMM run.
library(FDboost)
ntemp <- 300
m.gen.clu <- rep(100, 3)  #sum(m.gen.clu) = nrow(y.summer.m.transf.loess)
df.temp5 <- 3
system.time(
  a <- iterative.clustering.function(seed1=1, ngen5=ntemp, genelab.m5=rownames(y.summer.m.transf.loess)[1:ntemp], K5=3, m.gen.clu5=m.gen.clu, df.temp5=3,
 m.it5=300, Ytemp1=y.summer.m.transf.loess[1:ntemp, ], x1temp1=x1.spring.m.transf.loess[1:ntemp, ], x2temp1=x2.autumn.m.transf.loess[1:ntemp, ],
 x3temp1=x3.winter.m.transf.loess[1:ntemp, ], timep5=1:24, l=1, nT=24, norm="l2", bspl.bas.no1=12, bspl.bas.ord1=4,
 n.knots=6, seed2=5, ari.per.iteration=0, true.partit.for.ari=NA)
)
names(a)
a$final.clme.l[1:5, , drop=FALSE]
table(a$final.clme.l)

#           it300
# Ahg311263     2
# Ahg311286     2
# Ahg311317     3
# Ahg311324     3
# Ahg311338     3
# 1   2   3
# 115  72 113
# If using the 5378 genes:
# Ahg311263 Ahg311286 Ahg311317 Ahg311324 Ahg311338
#         2         2         3         2         2
#    1    2    3
# 1210 2296 1872
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
