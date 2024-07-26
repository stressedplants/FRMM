\name{binary.vectors.lth.partition.function}  %% binary_vectors_lth_partition_function.Rd
\alias{binary.vectors.lth.partition.function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Matrix with (individual) binary cluster membership variables by rows
%%  ~~function to do ... ~~
}
\description{
Given a partition ^P(l)=hatP.l2 of ngen2 elements into K2 clusters, find a binary vector v_k[i] with 1 ith observation (=gene) belongs to the kth cluster, 0 otherwise i.e. a "dummy" variable for the kth cluster. Store the K2 vectors in a matrix which is K2 x ngen2 and return this matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
binary.vectors.lth.partition.function(K2=K, ngen2=ngen, hatP.l2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K2}{ The number of clusters.}
  \item{ngen2}{ The number of observations (e.g. genes).}
  \item{hatP.l2}{ A vector of length = ngen2 with the cluster membership values for each element. These should be integers from 1:K2. % It can have names or it cannot
%%     ~~Describe \code{x} here~~
}

}
\details{ hatP.l2 can either have names attribute or not.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A binary matrix with dimensions K2 x ngen2 where the kth row has the cluster membership variable for the kth cluster, k=1, ..., K2.
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
part1 <- c(3, 1, 1, 1, 2, 1, 2, 3, 3, 1, 1, 1, 3, 3, 2) % , 3, 3, 3, 2, 2)
binary.vectors.lth.partition.function(K2=3, ngen2=15, hatP.l2=part1)

#  [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15]
#1    0    1    1    1    0    1    0    0    0     1     1     1     0     0     0
#2    0    0    0    0    1    0    1    0    0     0     0     0     0     0     1
#3    1    0    0    0    0    0    0    1    1     0     0     0     1     1     0
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
