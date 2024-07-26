\name{binary.matrix.truth.function}
\alias{binary.matrix.truth.function}
% binary_matrix_truth_function.Rd
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Binary Matrix: Clustered Pairs
%%  ~~function to do ... ~~
}
\description{
Given a partition true.tildetildeP of ngen elements and K clusters, find a binary matrix A with observations in rows and columns ordered by cluster such that those in the same cluster in true.tildetildeP appear together; 1 row and column genes are clustered together, 0 otherwise.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
binary.matrix.truth.function(K1=K, ngen1=ngen, true.tildetildeP1=true.tildetildeP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{K}{ The number of clusters.}
 \item{ngen}{ The number of observations.}
  \item{true.tildetildeP}{ A vector of length = ngen1 with the cluster membership values for each element. These should be integers from 1:K1.
%%     ~~Describe \code{x} here~~
}

}
\details{ If true.tildetildeP1 does not have names, the function generates them from 1:ngen1 and these are used in the matrix output as row and column names with the new order determined by imposing that all observations of the same cluster appear together.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A binary square symmetric matrix with dimensions ngen x ngen where 1 means that the row and column elements clustered together and 0 otherwise.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ Susana Conde
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
binary.matrix.truth.function(K1=3, ngen1=15, true.tildetildeP1=part1)

#    2 3 4 6 10 11 12 5 7 15 1 8 9 13 14
# 2  1 1 1 1  1  1  1 0 0  0 0 0 0  0  0
# 3  1 1 1 1  1  1  1 0 0  0 0 0 0  0  0
# 4  1 1 1 1  1  1  1 0 0  0 0 0 0  0  0
# 6  1 1 1 1  1  1  1 0 0  0 0 0 0  0  0
# 10 1 1 1 1  1  1  1 0 0  0 0 0 0  0  0
# 11 1 1 1 1  1  1  1 0 0  0 0 0 0  0  0
# 12 1 1 1 1  1  1  1 0 0  0 0 0 0  0  0
# 5  0 0 0 0  0  0  0 1 1  1 0 0 0  0  0
# 7  0 0 0 0  0  0  0 1 1  1 0 0 0  0  0
# 15 0 0 0 0  0  0  0 1 1  1 0 0 0  0  0
# 1  0 0 0 0  0  0  0 0 0  0 1 1 1  1  1
# 8  0 0 0 0  0  0  0 0 0  0 1 1 1  1  1
# 9  0 0 0 0  0  0  0 0 0  0 1 1 1  1  1
# 13 0 0 0 0  0  0  0 0 0  0 1 1 1  1  1
# 14 0 0 0 0  0  0  0 0 0  0 1 1 1  1  1


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
