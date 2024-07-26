\name{consensus.matrix.fun}
\alias{consensus.matrix.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Consensus Matrix
%%  ~~function to do ... ~~
}
\description{ Given some partition/s (= no. of columns of final.clme.per.run.mat3), compute the consensus matrix, whose entries indicate the proportion of time that the row and column observations cluster together in the partition/s.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
consensus.matrix.fun(M.run3 = 1, final.clme.per.run.mat3 = final.clme.per.run.mat, ngen3 = ngen, genelab.m3 = genelab.m, K3 = K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M.run3}{number of partitions we want to compute the consensus matrix from.
}
  \item{final.clme.per.run.mat3}{matrix with dimensions ngen3 x M.run3 containing the cluster membership of the M.run3 partitions by columns.
}
  \item{ngen3}{number of observations
}
  \item{genelab.m3}{character vector with the observation labels in the order as they appear in the rows in final.clme.per.run.mat3
}
  \item{K3}{number of clusters
}
}
\details{Call to binary.vectors.lth.partition.function.
}
\value{It returns a list a list with 3 components:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{A.l.binary.hatP}{a sublist with each of the binary matrices indicating whether a pair of observations belong to the same cluster in each partition}
\item{B.consensus.counts}{the consensus matrix in the form of counts}
\item{B.consensus.mat}{the consensus matrix in the form of proportions}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Susana Conde}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{binary.vectors.lth.partition.function}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
part1 <- matrix(c(3, 1, 1, 1, 2, 1, 2, 3, 3, 1, 1, 1, 3, 3, 2), ncol=1)
part2 <-  c(2, 2, 2, 1, 3, 1, 2, 3, 2, 3, 1, 1, 1, 3, 2)
b <- cbind(part1, part2)
a <- consensus.matrix.fun(M.run3 = 1, final.clme.per.run.mat3 = part1, ngen3 = nrow(part1), genelab.m3 =NULL, K3 = length(table(part1)))
a$B.consensus.mat
a <- consensus.matrix.fun(M.run3 = 2, final.clme.per.run.mat3 = cbind(part1, part2), ngen3 = nrow(part1), genelab.m3 =NULL, K3 = length(table(part1)))
a$B.consensus.mat
rownames(b) <- paste0("Ahg", 1:15)
a <- consensus.matrix.fun(M.run3 = 2, final.clme.per.run.mat3 = b, ngen3 = nrow(part1), genelab.m3 =NULL, K3 = length(table(part1)))
a$B.consensus.mat
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
