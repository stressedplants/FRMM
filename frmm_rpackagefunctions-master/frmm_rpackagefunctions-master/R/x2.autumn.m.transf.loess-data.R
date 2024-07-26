#' Arabidopsis halleri autumn functional gene expressions
#'
#' Data from Nature Plants paper, from an experiment done in
#' Arabidopsis halleri plants. Observations from autumn, during 48 hours
#' every other hour, starting from 4pm on. In total then, 24 time points.
#' It is a matrix with dimensions 5378 x 24 with the time points in columns
#' and genes in rows.
#'
#' @docType data
#'
#' @usage data(x2.autumn.m.transf.loess)
#'
#' @format An object of class \code{"matrix"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references Nagano et al. (2019) Nature Plants 5:74–83
#' (\href{https://www.nature.com/articles/s41477-018-0338-z}{PubMed})
#'
#' @source \href{https://www.nature.com/articles/s41477-018-0338-z}{QTL Archive}
#'
#' @examples
#' data(x2.autumn.m.transf.loess)
#' genelabels <- attr(x2.autumn.m.transf.loess, "dimnames")[[1]]
#' timepointlabels <- attr(x2.autumn.m.transf.loess, "dimnames")[[2]]
#' gene3spring <- x2.autumn.m.transf.loess[3, ]
"x2.autumn.m.transf.loess"
