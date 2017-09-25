\name{CandRowCol_plus}
\alias{CandRowCol_plus}

\title{
Seed finding
}
\description{
Seed finding
}
\usage{
CandRowCol_plus(i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
index
}
}

\author{
Taegyun Yun
}

\examples{

## The function is currently defined as
function (i) 
{
    sizechk <- 2
    sizechk2 <- 2
    tmprow <- which(clustermat_plus[, simatrow1[i]] == simatcol1[i])
    if (length(tmprow) >= sizechk) {
        tmpcol <- which(rowVars(t(clustermat_plus[tmprow, ])) == 
            0)
        if (length(tmpcol) >= sizechk2) {
            col_cnt <- length(tmpcol)
        }
        if (length(tmpcol) < sizechk2) {
            col_cnt <- 1
        }
    }
    if (length(tmprow) < sizechk) {
        tmpcol <- 1
        col_cnt <- 1
    }
    candrow <- tmprow[order(tmprow, decreasing = F)]
    candcol <- tmpcol[order(tmpcol, decreasing = F)]
    return(list(candrow = candrow, candcol = candcol, col_cnt = col_cnt))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
