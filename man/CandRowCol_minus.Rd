\name{CandRowCol_minus}
\alias{CandRowCol_minus}

\title{
Seed finding 2
}
\description{
Second version of seed finding
}
\usage{
CandRowCol_minus(i)
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
    tmprow <- which(clustermat_minus[, simatrow2[i]] == simatcol2[i])
    if (length(tmprow) >= sizechk) {
        tmpcol <- which(rowVars(t(clustermat_minus[tmprow, ])) == 
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
