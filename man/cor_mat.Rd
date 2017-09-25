\name{cor_mat}
\alias{cor_mat}

\title{
Measuring PCC
}
\description{
Measuring PCC in the matrix
}
\usage{
cor_mat(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
matrix to measure correlation
}
}

\author{
Taegyun Yun
}

\examples{


## The function is currently defined as
function (mat) 
{
    adjustzero <- intersect(which(rowMeans(mat) == 0), which(rowVars(mat) == 
        0))
    if (length(adjustzero) != 0) {
        mat[adjustzero, ] <- c(rep(0.001, dim(mat)[2]))
    }
    adjustwhere <- which(rowVars(mat) == 0)
    if (length(adjustwhere) != 0) {
        for (jj in 1:length(adjustwhere)) {
            mat[adjustwhere[jj], ][1] <- mat[adjustwhere[jj], 
                ][1] + abs(mat[adjustwhere[jj], ][1]/10000)
        }
    }
    cormat <- cor(t(mat))
    homoRow <- mean(cormat[lower.tri(cormat)])
    return(homoRow = homoRow)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
