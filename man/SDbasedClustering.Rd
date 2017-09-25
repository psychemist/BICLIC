\name{SDbasedClustering}
\alias{SDbasedClustering}

\title{
Individual dimension-based clustering
}
\description{
Individual dimension-based clustering to index genes
}
\usage{
SDbasedClustering(i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
index of each individual columns
}
}

\author{
Taegyun Yun
}

\examples{


## The function is currently defined as
function (i) 
{
    temp.col <- input[, i]
    sd_max <- sd(temp.col)
    if (length(table(temp.col)) == 1) {
        clvec <- c(rep(1, length(temp.col)))
        clvec2 <- c(rep(1, length(temp.col)))
    }
    if (length(table(temp.col)) != 1) {
        cl <- c(1)
        cnt <- 1
        indexed <- cbind(temp.col, c(1:length(temp.col)))
        indexedo <- indexed[order(indexed[, 1]), ]
        sentinel <- 0
        j <- 1
        while (sentinel != -1) {
            k <- 1
            sentinel2 <- 0
            if (j >= length(indexedo[, 1])) {
                sentinel <- -1
                sentinel2 <- -1
            }
            if ((j + k) == length(indexedo[, 1])) {
                cl <- c(cl, rep(cnt, (k + 1)))
                sentinel2 <- -1
                sentinel <- -1
                pass <- 1
            }
            if ((j + k) < length(indexedo[, 1])) {
                cl_target1 <- indexedo[j:(j + k), 1]
                sd_cl <- sd(cl_target1)
                if (sd_cl > sd_max) {
                  pass2 <- 0
                  if (j == 1) {
                    cnt <- cnt + 1
                    j <- j + 1
                    pass2 <- 1
                    sentinel2 <- -1
                  }
                  if ((j != 1) && (pass2 == 0)) {
                    cl <- c(cl, cnt)
                    cnt <- cnt + 1
                    j <- j + 1
                    pass2 <- 1
                    sentinel2 <- -1
                  }
                }
            }
            while (sentinel2 != -1) {
                pass <- 0
                k <- k + 1
                if ((j + k) == length(indexedo[, 1])) {
                  cl <- c(cl, rep(cnt, (k + 1)))
                  sentinel2 <- -1
                  sentinel <- -1
                  pass <- 1
                }
                cl_target2 <- indexedo[j:(j + k), 1]
                sd_cl2 <- sd(cl_target2)
                if (sd_cl2 <= sd_cl) {
                  sd_cl <- sd_cl2
                  pass <- 1
                }
                if ((sd_cl2 > sd_cl) && (pass == 0)) {
                  pass2 <- 0
                  if (j == 1) {
                    cl <- c(cl, rep(cnt, (k - 1)))
                    cnt <- cnt + 1
                    j <- j + k
                    pass2 <- 1
                  }
                  if ((j > 1) && (pass2 == 0)) {
                    cl <- c(cl, rep(cnt, k))
                    cnt <- cnt + 1
                    j <- j + k
                  }
                  sentinel2 <- -1
                }
            }
        }
        indexedo_tmp <- cbind(indexedo, cl)
        indexedo_tmpo <- indexedo_tmp[order(indexedo_tmp[, 2]), 
            ]
        clvec <- indexedo_tmpo[, 3]
    }
    return(list(clvec = clvec))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
