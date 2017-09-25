\name{DrawPattern}
\alias{DrawPattern}

\title{
Drawing patterns
}
\description{
Drawing correlated patterns of output biclusters
}
\usage{
DrawPattern(i, filename3, corstand)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
index of biclusters
}
  \item{filename3}{
the name of file to print
}
  \item{corstand}{
correlation threshold
}
}

\author{
Taegyun Yun
}

\examples{


## The function is currently defined as
function (i, filename3, corstand) 
{
    bi_filename3 <- paste("./result/", filename3, sep = "")
    bi_filename3 <- paste(bi_filename3, "_", sep = "")
    bi_filename3 <- paste(bi_filename3, corstand, sep = "")
    bi_filename3 <- paste(bi_filename3, "_bicluster_pattern#", 
        sep = "")
    bi_filename3 <- paste(bi_filename3, i, sep = "")
    bi_filename3 <- paste(bi_filename3, ".png", sep = "")
    patternname <- paste("bicluster_pattern #", i, sep = "")
    mat <- input[merged_rows_true[[i]], merged_cols_true[[i]]]
    x1 <- 1
    x2 <- dim(mat)[2]
    y1 <- min(mat)
    y2 <- max(mat)
    png(file = bi_filename3)
    matplot(c(x1, x2), c(y1, y2), type = "n", xlab = "Conditions", 
        ylab = "Expression level", main = patternname)
    for (j in 1:dim(mat)[1]) {
        matlines(as.numeric(mat[j, ]), col = 1)
    }
    dev.off()
    dev.next()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
