\name{CHECK_DUPLICATE2}
\alias{CHECK_DUPLICATE2}

\title{
Checking overlap level of biclusters in the merging function
}
\description{
Measures the overlap level of output biclusters in the mergeMatDynamic2 function
}
\usage{
CHECK_DUPLICATE2(order_vec, merged_rows_func, merged_cols_func, overlap_level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{order_vec}{
vector containing the decreasing order of bicluster size
}
  \item{merged_rows_func}{
expanded rows to check overlap
}
  \item{merged_cols_func}{
expanded columns to check overlap
}
  \item{overlap_level}{
overlap threshold
}
}

\author{
Taegyun Yun
}

\examples{

## The function is currently defined as
function (order_vec, merged_rows_func, merged_cols_func, overlap_level) 
{
    cnt <- 0
    chk_dup_vec <- c(rep(1, length(merged_rows_func)))
    for (i in 1:(length(merged_rows_func) - 1)) {
        if (chk_dup_vec[i] != 0) {
            row1 <- merged_rows_func[[order_vec[i]]]
            col1 <- merged_cols_func[[order_vec[i]]]
            area_vec <- c(rep(0, length(chk_dup_vec)))
            area1 <- length(row1) * length(col1)
            for (j in 1:length(merged_rows_func)) {
                if ((chk_dup_vec[j] != 0) && (i != j)) {
                  cnt <- cnt + 1
                  row2 <- merged_rows_func[[order_vec[j]]]
                  col2 <- merged_cols_func[[order_vec[j]]]
                  row_inter <- intersect(row1, row2)
                  col_inter <- intersect(col1, col2)
                  area_inter <- length(row_inter) * length(col_inter)
                  area2 <- length(row2) * length(col2)
                  if ((area1 > area2) || ((area1 == area2) && 
                    (i < j))) {
                    area_frac <- area_inter/area2
                    area_vec[j] <- area_frac
                    if (area_frac >= overlap_level) {
                      chk_dup_vec[j] <- 0
                    }
                  }
                }
            }
        }
    }
    return(list(chk_dup_vec = chk_dup_vec, area_vec = area_vec))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
