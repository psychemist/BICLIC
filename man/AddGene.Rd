\name{AddGene}
\alias{AddGene}

\title{
Gene-wise expansion
}
\description{
Gene-wise expansion in seed expansion step
}
\usage{
AddGene(mgenes, mconds, corstand, input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mgenes}{
candidate genes to expand
}
  \item{mconds}{
     candidate conditions to expand
}
  \item{corstand}{
     correlation threshold
}
  \item{input}{
     input data matrix
}
}

\author{
 Taegyun Yun
}

\examples{

## The function is currently defined as
function (mgenes, mconds, corstand, input) 
{
    new_mgenes <- c()
    seq_vec_row <- setdiff(c(1:dim(input)[1]), mgenes)
    seq_cor_row <- as.list(seq_vec_row)
    cor_result_row <- lapply(seq_cor_row, cor_mat_row_add, mgenes, 
        mconds)
    if (length(cor_result_row) != 0) {
        before_order_cor_row <- unlist(cor_result_row[which(cor_result_row >= 
            max(c(corstand/2, 0.5)))])
        if (length(before_order_cor_row) != 0) {
            before_order_row <- seq_vec_row[which(cor_result_row >= 
                max(c(corstand/2, 0.5)))]
            after_order_cor_row <- before_order_cor_row[order(before_order_cor_row, 
                decreasing = TRUE)]
            after_order_row <- before_order_row[order(before_order_cor_row, 
                decreasing = TRUE)]
            k <- 0
            loop_chk <- 0
            tmp_cor_vec <- c()
            candgenes_tmp <- mgenes
            while (loop_chk != -1) {
                k <- k + 1
                candgenes_tmp <- c(candgenes_tmp, after_order_row[k])
                candgenes <- candgenes_tmp[order(candgenes_tmp, 
                  decreasing = F)]
                m2 <- input[candgenes, mconds]
                tmp_cor <- cor_mat(m2)
                tmp_cor_vec[k] <- tmp_cor
                if ((k == length(after_order_row)) || (tmp_cor_vec[k] < 
                  corstand)) {
                  loop_chk <- -1
                }
            }
            if (k == length(after_order_row)) {
                row_thres_point <- k
                cor_tmp <- tmp_cor_vec[row_thres_point]
            }
            if (tmp_cor_vec[k] < corstand) {
                row_thres_point <- k - 1
                cor_tmp <- tmp_cor_vec[row_thres_point]
            }
            if (row_thres_point != 0) {
                new_mgenes <- after_order_row[1:row_thres_point]
            }
        }
    }
    return(new_mgenes = new_mgenes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
