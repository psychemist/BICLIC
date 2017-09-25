\name{AddCondition}
\alias{AddCondition}

\title{
Column-wise expansion
}
\description{
Column-wise expansion in seed expansion step
}
\usage{
AddCondition(mgenes, mconds, corstand, input)
}

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
    new_mconds <- c()
    seq_vec_column <- setdiff(c(1:dim(input)[2]), mconds)
    seq_cor_column <- as.list(seq_vec_column)
    cor_result_column <- lapply(seq_cor_column, cor_mat_col_add, 
        mgenes, mconds)
    if (length(cor_result_column) != 0) {
        before_order_cor_column <- unlist(cor_result_column[which(cor_result_column >= 
            max(c(corstand/2, 0.5)))])
        if (length(before_order_cor_column) != 0) {
            before_order_column <- seq_vec_column[which(cor_result_column >= 
                max(c(corstand/2, 0.5)))]
            after_order_cor_column <- before_order_cor_column[order(before_order_cor_column, 
                decreasing = TRUE)]
            after_order_column <- before_order_column[order(before_order_cor_column, 
                decreasing = TRUE)]
            k <- 0
            loop_chk <- 0
            tmp_cor_vec <- c()
            candconds_tmp <- mconds
            while (loop_chk != -1) {
                k <- k + 1
                candconds_tmp <- c(candconds_tmp, after_order_column[k])
                candconds <- candconds_tmp[order(candconds_tmp, 
                  decreasing = F)]
                m2 <- input[mgenes, candconds]
                tmp_cor <- cor_mat(m2)
                tmp_cor_vec[k] <- tmp_cor
                if ((k == length(after_order_column)) || (tmp_cor_vec[k] < 
                  corstand)) {
                  loop_chk <- -1
                }
            }
            if (k == length(after_order_column)) {
                column_thres_point <- k
                cor_tmp <- tmp_cor_vec[column_thres_point]
            }
            if (tmp_cor_vec[k] < corstand) {
                column_thres_point <- k - 1
                cor_tmp <- tmp_cor_vec[column_thres_point]
            }
            if (column_thres_point != 0) {
                new_mconds <- after_order_column[1:column_thres_point]
            }
        }
    }
    return(new_mconds = new_mconds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

