\name{MergeMat_dynamic2}
\alias{MergeMat_dynamic2}

\title{
Expanding seed biclusters
}
\description{
Expanding seed biclusters to larger size of biclusters and filtering less correlated genes and columns
}
\usage{
MergeMat_dynamic2(i, corstand, minrow, mincol, overlap_level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
index of seed biclusters
}
  \item{corstand}{
correlation threshold
}
  \item{minrow}{
minimum number of rows of biclusters to expand
}
  \item{mincol}{
minimum number of column of biclusters to expand
}
  \item{overlap_level}{
overlap level
}
}

\author{
Taegyun Yun
}

\examples{

## The function is currently defined as
function (i, corstand, minrow, mincol, overlap_level) 
{
    originalG <- initial_rows[[i]]
    columnP <- initial_cols[[i]]
    mgenes <- originalG
    mconds <- columnP
    mat0 <- input[originalG, columnP]
    mat0 <- as.matrix(mat0)
    cor0 <- cor_mat(mat0)
    mrow <- c()
    mcol <- c()
    dimInput <- dim(input)
    overstand <- 0
    mgenes_list <- list()
    mconds_list <- list()
    total_num <- 0
    mgenes_mconds_list <- list()
    unique_genes <- list()
    unique_conds <- list()
    genes <- list()
    conds <- list()
    cor_vec <- c()
    h_vec <- c()
    genes_tmp <- list()
    conds_tmp <- list()
    cor_vec_tmp <- c()
    if (cor0 >= corstand) {
        if (length(columnP) < mincol) {
            seq_vec_column <- setdiff(c(1:dimInput[2]), mconds)
            seq_cor_column <- as.list(seq_vec_column)
            cor_result_column <- lapply(seq_cor_column, cor_mat_col_add, 
                mgenes, mconds)
            before_order_cor_column <- unlist(cor_result_column[which(cor_result_column >= 
                max(c(corstand/2, 0.5)))])
            before_order_column <- seq_vec_column[which(cor_result_column >= 
                max(c(corstand/2, 0.5)))]
            after_order_cor_column <- before_order_cor_column[order(before_order_cor_column, 
                decreasing = TRUE)]
            after_order_column <- before_order_column[order(before_order_cor_column, 
                decreasing = TRUE)]
            num_to_col <- mincol - length(columnP)
            if (length(after_order_column) >= mincol - length(columnP)) {
                mconds_tmp <- c(columnP, after_order_column[1:num_to_col])
                mconds <- mconds_tmp[order(mconds_tmp, decreasing = F)]
                mat0 <- input[originalG, mconds]
                cor0 <- cor_mat(mat0)
            }
            else {
                cor0 <- 0
            }
        }
        if (cor0 >= corstand) {
            if ((length(mgenes) >= minrow) && (length(mconds) >= 
                mincol)) {
                total_num <- total_num + 1
                mgenes_list[[total_num]] <- mgenes
                mconds_list[[total_num]] <- mconds
                mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 
                  0, mconds_list[[total_num]])
            }
            cond_max_result <- AddCondition(mgenes, mconds, corstand, 
                input)
            gene_tmp <- mgenes
            cond_tmp <- c(mconds, cond_max_result)
            if ((length(gene_tmp) >= minrow) && (length(cond_tmp) >= 
                mincol)) {
                total_num <- total_num + 1
                mgenes_list[[total_num]] <- gene_tmp[order(gene_tmp, 
                  decreasing = F)]
                mconds_list[[total_num]] <- cond_tmp[order(cond_tmp, 
                  decreasing = F)]
                mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 
                  0, mconds_list[[total_num]])
            }
            gene_max_result <- AddGene(mgenes, mconds, corstand, 
                input)
            gene_tmp <- c(mgenes, gene_max_result)
            cond_tmp <- mconds
            if ((length(gene_tmp) >= minrow) && (length(cond_tmp) >= 
                mincol)) {
                total_num <- total_num + 1
                mgenes_list[[total_num]] <- gene_tmp[order(gene_tmp, 
                  decreasing = F)]
                mconds_list[[total_num]] <- cond_tmp[order(cond_tmp, 
                  decreasing = F)]
                mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 
                  0, mconds_list[[total_num]])
            }
            gene_tmp <- c(mgenes, gene_max_result)
            cond_tmp <- c(mconds, cond_max_result)
            gene_tmp <- gene_tmp[order(gene_tmp, decreasing = F)]
            cond_tmp <- cond_tmp[order(cond_tmp, decreasing = F)]
            tmp_mat <- input[gene_tmp, cond_tmp]
            gene_tmp2 <- gene_tmp
            cond_tmp2 <- cond_tmp
            tmp_mat <- input[gene_tmp, cond_tmp]
            tmp_mat2 <- tmp_mat
            cor_initial <- cor_mat(tmp_mat)
            if (cor_initial >= 0.999) {
                cor_initial <- 1
            }
            before_cor <- cor_initial
            sentinel <- 0
            roop_cnt <- 0
            if (cor_initial < corstand) {
                while (sentinel != -1) {
                  roop_cnt <- roop_cnt + 1
                  if (length(gene_tmp2) > 2) {
                    gene_cor_vec <- numeric(length(gene_tmp2))
                    cond_cor_vec <- numeric(length(cond_tmp2))
                    for (j in 1:length(gene_tmp2)) {
                      gene_cor_vec[j] <- cor_mat(input[gene_tmp2[-j], 
                        cond_tmp2])
                    }
                    gene_diff_vec <- gene_cor_vec - before_cor
                    cand_gene_remove <- which(gene_diff_vec > 
                      0)
                    if ((length(gene_tmp2) - length(cand_gene_remove)) < 
                      2) {
                      cand_gene_remove2_tmp <- gene_diff_vec[cand_gene_remove]
                      cand_gene_remove2 <- cand_gene_remove[order(cand_gene_remove2_tmp, 
                        decreasing = T)]
                      gene_cor_vec2 <- numeric(length(cand_gene_remove2))
                      for (k in 1:(length(gene_tmp2) - 2)) {
                        gene_cor_vec2[k] <- cor_mat(input[gene_tmp2[-cand_gene_remove2[k]], 
                          cond_tmp2])
                      }
                      after_cor_gene <- gene_cor_vec2[1]
                      diff_cor_gene <- after_cor_gene - before_cor
                      cand_gene_remove <- cand_gene_remove2[1]
                    }
                    else {
                      after_cor_gene <- cor_mat(input[gene_tmp2[-cand_gene_remove], 
                        cond_tmp2])
                      diff_cor_gene <- after_cor_gene - before_cor
                    }
                  }
                  else {
                    diff_cor_gene <- -1
                    after_cor_gene <- 0
                  }
                  for (j in 1:length(cond_tmp)) {
                    cond_cor_vec[j] <- cor_mat(input[gene_tmp2, 
                      cond_tmp2[-j]])
                  }
                  cond_diff_vec <- cond_cor_vec - before_cor
                  after_cor_cond <- max(cond_cor_vec)
                  diff_cor_cond <- max(cond_diff_vec)
                  cand_cond_remove <- which(cond_diff_vec == 
                    max(cond_diff_vec))
                  if ((diff_cor_gene <= 0) && (diff_cor_cond <= 
                    0)) {
                    sentinel <- -1
                  }
                  else if (((diff_cor_gene > 0) && (diff_cor_cond <= 
                    0)) || (after_cor_gene > after_cor_cond)) {
                    gene_tmp3 <- gene_tmp2
                    gene_tmp2 <- gene_tmp2[-cand_gene_remove]
                    tmp_mat2 <- input[gene_tmp2, cond_tmp2]
                    before_cor <- cor_mat(tmp_mat2)
                    if (before_cor >= corstand) {
                      cand_gene_remove2_tmp <- gene_diff_vec[cand_gene_remove]
                      cand_gene_remove2 <- cand_gene_remove[order(cand_gene_remove2_tmp, 
                        decreasing = T)]
                      gene_cor_vec2 <- numeric(length(cand_gene_remove2))
                      for (k in 1:length(cand_gene_remove2)) {
                        gene_cor_vec2[k] <- cor_mat(input[gene_tmp3[-cand_gene_remove2[1:k]], 
                          cond_tmp2])
                      }
                      gene_tmp2 <- gene_tmp3[-cand_gene_remove2[1:which(gene_cor_vec2 >= 
                        corstand)[1]]]
                      before_cor <- cor_mat(input[gene_tmp2, 
                        cond_tmp2])
                      sentinel <- -1
                    }
                  }
                  else if (((diff_cor_gene <= 0) && (diff_cor_cond > 
                    0)) || (after_cor_gene < after_cor_cond)) {
                    cond_tmp3 <- cond_tmp2
                    cond_tmp2 <- cond_tmp3[-cand_cond_remove]
                    tmp_mat2 <- input[gene_tmp2, cond_tmp2]
                    before_cor <- cor_mat(tmp_mat2)
                    if (before_cor >= corstand) {
                      sentinel <- -1
                    }
                  }
                }
            }
            else {
                if ((length(gene_tmp) >= minrow) && (length(cond_tmp) >= 
                  mincol)) {
                  total_num <- total_num + 1
                  mgenes_list[[total_num]] <- gene_tmp[order(gene_tmp, 
                    decreasing = F)]
                  mconds_list[[total_num]] <- cond_tmp[order(cond_tmp, 
                    decreasing = F)]
                  mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 
                    0, mconds_list[[total_num]])
                }
            }
            if (((length(gene_tmp2) >= minrow) && (length(cond_tmp2) >= 
                mincol)) && (before_cor > corstand)) {
                total_num <- total_num + 1
                mgenes_list[[total_num]] <- gene_tmp2[order(gene_tmp2, 
                  decreasing = F)]
                mconds_list[[total_num]] <- cond_tmp2[order(cond_tmp2, 
                  decreasing = F)]
                mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 
                  0, mconds_list[[total_num]])
            }
            unique_mgenes_mconds <- unique(mgenes_mconds_list)
            if (length(unique_mgenes_mconds) != 0) {
                gleng <- c()
                cleng <- c()
                for (j in 1:length(unique_mgenes_mconds)) {
                  vec1 <- unique_mgenes_mconds[[j]]
                  threspoint <- which(vec1 == 0)[1]
                  unique_genes[[j]] <- vec1[1:(threspoint - 1)]
                  unique_conds[[j]] <- vec1[(threspoint + 1):length(vec1)]
                  gleng[j] <- length(unique_genes[[j]])
                  cleng[j] <- length(unique_conds[[j]])
                  cor_vec_tmp[j] <- cor_mat(input[unique_genes[[j]], 
                    unique_conds[[j]]])
                }
                dimbi <- matrix(0, nrow = 3, ncol = length(gleng))
                dimbi[1, ] <- gleng
                dimbi[2, ] <- cleng
                dimbi[3, ] <- gleng * cleng
                num_tmp <- cbind(dimbi[3, ], dimbi[2, ], dimbi[1, 
                  ], c(1:length(dimbi[3, ])))
                if (length(dimbi[3, ]) != 1) {
                  num_order <- num_tmp[order(num_tmp[, 1], num_tmp[, 
                    2], num_tmp[, 3], decreasing = TRUE), ][, 
                    4]
                }
                if (length(dimbi[3, ]) == 1) {
                  num_order <- c(1)
                }
                if (length(gleng) != 1) {
                  overlap_result <- CHECK_DUPLICATE2(num_order, 
                    unique_genes, unique_conds, overlap_level)
                  temp_chk <- overlap_result$chk_dup_vec
                  if (length(temp_chk) != 0) {
                    true_index <- num_order[which(temp_chk == 
                      1)]
                    for (j in 1:length(true_index)) {
                      genes[[j]] <- unique_genes[[true_index[j]]]
                      conds[[j]] <- unique_conds[[true_index[j]]]
                      cor_vec[j] <- cor_vec_tmp[true_index[j]]
                    }
                  }
                }
                if (length(gleng) == 1) {
                  genes <- unique_genes
                  conds <- unique_conds
                  cor_vec <- cor_vec_tmp
                }
            }
        }
    }
    return(list(genes = genes, conds = conds, cor_vec = cor_vec))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
