\name{BICLIC}
\alias{BICLIC}

\title{
BICLIC function
}
\description{
The main function of BICLIC
}
\usage{
BICLIC(data, corstand, minrow, mincol, overlap_level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
input data matrix
}
  \item{corstand}{
correlation threshod
}
  \item{minrow}{
minimum number of rows
}
  \item{mincol}{
minimum number of columns
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
function (data, corstand, minrow, mincol, overlap_level) 
{
    input_tmp <- read.delim(data, sep = "\t")
    gname <- input_tmp[, 1]
    input <- input_tmp[, -1]
    input <<- as.matrix(input)
    dimInput <- dim(input)
    index <- c(1:dimInput[1])
    seq <- as.list(c(1:(dimInput[2])))
    sample_result1 <- lapply(seq, SDbasedClustering)
    sample_result2 <- lapply(seq, SDbasedClustering_minus)
    clustermat_plus <- c()
    numk_plus <- c()
    clustermat_minus <- c()
    numk_minus <- c()
    for (i in 1:(dimInput[2])) {
        clustermat_plus <- cbind(clustermat_plus, sample_result1[[i]]$clvec)
        numk_plus <- c(numk_plus, length(table(sample_result1[[i]]$clvec)))
        clustermat_minus <- cbind(clustermat_minus, sample_result2[[i]]$clvec)
        numk_minus <- c(numk_minus, length(table(sample_result2[[i]]$clvec)))
    }
    clustermat_plus <<- as.matrix(clustermat_plus)
    clustermat_minus <<- as.matrix(clustermat_minus)
    simatrow1 <- c()
    simatcol1 <- c()
    simatrow2 <- c()
    simatcol2 <- c()
    for (i in 1:dim(input)[2]) {
        simatrow1 <- c(simatrow1, rep(i, numk_plus[i]))
        simatcol1 <- c(simatcol1, c(1:numk_plus[i]))
        simatrow2 <- c(simatrow2, rep(i, numk_minus[i]))
        simatcol2 <- c(simatcol2, c(1:numk_minus[i]))
    }
    simatrow1 <<- simatrow1
    simatcol1 <<- simatcol1
    simatrow2 <<- simatrow2
    simatcol2 <<- simatcol2
    rowcol_result1 <- list()
    rowcol_result2 <- list()
    seq_simat1 <- as.list(c(1:length(simatrow1)))
    rowcol_result1 <- lapply(seq_simat1, CandRowCol_plus)
    rowcol_result2 <- lapply(seq_simat1, CandRowCol_minus)
    rowcol_result1 <- c(rowcol_result1, rowcol_result2)
    ptm <- proc.time()
    initial_rows <- list()
    initial_cols <- list()
    rownum <- c()
    colnum <- c()
    seed_cnt_tmp <- 0
    seed_tmp <- list()
    for (i in 1:length(rowcol_result1)) {
        if (rowcol_result1[[i]]$col_cnt != 1) {
            seed_cnt_tmp <- seed_cnt_tmp + 1
            seed_tmp[[seed_cnt_tmp]] <- c(rowcol_result1[[i]]$candrow, 
                0, rowcol_result1[[i]]$candcol)
        }
    }
    unique_seed <- unique(seed_tmp)
    if (length(unique_seed) != 0) {
        for (j in 1:length(unique_seed)) {
            vec1 <- unique_seed[[j]]
            threspoint <- which(vec1 == 0)
            initial_rows[[j]] <- vec1[1:(threspoint - 1)]
            initial_cols[[j]] <- vec1[(threspoint + 1):length(vec1)]
            rownum[j] <- length(initial_rows[[j]])
            colnum[j] <- length(initial_cols[[j]])
        }
    }
    initial_rows <<- initial_rows
    initial_cols <<- initial_cols
    seed_dim <- matrix(0, nrow = 3, ncol = length(unique_seed))
    seed_dim[1, ] <- rownum
    seed_dim[2, ] <- colnum
    seed_dim[3, ] <- rownum * colnum
    seed_h <- c()
    seed_cor <- c()
    seq_stat0 <- as.list(c(1:(length(initial_rows))))
    stat_result_cor0 <- lapply(seq_stat0, PCC_VEC_0)
    seed_cor <- unlist(stat_result_cor0)
    summary0 <- c(length(seed_dim[1, ]), mean(seed_dim[1, ]), 
        mean(seed_dim[2, ]), mean(seed_dim[3, ]), max(seed_dim[1, 
            ]), max(seed_dim[2, ]), max(seed_dim[3, ]), min(seed_dim[1, 
            ]), min(seed_dim[2, ]), min(seed_dim[3, ]), mean(seed_cor), 
        max(seed_cor), min(seed_cor))
    ptm2 <- proc.time() - ptm
    seed_covermat <- matrix(0, nrow = dim(input)[1], ncol = dim(input)[2])
    seed_num_tmp <- length(initial_rows)
    for (i in 1:length(initial_rows)) {
        seed_covermat[initial_rows[[i]], initial_cols[[i]]] <- seed_covermat[initial_rows[[i]], 
            initial_cols[[i]]] + 1
    }
    row_coverage_tmp <- rowMeans(seed_covermat)
    col_coverage_tmp <- colMeans(seed_covermat)
    row_coverage_num <- length(which(row_coverage_tmp != 0))
    col_coverage_num <- length(which(col_coverage_tmp != 0))
    cover_cnt <- 0
    for (i in 1:dim(seed_covermat)[2]) {
        cover_cnt <- cover_cnt + length(which(seed_covermat[, 
            i] == 0))
    }
    coverage_ratio <- 1 - cover_cnt/(dim(input)[1] * dim(input)[2])
    seed_cover_vec <- c(row_coverage_num/dim(input)[1], col_coverage_num/dim(input)[2], 
        coverage_ratio)
    seq_merge <- as.list(c(1:length(initial_rows)))
    ptm <- proc.time()
    merge_result <- lapply(seq_merge, MergeMat_dynamic2, corstand, 
        minrow, mincol, overlap_level)
    ptm3 <- proc.time() - ptm
    merge_num_tmp <- 0
    merged_rows_tmp <- list()
    merged_cols_tmp <- list()
    merged_rows_cols_tmp <- list()
    for (i in 1:length(merge_result)) {
        if (length(merge_result[[i]]$cor_vec) != 0) {
            for (j in 1:length(merge_result[[i]]$cor_vec)) {
                merge_num_tmp <- merge_num_tmp + 1
                merged_rows_tmp[[merge_num_tmp]] <- merge_result[[i]]$genes[[j]]
                merged_cols_tmp[[merge_num_tmp]] <- merge_result[[i]]$conds[[j]]
                merged_rows_cols_tmp[[merge_num_tmp]] <- c(merge_result[[i]]$genes[[j]], 
                  0, merge_result[[i]]$conds[[j]])
            }
        }
    }
    unique_merged_rows_cols <- unique(merged_rows_cols_tmp)
    merged_rows <<- list()
    merged_cols <<- list()
    bilist <- list()
    gleng <- c()
    cleng <- c()
    if (length(unique_merged_rows_cols) != 0) {
        for (j in 1:length(unique_merged_rows_cols)) {
            vec1 <- unique_merged_rows_cols[[j]]
            threspoint <- which(vec1 == 0)[1]
            merged_rows[[j]] <- vec1[1:(threspoint - 1)]
            merged_cols[[j]] <- vec1[(threspoint + 1):length(vec1)]
            bilist[[j]] <- input[merged_rows[[j]], merged_cols[[j]]]
            gleng[j] <- length(merged_rows[[j]])
            cleng[j] <- length(merged_cols[[j]])
        }
    }
    merged_rows <<- merged_rows
    merged_cols <<- merged_cols
    dimbi <- matrix(0, nrow = 3, ncol = length(gleng))
    dimbi[1, ] <- gleng
    dimbi[2, ] <- cleng
    dimbi[3, ] <- gleng * cleng
    num_tmp <- cbind(dimbi[3, ], dimbi[2, ], dimbi[1, ], c(1:length(dimbi[3, 
        ])))
    if (length(dimbi[3, ]) != 1) {
        num_order <- num_tmp[order(num_tmp[, 1], num_tmp[, 2], 
            num_tmp[, 3], decreasing = TRUE), ][, 4]
    }
    if (length(dimbi[3, ]) == 1) {
        num_order <- c(1)
    }
    h <- c()
    cor <- c()
    seq_stat1 <- as.list(c(1:(length(merged_rows))))
    stat_result_cor1 <- lapply(seq_stat1, PCC_VEC_1)
    cor <- unlist(stat_result_cor1)
    if (length(gleng) < 2) {
        summary1 <- c()
    }
    else {
        summary1 <- c(length(dimbi[1, ]), mean(dimbi[1, ]), mean(dimbi[2, 
            ]), mean(dimbi[3, ]), max(dimbi[1, ]), max(dimbi[2, 
            ]), max(dimbi[3, ]), min(dimbi[1, ]), min(dimbi[2, 
            ]), min(dimbi[3, ]), mean(cor), max(cor), min(cor))
    }
    covermat <- matrix(0, nrow = dim(input)[1], ncol = dim(input)[2])
    for (i in 1:length(merged_rows)) {
        covermat[merged_rows[[i]], merged_cols[[i]]] <- covermat[merged_rows[[i]], 
            merged_cols[[i]]] + 1
    }
    row_coverage_tmp <- rowMeans(covermat)
    col_coverage_tmp <- colMeans(covermat)
    row_coverage_num <- length(which(row_coverage_tmp != 0))
    col_coverage_num <- length(which(col_coverage_tmp != 0))
    cover_cnt <- 0
    for (i in 1:dim(covermat)[2]) {
        cover_cnt <- cover_cnt + length(which(covermat[, i] == 
            0))
    }
    coverage_ratio <- 1 - cover_cnt/(dim(input)[1] * dim(input)[2])
    tmp_cover_vec <- c(row_coverage_num/dim(input)[1], col_coverage_num/dim(input)[2], 
        coverage_ratio)
    order_vec <- num_order
    seq_duplicate <- as.list(c(1:(length(merged_rows) - 1)))
    order_vec <- order(dimbi[3, ], decreasing = TRUE)
    remove_duplicate_result <- lapply(seq_duplicate, CHECK_DUPLICATE, 
        order_vec, overlap_level)
    temp_chk <- c(rep(1, length(merged_rows)))
    for (i in 1:(length(merged_rows) - 1)) {
        temp_chk <- temp_chk * remove_duplicate_result[[i]]$chk_dup_vec
    }
    true_index <- order_vec[which(temp_chk == 1)]
    merged_rows_true <- list()
    merged_cols_true <- list()
    bilist_true <- list()
    cor_true <- c()
    gleng2 <- c()
    cleng2 <- c()
    for (i in 1:length(true_index)) {
        merged_rows_true[[i]] <- merged_rows[[true_index[i]]]
        merged_cols_true[[i]] <- merged_cols[[true_index[i]]]
        bilist_true[[i]] <- bilist[[true_index[i]]]
    }
    cor_true <- cor[true_index]
    gleng2 <- gleng[true_index]
    cleng2 <- cleng[true_index]
    dimbi_true <- dimbi[, true_index]
    if (length(gleng2) < 2) {
        summary2 <- c()
        summary_total <- c()
    }
    else {
        summary2 <- c(length(dimbi_true[1, ]), mean(dimbi_true[1, 
            ]), mean(dimbi_true[2, ]), mean(dimbi_true[3, ]), 
            max(dimbi_true[1, ]), max(dimbi_true[2, ]), max(dimbi_true[3, 
                ]), min(dimbi_true[1, ]), min(dimbi_true[2, ]), 
            min(dimbi_true[3, ]), mean(cor_true), max(cor_true), 
            min(cor_true))
        summary_total <- rbind(summary0, summary1, summary2)
    }
    covermat <- matrix(0, nrow = dim(input)[1], ncol = dim(input)[2])
    for (i in 1:length(merged_rows_true)) {
        covermat[merged_rows_true[[i]], merged_cols_true[[i]]] <- covermat[merged_rows_true[[i]], 
            merged_cols_true[[i]]] + 1
    }
    row_coverage_tmp <- rowMeans(covermat)
    col_coverage_tmp <- colMeans(covermat)
    row_coverage_num <- length(which(row_coverage_tmp != 0))
    col_coverage_num <- length(which(col_coverage_tmp != 0))
    cover_cnt <- 0
    for (i in 1:dim(covermat)[2]) {
        cover_cnt <- cover_cnt + length(which(covermat[, i] == 
            0))
    }
    coverage_ratio <- 1 - cover_cnt/(dim(input)[1] * dim(input)[2])
    cover_vec <- c(row_coverage_num/dim(input)[1], col_coverage_num/dim(input)[2], 
        coverage_ratio)
    which_max0 <- which(dimbi[3, ] == max(dimbi[3, ]))[1]
    which_min0 <- which(dimbi[3, ] == min(dimbi[3, ]))[1]
    bi_filename2 <- "summary_report.txt"
    cat("BICLIC_report", "\n\n", file = bi_filename2)
    cat("# of found biclusters: ", summary1[1], "\n", file = bi_filename2, 
        append = TRUE)
    cat("Maximum size is ", dimbi[1, which_max0], " * ", dimbi[2, 
        which_max0], " = ", dimbi[3, which_max0], "\n", file = bi_filename2, 
        append = TRUE)
    cat("Minimum size is ", dimbi[1, which_min0], " * ", dimbi[2, 
        which_min0], " = ", dimbi[3, which_min0], "\n", file = bi_filename2, 
        append = TRUE)
    cat("Average size is ", mean(dimbi[1, ]), " * ", mean(dimbi[2, 
        ]), " = ", mean(dimbi[3, ]), "\n", file = bi_filename2, 
        append = TRUE)
    cat("Average PCC: ", mean(cor), "\n\n", file = bi_filename2, 
        append = TRUE)
    which_max <- which(dimbi_true[3, ] == max(dimbi_true[3, ]))[1]
    which_min <- which(dimbi_true[3, ] == min(dimbi_true[3, ]))[1]
    cat("After overlapping removal", "\n\n", file = bi_filename2, 
        append = TRUE)
    cat("# of found biclusters: ", summary2[1], "\n", file = bi_filename2, 
        append = TRUE)
    cat("Maximum size is ", dimbi_true[1, which_max], " * ", 
        dimbi_true[2, which_max], " = ", dimbi_true[3, which_max], 
        "\n", file = bi_filename2, append = TRUE)
    cat("Minimum size is ", dimbi_true[1, which_min], " * ", 
        dimbi_true[2, which_min], " = ", dimbi_true[3, which_min], 
        "\n", file = bi_filename2, append = TRUE)
    cat("Average size is ", mean(dimbi_true[1, ]), " * ", mean(dimbi_true[2, 
        ]), " = ", mean(dimbi_true[3, ]), "\n", file = bi_filename2, 
        append = TRUE)
    cat("Average PCC: ", mean(cor_true), "\n", file = bi_filename2, 
        append = TRUE)
    cat("# of found seeds: ", summary0[1], "\n", file = bi_filename2, 
        append = TRUE)
    cat("Maximum size is ", seed_dim[1, which(seed_dim[3, ] == 
        max(seed_dim[3, ]))[1]], " * ", seed_dim[2, which(seed_dim[3, 
        ] == max(seed_dim[3, ]))[1]], " = ", seed_dim[3, which(seed_dim[3, 
        ] == max(seed_dim[3, ]))[1]], "\n", file = bi_filename2, 
        append = TRUE)
    cat("Minimum size is ", seed_dim[1, which(seed_dim[3, ] == 
        min(seed_dim[3, ]))[1]], " * ", seed_dim[2, which(seed_dim[3, 
        ] == min(seed_dim[3, ]))[1]], " = ", seed_dim[3, which(seed_dim[3, 
        ] == min(seed_dim[3, ]))[1]], "\n", file = bi_filename2, 
        append = TRUE)
    cat("Average size is ", mean(seed_dim[1, ]), " * ", mean(seed_dim[2, 
        ]), " = ", mean(seed_dim[3, ]), "\n", file = bi_filename2, 
        append = TRUE)
    cat("Average PCC: ", mean(seed_cor), "\n", file = bi_filename2, 
        append = TRUE)
    gname2 <- as.character(gname)
    for (i in 1:summary2[1]) {
        cat("Bicluster ", i, "\nsize = ", dimbi_true[3, i], ",  of gene = ", 
            dimbi_true[1, i], ",  of condition = ", dimbi_true[2, 
                i], ", PCC = ", cor_true[i], "\n", file = bi_filename2, 
            append = TRUE)
        cat("Genes are: ", gname2[merged_rows_true[[i]]], "\n", 
            file = bi_filename2, append = TRUE)
        cat("Conditions are: ", merged_cols_true[[i]], "\n\n", 
            file = bi_filename2, append = TRUE)
    }
    return(list(seed_cover_vec = seed_cover_vec, cover_vec = cover_vec, 
        tmp_cover_vec = tmp_cover_vec, gleng = gleng, cleng = cleng, 
        gleng2 = gleng2, cleng2 = cleng2, cor_true = cor_true, 
        summary_total = summary_total, bilist_true = bilist_true, 
        merged_rows_true = merged_rows_true, merged_cols_true = merged_cols_true, 
        initial_rows = initial_rows, initial_cols = initial_cols))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
