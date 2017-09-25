#The function performs BICLIC algorithm by by using clustered seeds and their expansion with the correlation of gene expression
#arguments:
#
#	data 				    The input matrix to be biclustered. maxtrix should be tab-delimited text file.
# corstand        The correlation threshold to find correlated expression
# minrow          The minimum number of row for final bicluster matrix
# mincon          The minimum number of condition for final bicluster matrix. This should be larger than 2 
# overlap_level   The threshold to filter out overlapped biclusters, 1 means filtering 100% overlapped biclusters
#
#
# example run:
# biclic_result <- BICLIC("YeastGalactose.txt", 0.9, 5, 5, 1)
#
#
# output file: "summary_report.txt" file 
# 
#
# output objects:
#  seed_cover_vec      Coverage of seeds in gene, condition, and cell dimension
#  cover_vec           Coverage of biclusters in gene, condition, and cell dimension
#  tmp_cover_Vec       Coverage of biclusters before overlapped bicluster removal in gene, condition, and cell dimension      
#  gleng2               Vector contating the number of genes of output biclusters
#  cleng2 			   Vector containing the number of conditions of output biclusters 
#  cor_true        Vectors containing the PCC of output biclusters
#  summary_total   Matrix containing the summary statistics of BICLIC output
#  blilist         List containing the output biclusters
#  merged_rows_true   List containing the rows of output biclusters
#  merged_cols_true   List contaning the columns of output biclusters

library(matrixStats)


BICLIC <- function(data, corstand, minrow, mincol, overlap_level)
{	
	input_tmp <- read.delim(data, sep="\t")
	
  gname <- input_tmp[,1]
  input <- input_tmp[,-1] 
  input <<- as.matrix(input)
  
 	dimInput <- dim(input)
	index <- c(1:dimInput[1])
	seq<-as.list(c(1:(dimInput[2])))
             
  ## one dimensional clustering
  sample_result1<-lapply(seq, SDbasedClustering)
  sample_result2<-lapply(seq, SDbasedClustering_minus)
  
 	clustermat_plus<-c()
	numk_plus <- c()
	clustermat_minus<-c()
	numk_minus <- c()
	
	for(i in 1:(dimInput[2]))
	{
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

	for(i in 1:dim(input)[2])
	{
		simatrow1 <- c(simatrow1, rep(i,numk_plus[i]))
		simatcol1 <- c(simatcol1, c(1:numk_plus[i]))
		simatrow2 <- c(simatrow2, rep(i,numk_minus[i]))
		simatcol2 <- c(simatcol2, c(1:numk_minus[i]))
	}
	
	simatrow1 <<- simatrow1
	simatcol1 <<- simatcol1
	simatrow2 <<- simatrow2
	simatcol2 <<- simatcol2
		

	## seed generation
	rowcol_result1 <- list()
	rowcol_result2 <- list()
  seq_simat1<-as.list(c(1:length(simatrow1)))
  rowcol_result1<-lapply(seq_simat1, CandRowCol_plus)
  rowcol_result2<-lapply(seq_simat1, CandRowCol_minus)
  rowcol_result1 <- c(rowcol_result1, rowcol_result2)

  ptm <- proc.time()
  initial_rows <- list()
	initial_cols <- list()
	rownum <- c()
  colnum <- c()
  
	seed_cnt_tmp <- 0
	seed_tmp <- list()
	for(i in 1:length(rowcol_result1))
	{
		if (rowcol_result1[[i]]$col_cnt != 1)
		{
		  seed_cnt_tmp <- seed_cnt_tmp + 1
		  seed_tmp[[seed_cnt_tmp]] <- c(rowcol_result1[[i]]$candrow, 0, rowcol_result1[[i]]$candcol)
		}
	}
	
	## unique seed finding
	unique_seed <- unique(seed_tmp)
  if(length(unique_seed)!=0)
  {
    for (j in 1: length(unique_seed))
    {
      vec1 <- unique_seed[[j]]
      threspoint <- which(vec1==0)
      initial_rows[[j]] <- vec1[1:(threspoint-1)]
      initial_cols[[j]] <- vec1[(threspoint+1):length(vec1)]
      rownum[j] <- length(initial_rows[[j]])
			colnum[j] <- length(initial_cols[[j]])	
		 }
  }
  
  initial_rows <<- initial_rows
  initial_cols <<- initial_cols
  
  
 	seed_dim <- matrix(0, nrow=3, ncol=length(unique_seed))
	seed_dim[1,] <- rownum
	seed_dim[2,] <- colnum
	seed_dim[3,] <- rownum * colnum
	
	## stat of seeds
	seed_h <- c()
	seed_cor <- c()
	seq_stat0 <- as.list(c(1:(length(initial_rows))))
	stat_result_cor0 <- lapply(seq_stat0, PCC_VEC_0)
	
  seed_cor <- unlist(stat_result_cor0)
	summary0 <- c(length(seed_dim[1,]),mean(seed_dim[1,]), mean(seed_dim[2,]), mean(seed_dim[3,]), max(seed_dim[1,]), max(seed_dim[2,]), max(seed_dim[3,]), min(seed_dim[1,]), min(seed_dim[2,]), min(seed_dim[3,]), mean(seed_cor), max(seed_cor), min(seed_cor))

  ptm2 <- proc.time() - ptm

	## seed coverage calculation
  seed_covermat <- matrix(0, nrow=dim(input)[1], ncol=dim(input)[2])
  seed_num_tmp <- length(initial_rows)
  for(i in 1:length(initial_rows))
  {
    seed_covermat[initial_rows[[i]], initial_cols[[i]]]  <- seed_covermat[initial_rows[[i]], initial_cols[[i]]] + 1
  }
  
  row_coverage_tmp <- rowMeans(seed_covermat)
  col_coverage_tmp <- colMeans(seed_covermat)
  row_coverage_num <- length(which(row_coverage_tmp!=0))
  col_coverage_num <- length(which(col_coverage_tmp!=0))

  cover_cnt <- 0
  for(i in 1:dim(seed_covermat)[2])
  {
    cover_cnt <- cover_cnt + length(which(seed_covermat[,i]==0))  
  }
  coverage_ratio <- 1 - cover_cnt / (dim(input)[1] * dim(input)[2])
  seed_cover_vec <- c(row_coverage_num/dim(input)[1], col_coverage_num/dim(input)[2], coverage_ratio)
  
  ## seed expansion
  seq_merge <-as.list(c(1:length(initial_rows)))
  ptm <- proc.time()
  
	merge_result <-lapply(seq_merge, MergeMat_dynamic2, corstand, minrow, mincol, overlap_level)
  
  ptm3 <- proc.time() - ptm
  
	merge_num_tmp <- 0
	merged_rows_tmp <- list()
	merged_cols_tmp <- list()
	merged_rows_cols_tmp <- list()
	
  for(i in 1:length(merge_result))
	{
	 if(length(merge_result[[i]]$cor_vec)!=0)
	 {
  	 for(j in 1:length(merge_result[[i]]$cor_vec))
  	 {
  	   merge_num_tmp <- merge_num_tmp + 1
       merged_rows_tmp[[merge_num_tmp]] <- merge_result[[i]]$genes[[j]]
       merged_cols_tmp[[merge_num_tmp]] <- merge_result[[i]]$conds[[j]]
       merged_rows_cols_tmp[[merge_num_tmp]] <- c(merge_result[[i]]$genes[[j]],0, merge_result[[i]]$conds[[j]])       
     }
	 }
	}
	
	## unique bicluster finding
	unique_merged_rows_cols <- unique(merged_rows_cols_tmp)
	merged_rows <<- list()
	merged_cols <<- list()
	bilist <- list()
	gleng <- c()
	cleng <- c()

  if(length(unique_merged_rows_cols)!=0)
  {
    for (j in 1: length(unique_merged_rows_cols))
    {
      vec1 <- unique_merged_rows_cols[[j]]
      threspoint <- which(vec1 == 0)[1]
      merged_rows[[j]] <- vec1[1:(threspoint-1)]
      merged_cols[[j]] <- vec1[(threspoint+1):length(vec1)]
      bilist[[j]] <- input[merged_rows[[j]], merged_cols[[j]]]
      gleng[j] <- length(merged_rows[[j]])
      cleng[j] <- length(merged_cols[[j]])
    }
  }
  merged_rows <<- merged_rows
  merged_cols <<- merged_cols
  
  dimbi <- matrix(0,nrow=3,ncol=length(gleng))
  dimbi[1,] <- gleng
  dimbi[2,] <- cleng
  dimbi[3,] <- gleng * cleng
  
 	num_tmp <- cbind(dimbi[3,], dimbi[2,], dimbi[1,], c(1:length(dimbi[3,])))
 	
	if(length(dimbi[3,]) != 1)
	{
		num_order <- num_tmp[order(num_tmp[,1], num_tmp[,2], num_tmp[,3], decreasing=TRUE),][,4]
	}
	
	if(length(dimbi[3,]) == 1)
	{
		num_order <- c(1)
	}
	
	## stat of bicluster before overlap removal
  h <- c()
	cor <- c()
	seq_stat1 <- as.list(c(1:(length(merged_rows))))
	
  stat_result_cor1 <- lapply(seq_stat1, PCC_VEC_1)
	cor <- unlist(stat_result_cor1)
	
  if (length(gleng)<2)
	{
	  summary1 <- c()
	} else
	{
    summary1 <- c(length(dimbi[1,]),mean(dimbi[1,]), mean(dimbi[2,]), mean(dimbi[3,]), max(dimbi[1,]), max(dimbi[2,]), max(dimbi[3,]), min(dimbi[1,]), min(dimbi[2,]), min(dimbi[3,]), mean(cor), max(cor), min(cor))
	}
	## coverage of bicluster before overlapping removal
	
	covermat <- matrix(0, nrow=dim(input)[1], ncol=dim(input)[2])

  for(i in 1:length(merged_rows))
  {
    covermat[merged_rows[[i]], merged_cols[[i]]]  <- covermat[merged_rows[[i]], merged_cols[[i]]] + 1
  }

  row_coverage_tmp <- rowMeans(covermat)
  col_coverage_tmp <- colMeans(covermat)
  row_coverage_num <- length(which(row_coverage_tmp!=0))
  col_coverage_num <- length(which(col_coverage_tmp!=0))

  cover_cnt <- 0
  for(i in 1:dim(covermat)[2])
  {
    cover_cnt <- cover_cnt + length(which(covermat[,i]==0))
  }
  coverage_ratio <- 1 - cover_cnt / (dim(input)[1] * dim(input)[2])
  tmp_cover_vec <- c(row_coverage_num/dim(input)[1], col_coverage_num/dim(input)[2], coverage_ratio)

  ## remove overlapped bicluster
	order_vec <- num_order
	seq_duplicate<-as.list(c(1:(length(merged_rows)-1)))
	order_vec <- order(dimbi[3,], decreasing=TRUE)
	remove_duplicate_result <- lapply(seq_duplicate, CHECK_DUPLICATE, order_vec, overlap_level) 
	
  temp_chk <- c(rep(1,length(merged_rows)))
	for(i in 1:(length(merged_rows)-1))
	{
    temp_chk <- temp_chk * remove_duplicate_result[[i]]$chk_dup_vec 
	}
	true_index <- order_vec[which(temp_chk==1)]

	merged_rows_true <- list()
	merged_cols_true <- list()
	bilist_true <- list()
	cor_true <- c()
	gleng2 <- c()
	cleng2 <- c()
		
	for(i in 1:length(true_index))
	{
    merged_rows_true[[i]] <- merged_rows[[true_index[i]]]
    merged_cols_true[[i]] <- merged_cols[[true_index[i]]] 
    bilist_true[[i]] <- bilist[[true_index[i]]]
	}
  cor_true <- cor[true_index]
  gleng2 <- gleng[true_index]
  cleng2 <- cleng[true_index]
  dimbi_true <- dimbi[,true_index]
	 	
	## stat of biclusters after removal
	if (length(gleng2)<2)
	{
	 summary2 <- c()
   summary_total <- c()
	} else
  {
	 summary2 <- c(length(dimbi_true[1,]),mean(dimbi_true[1,]), mean(dimbi_true[2,]), mean(dimbi_true[3,]), max(dimbi_true[1,]), max(dimbi_true[2,]), max(dimbi_true[3,]), min(dimbi_true[1,]), min(dimbi_true[2,]), min(dimbi_true[3,]), mean(cor_true), max(cor_true), min(cor_true))
	 summary_total <- rbind(summary0, summary1, summary2)
	}
                                       
	## coverage of biclusters
	covermat <- matrix(0, nrow=dim(input)[1], ncol=dim(input)[2])

  for(i in 1:length(merged_rows_true))
  {
    covermat[merged_rows_true[[i]], merged_cols_true[[i]]]  <- covermat[merged_rows_true[[i]], merged_cols_true[[i]]] + 1
  }
  
  row_coverage_tmp <- rowMeans(covermat)
  col_coverage_tmp <- colMeans(covermat)
  row_coverage_num <- length(which(row_coverage_tmp!=0))
  col_coverage_num <- length(which(col_coverage_tmp!=0))

  cover_cnt <- 0
  for(i in 1:dim(covermat)[2])
  {
    cover_cnt <- cover_cnt + length(which(covermat[,i]==0))  
  }
  coverage_ratio <- 1 - cover_cnt / (dim(input)[1] * dim(input)[2])
  cover_vec <- c(row_coverage_num/dim(input)[1], col_coverage_num/dim(input)[2], coverage_ratio)
  
	
	
  ## write results to files
	
	which_max0 <- which(dimbi[3,]==max(dimbi[3,]))[1]
	which_min0 <- which(dimbi[3,]==min(dimbi[3,]))[1]
	bi_filename2 <- "summary_report.txt"
	cat("BICLIC_report", "\n\n",file=bi_filename2)
	cat("# of found biclusters: ", summary1[1], "\n",file=bi_filename2, append=TRUE)	
	cat("Maximum size is ", dimbi[1,which_max0], " * " , dimbi[2,which_max0], " = ", dimbi[3,which_max0], "\n",file=bi_filename2, append=TRUE)
	cat("Minimum size is ", dimbi[1,which_min0], " * " , dimbi[2,which_min0], " = ", dimbi[3,which_min0], "\n",file=bi_filename2, append=TRUE)
	cat("Average size is ", mean(dimbi[1,]), " * " , mean(dimbi[2,]), " = ", mean(dimbi[3,]), "\n",file=bi_filename2, append=TRUE)
	cat("Average PCC: ", mean(cor), "\n\n",file=bi_filename2, append=TRUE)
	
	which_max <- which(dimbi_true[3,]==max(dimbi_true[3,]))[1]
	which_min <- which(dimbi_true[3,]==min(dimbi_true[3,]))[1]
	cat("After overlapping removal", "\n\n",file=bi_filename2, append=TRUE)
	cat("# of found biclusters: ", summary2[1], "\n",file=bi_filename2, append=TRUE)	
	cat("Maximum size is ", dimbi_true[1,which_max], " * " , dimbi_true[2,which_max], " = ", dimbi_true[3,which_max], "\n",file=bi_filename2, append=TRUE)
	cat("Minimum size is ", dimbi_true[1,which_min], " * " , dimbi_true[2,which_min], " = ", dimbi_true[3,which_min], "\n",file=bi_filename2, append=TRUE)
	cat("Average size is ", mean(dimbi_true[1,]), " * " , mean(dimbi_true[2,]), " = ", mean(dimbi_true[3,]), "\n",file=bi_filename2, append=TRUE)
	cat("Average PCC: ", mean(cor_true), "\n",file=bi_filename2, append=TRUE)
	
	cat("# of found seeds: ", summary0[1], "\n",file=bi_filename2, append=TRUE)	
	cat("Maximum size is ", seed_dim[1,which(seed_dim[3,]==max(seed_dim[3,]))[1]], " * " , seed_dim[2,which(seed_dim[3,]==max(seed_dim[3,]))[1]], " = ", seed_dim[3,which(seed_dim[3,]==max(seed_dim[3,]))[1]], "\n",file=bi_filename2, append=TRUE)
	cat("Minimum size is ", seed_dim[1,which(seed_dim[3,]==min(seed_dim[3,]))[1]], " * " , seed_dim[2,which(seed_dim[3,]==min(seed_dim[3,]))[1]], " = ", seed_dim[3,which(seed_dim[3,]==min(seed_dim[3,]))[1]], "\n",file=bi_filename2, append=TRUE)
	cat("Average size is ", mean(seed_dim[1,]), " * " , mean(seed_dim[2,]), " = ", mean(seed_dim[3,]), "\n",file=bi_filename2, append=TRUE)
	cat("Average PCC: ", mean(seed_cor), "\n",file=bi_filename2, append=TRUE)
	
	gname2 <- as.character(gname)
	for(i in 1:summary2[1])
	{
		cat("Bicluster ", i,  "\nsize = ", dimbi_true[3,i], ",  of gene = ", dimbi_true[1,i], ",  of condition = ", dimbi_true[2,i], ", PCC = ", cor_true[i], "\n",file=bi_filename2, append=TRUE)
		cat("Genes are: ", gname2[merged_rows_true[[i]]], "\n", file=bi_filename2, append=TRUE)
		cat("Conditions are: ", merged_cols_true[[i]], "\n\n", file=bi_filename2, append=TRUE)	
	}
	

	
 return(list(seed_cover_vec=seed_cover_vec, cover_vec=cover_vec, tmp_cover_vec=tmp_cover_vec, gleng=gleng, cleng=cleng, gleng2=gleng2, cleng2=cleng2, cor_true=cor_true, summary_total=summary_total, bilist_true=bilist_true, merged_rows_true=merged_rows_true, merged_cols_true=merged_cols_true, initial_rows=initial_rows, initial_cols=initial_cols))
}

################################################
#CandRowCol_plus <- function(i, clustermat_plus, simatrow1, simatcol1)
CandRowCol_plus <- function(i)
{		
  #sizechk is minimum bicluster size threshold
  sizechk <- 2
  sizechk2 <- 2

	tmprow <- which(clustermat_plus[,simatrow1[i]]==simatcol1[i])

	if (length(tmprow) >= sizechk)
	{
		tmpcol <- which(rowVars(t(clustermat_plus[tmprow,]))==0)
		if (length(tmpcol) >= sizechk2)
		{
			col_cnt <- length(tmpcol)
		}
		if (length(tmpcol) < sizechk2)
		{
			col_cnt <- 1
		}
	}
	
	if (length(tmprow) < sizechk)
	{
		tmpcol <- 1
		col_cnt <- 1
	}
	
  candrow <- tmprow[order(tmprow,decreasing=F)]
  candcol <- tmpcol[order(tmpcol,decreasing=F)]
	return(list(candrow=candrow, candcol=candcol, col_cnt=col_cnt))
}

#CandRowCol_minus <- function(i, clustermap_minus, simatrow2, simatcol2)
CandRowCol_minus <- function(i)
{		
  #sizechk is minimum bicluster size threshold
  sizechk <- 2
  sizechk2 <- 2

	tmprow <- which(clustermat_minus[,simatrow2[i]]==simatcol2[i])

	if (length(tmprow) >= sizechk)
	{
		tmpcol <- which(rowVars(t(clustermat_minus[tmprow,]))==0)
		if (length(tmpcol) >= sizechk2)
		{
			col_cnt <- length(tmpcol)
		}
		if (length(tmpcol) < sizechk2)
		{
			col_cnt <- 1
		}
	}
	
	if (length(tmprow) < sizechk)
	{
		tmpcol <- 1
		col_cnt <- 1
	}
	
  candrow <- tmprow[order(tmprow,decreasing=F)]
  candcol <- tmpcol[order(tmpcol,decreasing=F)]
	return(list(candrow=candrow, candcol=candcol, col_cnt=col_cnt))
}

##################################################################

AddGene <- function(mgenes, mconds, corstand, input)
{
  new_mgenes <- c()
  seq_vec_row <- setdiff(c(1:dim(input)[1]),mgenes)
  seq_cor_row <- as.list(seq_vec_row)
  cor_result_row <- lapply(seq_cor_row, cor_mat_row_add, mgenes, mconds)
	
	if (length(cor_result_row)!=0)
	{
    before_order_cor_row <- unlist(cor_result_row[which(cor_result_row >= max(c(corstand/2,0.5)))])
    if(length(before_order_cor_row)!=0)
    {
    	before_order_row <- seq_vec_row[which(cor_result_row >= max(c(corstand/2,0.5)))]
    	after_order_cor_row <- before_order_cor_row[order(before_order_cor_row,decreasing=TRUE)]
    	after_order_row <- before_order_row[order(before_order_cor_row,decreasing=TRUE)]
      
       k <- 0
      loop_chk <- 0
      tmp_cor_vec <- c()
      
      candgenes_tmp <- mgenes
      while (loop_chk != -1)
      {
        k <- k + 1
        candgenes_tmp <- c(candgenes_tmp, after_order_row[k])
        candgenes <- candgenes_tmp[order(candgenes_tmp, decreasing=F)]
        m2 <- input[candgenes, mconds]
        tmp_cor <- cor_mat(m2)
        tmp_cor_vec[k] <- tmp_cor
    
        if ((k == length(after_order_row)) || (tmp_cor_vec[k] < corstand))
        {
          loop_chk <- -1
        }  
      }
      
      if (k == length(after_order_row))
      {
        row_thres_point <- k 
        cor_tmp <- tmp_cor_vec[row_thres_point]
      }
      
      if (tmp_cor_vec[k] < corstand)
      {
        row_thres_point <- k - 1
        cor_tmp <- tmp_cor_vec[row_thres_point]
      }
      
      if(row_thres_point!=0)
      {
        new_mgenes <- after_order_row[1:row_thres_point]
      }
    }
  }

  return (new_mgenes=new_mgenes)

}

##################################################################

AddCondition <- function(mgenes, mconds, corstand, input)
{
  new_mconds <- c()
  seq_vec_column <- setdiff(c(1:dim(input)[2]),mconds)
  seq_cor_column <- as.list(seq_vec_column)
  cor_result_column <- lapply(seq_cor_column, cor_mat_col_add, mgenes, mconds)
  if (length(cor_result_column)!=0)
  {
  	before_order_cor_column <- unlist(cor_result_column[which(cor_result_column >= max(c(corstand/2,0.5)))])
  	if(length(before_order_cor_column)!=0)
  	{
    	before_order_column <- seq_vec_column[which(cor_result_column >= max(c(corstand/2,0.5)))]
    	
    	after_order_cor_column <- before_order_cor_column[order(before_order_cor_column,decreasing=TRUE)]
    	after_order_column <- before_order_column[order(before_order_cor_column,decreasing=TRUE)]
    
      k <- 0
      loop_chk <- 0
      tmp_cor_vec <- c()
      
      candconds_tmp <- mconds
      while (loop_chk != -1)
      {
        k <- k + 1
        candconds_tmp <- c(candconds_tmp, after_order_column[k])
        candconds <- candconds_tmp[order(candconds_tmp, decreasing=F)]
        m2 <- input[mgenes, candconds]
        tmp_cor <- cor_mat(m2)
        tmp_cor_vec[k] <- tmp_cor
    
        if ((k == length(after_order_column)) || (tmp_cor_vec[k] < corstand))
        {
          loop_chk <- -1
        }  
      }
      
      if (k == length(after_order_column))
      {
        column_thres_point <- k 
        cor_tmp <- tmp_cor_vec[column_thres_point]
      }
      
      if (tmp_cor_vec[k] < corstand)
      {
        column_thres_point <- k - 1
        cor_tmp <- tmp_cor_vec[column_thres_point]
      }
      
      
      if(column_thres_point!=0)
      {
        new_mconds <- after_order_column[1:column_thres_point]
      }
    }
  }

  return (new_mconds=new_mconds)
}

######################################################################

MergeMat_dynamic2<- function(i, corstand, minrow, mincol, overlap_level)
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

	if(cor0 >= corstand)
	{
	  ### if column length is 2
		if (length(columnP)<mincol)
		{
      seq_vec_column <- setdiff(c(1:dimInput[2]),mconds)
	    seq_cor_column <- as.list(seq_vec_column)
      cor_result_column <- lapply(seq_cor_column, cor_mat_col_add, mgenes, mconds)
		
			before_order_cor_column <- unlist(cor_result_column[which(cor_result_column >= max(c(corstand/2,0.5)))])
			before_order_column <- seq_vec_column[which(cor_result_column >= max(c(corstand/2,0.5)))]
			
			after_order_cor_column <- before_order_cor_column[order(before_order_cor_column,decreasing=TRUE)]
			after_order_column <- before_order_column[order(before_order_cor_column,decreasing=TRUE)]
    
      num_to_col <- mincol - length(columnP)
      if (length(after_order_column) >= mincol -length(columnP))
      {
  	    mconds_tmp <- c(columnP, after_order_column[1:num_to_col])
  	    mconds <- mconds_tmp[order(mconds_tmp, decreasing=F)]
  	    mat0 <- input[originalG, mconds]
  	    cor0 <- cor_mat(mat0)
	    } else
	    {
	     cor0 <- 0
	    }
		}
	
		if (cor0 >= corstand)
		{
		  ## initial seed
      if ((length(mgenes)>= minrow) && (length(mconds) >= mincol))
      {
        total_num <- total_num + 1
  		  mgenes_list[[total_num]] <- mgenes
  		  mconds_list[[total_num]] <- mconds
  		  mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 0, mconds_list[[total_num]])
		  }
		  
		  ## condition expanded bicluster
		  cond_max_result <- AddCondition(mgenes, mconds, corstand, input)
      gene_tmp <- mgenes
      cond_tmp <- c(mconds,cond_max_result)
      if ((length(gene_tmp)>= minrow) && (length(cond_tmp) >= mincol))
      {
        total_num <- total_num + 1
        mgenes_list[[total_num]] <- gene_tmp[order(gene_tmp, decreasing=F)]
        mconds_list[[total_num]] <- cond_tmp[order(cond_tmp, decreasing=F)]
        mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 0, mconds_list[[total_num]])
      }

      ## gene expanded bicluster
      gene_max_result <- AddGene(mgenes, mconds, corstand, input)
		  gene_tmp <- c(mgenes,gene_max_result)
      cond_tmp <- mconds
      if ((length(gene_tmp)>= minrow) && (length(cond_tmp) >= mincol))
      {
        total_num <- total_num + 1
        mgenes_list[[total_num]] <- gene_tmp[order(gene_tmp, decreasing=F)]
        mconds_list[[total_num]] <- cond_tmp[order(cond_tmp, decreasing=F)]
        mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 0, mconds_list[[total_num]])
      }
      
      ## backward eliminated bicluster
      gene_tmp <- c(mgenes, gene_max_result)
      cond_tmp <- c(mconds, cond_max_result)
      gene_tmp <- gene_tmp[order(gene_tmp,decreasing=F)]
      cond_tmp <- cond_tmp[order(cond_tmp,decreasing=F)]
      tmp_mat <- input[gene_tmp, cond_tmp]
      
      gene_tmp2 <- gene_tmp
      cond_tmp2 <- cond_tmp    
      tmp_mat <- input[gene_tmp, cond_tmp]
      tmp_mat2 <- tmp_mat
      cor_initial <- cor_mat(tmp_mat)
      if (cor_initial >= 0.999)
      {
       cor_initial <- 1
      }
      before_cor <- cor_initial
      sentinel <- 0
      roop_cnt <- 0
      
      if(cor_initial < corstand)
      {
        while (sentinel != -1)
        {
          roop_cnt <- roop_cnt + 1
          if (length(gene_tmp2) > 2)
          {
            gene_cor_vec <- numeric(length(gene_tmp2))
            cond_cor_vec <- numeric(length(cond_tmp2))
            for(j in 1:length(gene_tmp2))
            {
              gene_cor_vec[j] <- cor_mat(input[gene_tmp2[-j], cond_tmp2])
            }
            
            gene_diff_vec <- gene_cor_vec - before_cor
            cand_gene_remove <- which(gene_diff_vec > 0)
            if((length(gene_tmp2) - length(cand_gene_remove)) < 2)
            {
              cand_gene_remove2_tmp <- gene_diff_vec[cand_gene_remove]
              cand_gene_remove2 <- cand_gene_remove[order(cand_gene_remove2_tmp, decreasing=T)]
              gene_cor_vec2 <- numeric(length(cand_gene_remove2))
          
              for(k in 1:(length(gene_tmp2)-2))
              {
                gene_cor_vec2[k] <- cor_mat(input[gene_tmp2[-cand_gene_remove2[k]], cond_tmp2])
              }
              after_cor_gene <- gene_cor_vec2[1]
              diff_cor_gene <- after_cor_gene - before_cor
              cand_gene_remove <- cand_gene_remove2[1]
              
            } else
            {
              after_cor_gene <- cor_mat(input[gene_tmp2[-cand_gene_remove], cond_tmp2]) 
              diff_cor_gene <- after_cor_gene - before_cor
            }
            
          } else
          {
            diff_cor_gene <- -1
            after_cor_gene <- 0  
          }
          
          for(j in 1:length(cond_tmp))
          {
            cond_cor_vec[j] <- cor_mat(input[gene_tmp2, cond_tmp2[-j]])
          }
          cond_diff_vec <- cond_cor_vec - before_cor
          after_cor_cond <- max(cond_cor_vec)
          diff_cor_cond <- max(cond_diff_vec)
          cand_cond_remove <- which(cond_diff_vec==max(cond_diff_vec))
          
          if ((diff_cor_gene <= 0) && (diff_cor_cond <= 0))
          {
            sentinel <- -1
          } else if (((diff_cor_gene > 0) && (diff_cor_cond <= 0)) || (after_cor_gene > after_cor_cond))
          {
            gene_tmp3 <- gene_tmp2
            gene_tmp2 <- gene_tmp2[-cand_gene_remove]

            tmp_mat2 <- input[gene_tmp2, cond_tmp2]
            before_cor <- cor_mat(tmp_mat2)
            if (before_cor >= corstand)
            {
              cand_gene_remove2_tmp <- gene_diff_vec[cand_gene_remove]
              cand_gene_remove2 <- cand_gene_remove[order(cand_gene_remove2_tmp, decreasing=T)]
              gene_cor_vec2 <- numeric(length(cand_gene_remove2))
          
              for(k in 1:length(cand_gene_remove2))
              {
                gene_cor_vec2[k] <- cor_mat(input[gene_tmp3[-cand_gene_remove2[1:k]], cond_tmp2])
              }
              gene_tmp2 <- gene_tmp3[-cand_gene_remove2[1:which(gene_cor_vec2 >= corstand)[1]]]
              before_cor <- cor_mat(input[gene_tmp2,cond_tmp2]) 
              sentinel <- -1
            }
          } else if (((diff_cor_gene <= 0) && (diff_cor_cond > 0)) || (after_cor_gene < after_cor_cond))
          {
            cond_tmp3 <- cond_tmp2
            cond_tmp2 <- cond_tmp3[-cand_cond_remove]
            tmp_mat2 <- input[gene_tmp2, cond_tmp2]
            before_cor <- cor_mat(tmp_mat2)
            if (before_cor >= corstand)
            {
              sentinel <- -1
            }
          }
          
        }
      } else
      {
        if ((length(gene_tmp)>= minrow) && (length(cond_tmp) >= mincol))
        {
          total_num <- total_num + 1
          mgenes_list[[total_num]] <- gene_tmp[order(gene_tmp, decreasing=F)]
          mconds_list[[total_num]] <- cond_tmp[order(cond_tmp, decreasing=F)]
          mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 0, mconds_list[[total_num]])
        }
      }
      
      if (((length(gene_tmp2)>= minrow) && (length(cond_tmp2) >= mincol)) && (before_cor > corstand))
      {
        total_num <- total_num + 1
        mgenes_list[[total_num]] <- gene_tmp2[order(gene_tmp2, decreasing=F)]
        mconds_list[[total_num]] <- cond_tmp2[order(cond_tmp2, decreasing=F)]
        mgenes_mconds_list[[total_num]] <- c(mgenes_list[[total_num]], 0, mconds_list[[total_num]])
      }
      
      ## unique gene
      unique_mgenes_mconds <- unique(mgenes_mconds_list)
      if(length(unique_mgenes_mconds)!=0)
      {
        gleng <- c()
        cleng <- c()
        for (j in 1: length(unique_mgenes_mconds))
        {
          vec1 <- unique_mgenes_mconds[[j]]
          threspoint <- which(vec1 == 0)[1]
          unique_genes[[j]] <- vec1[1:(threspoint-1)]
          unique_conds[[j]] <- vec1[(threspoint+1):length(vec1)]
          gleng[j] <- length(unique_genes[[j]])
          cleng[j] <- length(unique_conds[[j]])
          cor_vec_tmp[j] <- cor_mat(input[unique_genes[[j]],unique_conds[[j]]])
        }
        
        dimbi <- matrix(0,nrow=3,ncol=length(gleng))
        dimbi[1,] <- gleng
        dimbi[2,] <- cleng
        dimbi[3,] <- gleng * cleng
        
       	num_tmp <- cbind(dimbi[3,], dimbi[2,], dimbi[1,], c(1:length(dimbi[3,])))
       	
      	if(length(dimbi[3,]) != 1)
      	{
      		num_order <- num_tmp[order(num_tmp[,1], num_tmp[,2], num_tmp[,3], decreasing=TRUE),][,4]
      	}
      	
      	if(length(dimbi[3,]) == 1)
      	{
      		num_order <- c(1)
      	}
      
        ## overlap removal
        if (length(gleng)!=1)
        {
          overlap_result <- CHECK_DUPLICATE2(num_order, unique_genes, unique_conds, overlap_level)
          temp_chk <- overlap_result$chk_dup_vec
          if(length(temp_chk)!=0)
          {
            true_index <- num_order[which(temp_chk==1)]
    	      for(j in 1:length(true_index))
    	      {
    	       genes[[j]] <- unique_genes[[true_index[j]]]
    	       conds[[j]] <- unique_conds[[true_index[j]]]
    	       cor_vec[j] <- cor_vec_tmp[true_index[j]]
    	       
    	      }
  	      }
	      }
	      
	      if (length(gleng)==1)
        {
          genes <- unique_genes
	        conds <- unique_conds
	        cor_vec <- cor_vec_tmp
        }
      }   
    }
  }
  
  return(list(genes=genes, conds=conds, cor_vec=cor_vec))
}

#########################################
cor_mat <- function(mat)
{
	adjustzero <- intersect(which(rowMeans(mat)==0), which(rowVars(mat)==0)) 
	if (length(adjustzero) != 0)
	{
		mat[adjustzero,] <- c(rep(0.001, dim(mat)[2]))
	}

	adjustwhere <- which(rowVars(mat)==0)
	if (length(adjustwhere)!=0)
	{
		for(jj in 1:length(adjustwhere))
		{
			mat[adjustwhere[jj],][1] <- mat[adjustwhere[jj],][1] + abs(mat[adjustwhere[jj],][1]/10000) 
		}
	}

	cormat <- cor(t(mat))
	homoRow <- mean(cormat[lower.tri(cormat)])

	return (homoRow=homoRow)
}



cor_mat_row <- function(jth, mat0, columnP, corstand2)
{
	mat0 <- colMeans(mat0)
	mat <- rbind(mat0, input[jth,columnP])
	cor_chk <- 0
	corchk2 <- 1
	

	homoRow <- cor_mat(mat)
	if (homoRow >= corstand2)
	{
		cor_chk <- 1
	}
	
	return(cor_chk=cor_chk)
}

cor_mat_row_add <- function(jth, row_tmp, column_tmp)
{
  temp_row <- c(row_tmp, jth)
  merged_mat <- input[temp_row, column_tmp]
	homoRow <- cor_mat(merged_mat)
	
	return(homoRow=homoRow)
}

cor_mat_col_add <- function(jth, row_tmp, column_tmp)
{
  temp_col <- c(column_tmp, jth)
	merged_mat <- input[row_tmp, temp_col]
	
	homoRow <- cor_mat(merged_mat)
	
	return(homoRow=homoRow)
}


#############################################
SDbasedClustering <- function(i)
{	
	temp.col <- input[,i]
	sd_max <- sd(temp.col)
	
	if (length(table(temp.col))==1)
  {
  	clvec <- c(rep(1, length(temp.col)))
  	clvec2 <- c(rep(1, length(temp.col)))
  }
  
  if (length(table(temp.col))!=1)
  {
  	cl <- c(1)
  	cnt <- 1
    indexed <- cbind(temp.col, c(1:length(temp.col)))
    indexedo <- indexed[order(indexed[,1]),]

    sentinel <- 0
    j <- 1
    
    while(sentinel != -1)
    {
    	k <- 1
    	sentinel2 <- 0
    	
    	if (j >= length(indexedo[,1]))
    	{
    		sentinel <- -1
    		sentinel2 <- -1
    	}
    	
    	if ((j + k) == length(indexedo[,1]))
  		{
  			cl <- c(cl, rep(cnt, (k+1)))
  			sentinel2 <- -1
  			sentinel <- -1
  			pass <- 1
  		}
  		if ((j + k) < length(indexedo[,1]))
  		{
	    	cl_target1 <- indexedo[j:(j+k),1]
	    	sd_cl <- sd(cl_target1)
	    	if (sd_cl > sd_max)
	    	{
	    		pass2 <- 0
	    		if(j == 1)
	    		{
	    			cnt <- cnt + 1
	    			j <- j + 1
	    			pass2 <- 1
	    			sentinel2 <- -1
	    		}
	    		if ((j!=1) && (pass2 == 0))
	    		{
	    			cl <- c(cl, cnt)
	    			cnt <- cnt + 1
	    			j <- j + 1
	    			pass2 <- 1
	    			sentinel2 <- -1
	    		}
	    	}
	    }
    	
    	while(sentinel2 != -1)
    	{
    		pass <- 0
    		k <- k + 1
    		if ((j + k) == length(indexedo[,1]))
    		{
    			cl <- c(cl, rep(cnt, (k+1)))
    			sentinel2 <- -1
    			sentinel <- -1
    			pass <- 1
    		}
    		cl_target2 <- indexedo[j:(j+k),1]
    		
    		sd_cl2 <- sd(cl_target2)
    		
    		if (sd_cl2 <= sd_cl)
    		{
    			sd_cl <- sd_cl2
    			pass <- 1
    		}
    		if ((sd_cl2 > sd_cl) && (pass == 0))
    		{
    			pass2 <- 0
    			if (j == 1)
    			{
  					cl <- c(cl, rep(cnt,(k-1)))
  					cnt <- cnt + 1
  					j <- j + k 
  					pass2 <- 1
  				}
  				if ((j > 1) && (pass2 == 0))
    			{
  					cl <- c(cl, rep(cnt, k))
  					cnt <- cnt + 1
  					j <- j + k 
  				}
  				sentinel2 <- -1
    		}
    	}
    }
    
    indexedo_tmp <- cbind(indexedo, cl)
  	indexedo_tmpo <- indexedo_tmp[order(indexedo_tmp[,2]),]
  	clvec <- indexedo_tmpo[,3]
  }

  return(list(clvec=clvec))
}

##############################################################
SDbasedClustering_minus <- function(i)
{	
	temp.col <- input[,i]
	sd_max <- sd(temp.col)
	
	if (length(table(temp.col))==1)
  {
  	clvec <- c(rep(1, length(temp.col)))
  	clvec2 <- c(rep(1, length(temp.col)))
  }
  
  if (length(table(temp.col))!=1)
  {
  	cl <- c(1)
  	cnt <- 1
    indexed <- cbind(temp.col, c(1:length(temp.col)))
    indexedo <- indexed[order(indexed[,1], decreasing=TRUE),]

    sentinel <- 0
    j <- 1
    
    while(sentinel != -1)
    {
    	k <- 1
    	sentinel2 <- 0
    	
    	if (j >= length(indexedo[,1]))
    	{
    		sentinel <- -1
    		sentinel2 <- -1
    	}
    	
    	if ((j + k) == length(indexedo[,1]))
  		{
  			cl <- c(cl, rep(cnt, (k+1)))
  			sentinel2 <- -1
  			sentinel <- -1
  			pass <- 1
  		}
  		if ((j + k) < length(indexedo[,1]))
  		{
	    	cl_target1 <- indexedo[j:(j+k),1]
	    	sd_cl <- sd(cl_target1)
	    	if (sd_cl > sd_max)
	    	{
	    		pass2 <- 0
	    		if(j == 1)
	    		{
	    			cnt <- cnt + 1
	    			j <- j + 1
	    			pass2 <- 1
	    			sentinel2 <- -1
	    		}
	    		if ((j!=1) && (pass2 == 0))
	    		{
	    			cl <- c(cl, cnt)
	    			cnt <- cnt + 1
	    			j <- j + 1
	    			pass2 <- 1
	    			sentinel2 <- -1
	    		}
	    	}
	    }
    	
    	while(sentinel2 != -1)
    	{
    		pass <- 0
    		k <- k + 1
    		if ((j + k) == length(indexedo[,1]))
    		{
    			cl <- c(cl, rep(cnt, (k+1)))
    			sentinel2 <- -1
    			sentinel <- -1
    			pass <- 1
    		}
    		cl_target2 <- indexedo[j:(j+k),1]
    		
    		sd_cl2 <- sd(cl_target2)
    		
    		if (sd_cl2 <= sd_cl)
    		{
    			sd_cl <- sd_cl2
    			pass <- 1
    		}
    		if ((sd_cl2 > sd_cl) && (pass == 0))
    		{
    			pass2 <- 0
    			if (j == 1)
    			{
  					cl <- c(cl, rep(cnt,(k-1)))
  					cnt <- cnt + 1
  					j <- j + k 
  					pass2 <- 1
  				}
  				if ((j > 1) && (pass2 == 0))
    			{
  					cl <- c(cl, rep(cnt, k))
  					cnt <- cnt + 1
  					j <- j + k 
  				}
  				sentinel2 <- -1
    		}
    	}
    }
    
    indexedo_tmp <- cbind(indexedo, cl)
  	indexedo_tmpo <- indexedo_tmp[order(indexedo_tmp[,2], decreasing=FALSE),]
  	clvec <- indexedo_tmpo[,3]
  }

  return(list(clvec=clvec))
}
  

DrawPattern <- function(i, filename3, corstand)
{
	bi_filename3 <- paste("./result/",filename3, sep="")
	bi_filename3 <- paste(bi_filename3,"_", sep="")
	bi_filename3 <- paste(bi_filename3,corstand, sep="")
	bi_filename3 <- paste(bi_filename3,"_bicluster_pattern#", sep="")
	bi_filename3 <- paste(bi_filename3,i, sep="")
	bi_filename3 <- paste(bi_filename3,".png", sep="")
	
	patternname <- paste("bicluster_pattern #", i,sep="")
	mat <- input[merged_rows_true[[i]], merged_cols_true[[i]]]
	x1 <- 1
	x2 <- dim(mat)[2]
	y1 <- min(mat)
	y2 <- max(mat)
  
  png(file=bi_filename3)		
	matplot(c(x1,x2), c(y1,y2), type= "n", xlab = "Conditions", ylab = "Expression level", main = patternname)
		
	for(j in 1:dim(mat)[1])
	{
		matlines(as.numeric(mat[j,]),col = 1)
	}
	
	dev.off()
	dev.next()
}


CHECK_DUPLICATE <- function(i, order_vec, overlap_level)
{
  chk_dup_vec <- c(rep(1,length(merged_rows)))
	row1 <- merged_rows[[order_vec[i]]]
  col1 <- merged_cols[[order_vec[i]]]
  area_vec <- c(rep(0, length(chk_dup_vec)))
  area1 <- length(row1) * length(col1)
  for(j in 1:length(merged_rows))
  {
    if (i!=j)
    {
      row2 <- merged_rows[[order_vec[j]]]
      col2 <- merged_cols[[order_vec[j]]]
     
      row_inter <- intersect(row1, row2)
      col_inter <- intersect(col1, col2)
      area_inter <- length(row_inter) * length(col_inter)
      area2 <- length(row2) * length(col2)
      
      if ((area1 > area2) || ((area1 == area2) && (i < j)))
      { 
        area_frac <- area_inter / area2
        area_vec[j] <- area_frac
        if(area_frac >= overlap_level)
        {
          chk_dup_vec[j] <- 0
        }
      } 
    }
   }
   
   return(list(chk_dup_vec=chk_dup_vec, area_vec=area_vec))
 }
 
CHECK_DUPLICATE2 <- function(order_vec, merged_rows_func, merged_cols_func, overlap_level)
{
  cnt <- 0
  chk_dup_vec <- c(rep(1,length(merged_rows_func)))
  for(i in 1:(length(merged_rows_func)-1))
  {
    if(chk_dup_vec[i]!=0)
    {
    	row1 <- merged_rows_func[[order_vec[i]]]
      col1 <- merged_cols_func[[order_vec[i]]]
      area_vec <- c(rep(0, length(chk_dup_vec)))
      area1 <- length(row1) * length(col1)
      for(j in 1:length(merged_rows_func))
      {
        if((chk_dup_vec[j] != 0) && (i != j))
        {
          cnt <- cnt + 1
          row2 <- merged_rows_func[[order_vec[j]]]
          col2 <- merged_cols_func[[order_vec[j]]]
         
          row_inter <- intersect(row1, row2)
          col_inter <- intersect(col1, col2)
          area_inter <- length(row_inter) * length(col_inter)
          area2 <- length(row2) * length(col2)
          
          if ((area1 > area2) || ((area1 == area2) && (i < j)))
          {
            area_frac <- area_inter / area2
            area_vec[j] <- area_frac
            if(area_frac >= overlap_level)
            {
              chk_dup_vec[j] <- 0
            }
          } 
         }
       }
     }
   }
     
   return(list(chk_dup_vec=chk_dup_vec, area_vec=area_vec))
 }
      

 #PCC_VEC_0 <- function(i, input, rows, cols)
 PCC_VEC_0 <- function(i)
 {
  cors <- c()
  mat <- input[initial_rows[[i]], initial_cols[[i]]]
  cors <- cor_mat(mat)

  return(cors=cors)
 }
 
 #PCC_VEC_1 <- function(i, input, merged_rows, merged_cols)
 PCC_VEC_1 <- function(i)
 {
  cors <- c()
  mat <- input[merged_rows[[i]], merged_cols[[i]]]
  cors <- cor_mat(mat)

  return(cors=cors)
 }
 
  	



