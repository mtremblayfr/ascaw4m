asca_w4m <- function(datamatrix, samplemetadata, factors, variablemetadata, threshold, scaling="none", nPerm)
{
  ## Transpose
#  datamatrix <- t(datamatrix)

  # Check sample ID's
  rownames(datamatrix) <- make.names(rownames(datamatrix), unique = TRUE)
  colnames(datamatrix) <- make.names(colnames(datamatrix), unique = TRUE)
  rownames(samplemetadata) <- make.names(rownames(samplemetadata), unique = TRUE)
  rownames(variablemetadata) <- make.names(rownames(variablemetadata), unique = TRUE)
  
  if(!identical(rownames(datamatrix), rownames(samplemetadata))) 
  {
    if(identical(sort(rownames(datamatrix)), sort(rownames(samplemetadata)))) 
    {
      cat("\n\nMessage: Re-ordering dataMatrix sample names to match sampleMetadata\n")
      datamatrix <- datamatrix[rownames(samplemetadata), , drop = FALSE]
      stopifnot(identical(sort(rownames(datamatrix)), sort(rownames(samplemetadata))))
    }else {
      
      cat("\n\nStop: The sample names of dataMatrix and sampleMetadata do not match:\n")
      print(cbind.data.frame(indice = 1:nrow(datamatrix),
                             dataMatrix=rownames(datamatrix),
                             sampleMetadata=rownames(samplemetadata))[rownames(datamatrix) != rownames(samplemetadata), , drop = FALSE])
    }
  }
  
  # Check feature ID's
  if(!identical(colnames(datamatrix), rownames(variablemetadata))) 
  {
   if(identical(sort(colnames(datamatrix)), sort(rownames(variablemetadata)))) 
   {
      cat("\n\nMessage: Re-ordering dataMatrix variable names to match variableMetadata\n")
      datamatrix <- datamatrix[, rownames(variablemetadata), drop = FALSE]
      stopifnot(identical(sort(colnames(datamatrix)), sort(rownames(variablemetadata))))
    }else {
      cat("\n\nStop: The variable names of dataMatrix and variableMetadata do not match:\n")
      print(cbind.data.frame(indice = 1:ncol(datamatrix),
                             dataMatrix=colnames(datamatrix),
                             variableMetadata=rownames(variablemetadata))[colnames(datamatrix) != rownames(variablemetadata), , drop = FALSE])
    }
   }

  # Design
  design <- data.matrix(samplemetadata[, colnames(samplemetadata) %in% factors])
  
  # Scaling if scaling!=none
  datamatrix <- prep(datamatrix, scaling)

  # Computation of the A-SCA model
  data.asca <- ASCA.Calculate_w4m(datamatrix, design, scaling=scaling)

  # Permutation test
  data.asca.permutation <- ASCA.DoPermutationTest(data.asca[[1]], perm=nPerm)
  p <- c(data.asca.permutation, 0)
  

  # % of explained variance
  ssq <- (data.asca[[2]]$summary.ssq)
  ssq <- cbind(round(rbind(ssq[2], ssq[3],ssq[4],ssq[5])*100, 2), p)
  rownames(ssq) <- c(factors[1], factors[2], "Interaction", "Residuals")
  colnames(ssq) <- c("% of explained variance", "Permutation p-value")

  # Add Scores and loadings at the end of meatadata files
  noms <- colnames(samplemetadata)
  samplemetadata <- cbind(samplemetadata, (data.asca[[1]]$'1'$means.matrix + data.asca[[1]]$remainder) %*% data.asca[[1]]$'1'$svd$v[, 1:2], 
                  (data.asca[[1]]$'2'$means.matrix + data.asca[[1]]$remainder) %*% data.asca[[1]]$'2'$svd$v[, 1:2], 
                  (data.asca[[1]]$'12'$means.matrix + data.asca[[1]]$remainder) %*% data.asca[[1]]$'12'$svd$v[, 1:2])
  colnames(samplemetadata) <- c(noms, paste(factors[1],"XSCOR-p1", sep="_"), paste(factors[1],"XSCOR-p2", sep="_"), 
                        paste(factors[2],"XSCOR-p1", sep="_"), paste(factors[2],"XSCOR-p2", sep="_"), 
                        "Interact_XSCOR-p1", "Interact_XSCOR-p2")
  
  noms <- colnames(variablemetadata)
  variablemetadata <- cbind(variablemetadata, data.asca[[1]]$'1'$svd$v[, 1:2], data.asca[[1]]$'2'$svd$v[, 1:2], data.asca[[1]]$'12'$svd$v[, 1:2])
  colnames(variablemetadata) <- c(noms, paste(factors[1],"XLOAD-p1", sep="_"), paste(factors[1],"XLOAD-p2", sep="_"), 
                        paste(factors[2],"XLOAD-p1", sep="_"), paste(factors[2],"XLOAD-p2", sep="_"), 
                        "Interact_XLOAD-p1", "Interact_XLOAD-p2")
  
  l <- list(data.asca[[1]], data.asca.permutation, ssq, samplemetadata, variablemetadata)
  names(l) <- c("ASCA","p-values", "ssq", "samplemetadata", "variablemetadata")
  return(l)
}
