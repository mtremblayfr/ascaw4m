#!/usr/bin/env Rscript

###################################################################################################
#
# MetStaT ASCA.calculate function
#
#
# R-Package: MetStaT
#
# Version: 1.0
#
# Author (asca.calculate): Tim Dorscheidt
# Author (wrapper & .r adaptation for workflow4metabolomics.org): M. Tremblay-Franco & Y. Guitton #
# 
# Expected parameters from the commandline
# input files:
#             dataMatrix
#             sampleMetadata
#             variableMetadata
# params:
#             Factors (Factor1 & Factor2)
#             scaling
#             Number of permutations
#             Significance threshold
# output files:
#             sampleMetadata
#             variableMetadata
#			  Graphical outputs
#			  Information text
###################################################################################################
pkgs=c("MetStaT","batch","pcaMethods")
for(pkg in pkgs) {
  suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))
  cat(pkg,"\t",as.character(packageVersion(pkg)),"\n",sep="")
}


listArguments = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects

#Redirect all stdout to the log file
sink(listArguments$information)

# ----- PACKAGE -----
cat("\tPACKAGE INFO\n")
sessionInfo()

source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}

#load asca_w4m function
source_local("asca_w4m.R")
source_local("ASCA.Calculate_w4m.R")
source_local("ASCA.PlotScoresPerLevel_w4m.R")
print("first loadings OK")

## libraries
##----------

cat('\n\nRunning ASCA.calculate\n');
options(warn=-1);
#remove rgl warning
options(rgl.useNULL = TRUE);


## constants
##----------

modNamC <- "asca" ## module name

topEnvC <- environment()
flgC <- "\n"

## functions
##----------For manual input of function
##--end function

flgF <- function(tesC,
                 envC = topEnvC,
                 txtC = NA) { ## management of warning and error messages

    tesL <- eval(parse(text = tesC), envir = envC)

    if(!tesL) {

        sink(NULL)
        stpTxtC <- ifelse(is.na(txtC),
                          paste0(tesC, " is FALSE"),
                          txtC)

        stop(stpTxtC,
             call. = FALSE)

    }

} ## flgF


## log file
##---------
cat("\nStart of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")


## arguments
##----------
## loading files and checks
xMN <- t(as.matrix(read.table(listArguments[["dataMatrix_in"]],
                              check.names = FALSE,
                              header = TRUE,
                              row.names = 1,
                              sep = "\t")))

samDF <- read.table(listArguments[["sampleMetadata_in"]],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
					sep = "\t")

varDF <- read.table(listArguments[["variableMetadata_in"]],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
					sep = "\t")

result <- asca_w4m(xMN, samDF, c(listArguments[["factor1"]],listArguments[["factor2"]]), varDF, as.numeric(listArguments[["threshold"]]), 
                   scaling=listArguments[["scaling"]], listArguments[["nPerm"]])


##saving

if (exists("result")) {
	## writing output files
	cat("\n\nWriting output files\n\n");
	write.table(result[[4]],
			file = listArguments$sampleMetadata_out,
			quote = FALSE,
			row.names = TRUE,
			sep = "\t")

	write.table(result[[5]],
			file = listArguments$variableMetadata_out,
			quote = FALSE,
			row.names = TRUE,
			sep = "\t")

	# Graphical display for each significant parameter
	print(result[[3]])
	cat("\n p-value of Residuals must not be taken into account\n")
	
	if (any(result[[2]] < as.numeric(listArguments[["threshold"]])))
	{
		data.asca.permutation <- result[[2]]
		design <- data.matrix(samDF[, colnames(samDF) %in% c(listArguments[["factor1"]],listArguments[["factor2"]])])
		
		pdf(listArguments$figure, onefile=TRUE)
		par(mfrow=c(1,3))
		if (data.asca.permutation[1] < as.numeric(listArguments[["threshold"]]))
		{
			eigenvalues <- data.frame(1:length(unique(design[,1])), result[[1]]$'1'$svd$var.explained[1:length(unique(design[,1]))])
			colnames(eigenvalues) <- c("PC", "explainedVariance")
			barplot(eigenvalues[,2], names.arg=eigenvalues[,1], ylab="% of explained variance", xlab="Principal component")
			noms <- levels(as.factor(samDF[, listArguments$factor1]))
			ASCA.PlotScoresPerLevel_w4m(result[[1]], ee="1", interaction=0, factorName=listArguments$factor1, factorModalite=noms)
			Date.loadings <- data.matrix(result[[5]][,2:3])
			Date.loadings.leverage <- diag(Date.loadings%*%t(Date.loadings))
			names(Date.loadings.leverage) <- colnames(xMN)
			Date.loadings.leverage <- sort(Date.loadings.leverage, decreasing=TRUE)
			barplot(Date.loadings.leverage[Date.loadings.leverage > 0.001], main="PC1 loadings")
		}
		if (data.asca.permutation[2] < as.numeric(listArguments[["threshold"]]))
		{
			eigenvalues <- data.frame(1:length(unique(design[,2])), result[[1]]$'2'$svd$var.explained[1:length(unique(design[,2]))])
			colnames(eigenvalues) <- c("PC", "explainedVariance")
			barplot(eigenvalues[,2], names.arg=eigenvalues[,1], ylab="% of explained variance", xlab="Principal component")    
			noms <- levels(as.factor(samDF[, listArguments$factor2]))
			ASCA.PlotScoresPerLevel_w4m(result[[1]], ee="2", interaction=0, factorName=listArguments$factor2, factorModalite=noms)
			Date.loadings <- data.matrix(result[[5]][,4:5])
			Date.loadings.leverage <- diag(Date.loadings%*%t(Date.loadings))
			names(Date.loadings.leverage) <- colnames(xMN)
			Date.loadings.leverage <- sort(Date.loadings.leverage, decreasing=TRUE)
			barplot(Date.loadings.leverage[Date.loadings.leverage > 0.001], main="PC1 loadings")
		}
  	if (data.asca.permutation[3] < as.numeric(listArguments[["threshold"]]))
  	{
  	  eigenvalues <- data.frame(1:(length(unique(design[,1]))*length(unique(design[,2]))), result[[1]]$'12'$svd$var.explained[1:(length(unique(design[,1]))*length(unique(design[,2])))])
  	  colnames(eigenvalues) <- c("PC", "explainedVariance")
  	  barplot(eigenvalues[,2], names.arg=eigenvalues[,1], ylab="% of explained variance", xlab="Principal component")
  	  noms1 <- data.matrix(levels(as.factor(samDF[, listArguments$factor1])))
  	  noms2 <- data.matrix(levels(as.factor(samDF[, listArguments$factor2])))
  	  noms <- apply(noms1, 1, FUN=function(x){paste(x, "-", noms2, sep="")})
  	  noms <- apply(noms, 1, FUN=function(x){c(noms)})
  	  ASCA.PlotScoresPerLevel_w4m(result[[1]], ee="12", interaction=1, factorModalite=noms[,1])
  	  Date.loadings <- data.matrix(result[[5]][,6:7])
  	  Date.loadings.leverage <- diag(Date.loadings%*%t(Date.loadings))
  	  names(Date.loadings.leverage) <- colnames(xMN)
  	  Date.loadings.leverage <- sort(Date.loadings.leverage, decreasing=TRUE)
  	  barplot(Date.loadings.leverage[Date.loadings.leverage > 0.001], main="PC1 loadings")
  	}
    dev.off()
	}

	tryCatch({
	save(result, file="asca.RData");
	}, warning = function(w) {
	print(paste("Warning: ", w));
	}, error = function(err) {
	stop(paste("ERROR saving result RData object:", err));
	});
}

## ending
##-------

cat("\nEnd of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")

sink()

# options(stringsAsFactors = strAsFacL)


rm(list = ls())
