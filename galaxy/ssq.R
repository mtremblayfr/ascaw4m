library(lmdme)
library(MetStaT)

## Data : attention standardiser au prealable
data <- read.table("E:/PROJETS/Asca_W4M/Test_Matlab_data.txt", sep="\t", dec=",", header=TRUE, row.names=1)
design <- read.table("E:/PROJETS/Asca_W4M/Test_Matlab_design.txt", sep="\t", dec=",", header=TRUE)
design[,1] <- as.factor(design[,1])
design[,2] <- as.factor(design[,2])

## Verifier noms


fit <- lmdme(model=~F1 + F2 + F1:F2, data=data, design=design)

permuted <- permutation(model=~F1*F2, data=data, design=design, NPermutations=100, nCpus=3)

decomposition(fit, decomposition = "pca", scale="none", type="coefficient")




ssq <- function(fit)
{
  Overall_means <- sum(sum(fitted.values(fit)$'(Intercept)'^2))/sum(sum(data^2))
  Factors <- c(sum(sum(fitted.values(fit)$'F1'^2))/sum(sum(data^2)), sum(sum(fitted.values(fit)$'F2'^2))/sum(sum(data^2)))
  Interactions <- sum(sum(fitted.values(fit)$'F1:F2'^2))/sum(sum(data^2))
  Residuals <- 1 - Overall_means - Factors[1] - Factors[2] - Interactions
  
  return(list(Overall_means, Factors, Interactions, Residuals))
}


par(mfrow=c(2,2))
biplot(fit, xlabs="o", mfcol=NULL)
##Just the term of interest
biplot(fit, xlabs="o", term="F1")
##In separate graphics
biplot(fit, xlabs="o", term=c("F1", "F2"), mfcol=c(1,1))
##All terms in the same graphic
biplot(fit, xlabs="o", mfcol=c(1,3))

test <-lapply(permuted, FUN = ssq)
test1 <- matrix(unlist(test), ncol=5, byrow=TRUE)
test[2:101]

apply(apply(test1, 2, FUN = function(x){x > x[1]})[-1,], 2 , sum) / (length(test)-1)

score_moyen <- data.frame(fit@components$F1$rotation)

score <- data.frame(cbind(design, t(fit@residuals$'F1:F2')%*%fit@components$F1$rotation))
  
pc <- fit@components$F1$sdev / sum(fit@components$F1$sdev)

sp <- ggplot(score_moyen, aes(x=PC1, y=PC2))
sp + geom_point(size=2) + xlab(paste("PC1", round(pc[1]*100,1), "%")) + ylab(paste("PC2", round(pc[2]*100,1), "%")) + 
  geom_text(data=score, aes(PC1, PC2, label=Ind, col=F1), size=4, hjust=0, nudge_x=0.05, vjust=0, nudge_y=0.5)


