ASCA.PlotScoresPerLevel_w4m <- function (asca, ee, pcs="1,2", interaction=0, factorName="", factorModalite) 
{
  pcs <- as.numeric(strsplit(pcs, split=",")[[1]])
  y <- (asca[[ee]]$means.matrix + asca$remainder) %*% asca[[ee]]$svd$v
  t.list.x <- list()
  t.list.y <- list()
  list.color.pattern <- list()
  color.per.variable <- rep(0, dim(asca$data)[1])
  pattern.per.variable <- rep(0, dim(asca$data)[1])
  kColOptions <- c(24, 552, 254, 26, 84, 51, 652, 68, 76, 96, 
                   10, 60, 33, 245, 147, 12, 26, 164, 181, 52, 512, 344, 
                   201, 111)
  kPointOptions <- 1:30
  for (p in 1:dim(asca[[ee]]$level.combinations$row.pattern)[1]) {
    if (length(asca[[ee]]$level.combinations$row.pattern[p, 
                                                         ]) == 1) {
      list.color.pattern[[p]] <- c(kColOptions[p], kPointOptions[p])
    }
    else if (length(asca[[ee]]$level.combinations$row.pattern[p, 
                                                              ]) == 2) {
      list.color.pattern[[p]] <- c(kColOptions[asca[[ee]]$level.combinations$row.pattern[p, 
                                                                                         1]], kPointOptions[asca[[ee]]$level.combinations$row.pattern[p, 
                                                                                                                                                      2]])
    }
    else {
      list.color.pattern[[p]] <- c(kColOptions[asca[[ee]]$level.combinations$row.pattern[p, 
                                                                                         1]]%%9, floor(kPointOptions[asca[[ee]]$level.combinations$row.pattern[p, 
                                                                                                                                                               2]]/9))
    }
    color.per.variable[asca[[ee]]$level.combinations$indices.per.pattern[[p]]] <- list.color.pattern[[p]][1]
    pattern.per.variable[asca[[ee]]$level.combinations$indices.per.pattern[[p]]] <- list.color.pattern[[p]][2]
    t.list.x[[p]] <- y[asca[[ee]]$level.combinations$indices.per.pattern[[p]], 
                       pcs[1]]
    t.list.y[[p]] <- y[asca[[ee]]$level.combinations$indices.per.pattern[[p]], 
                       pcs[2]]
  }
  legend.colors.patterns <- do.call(rbind, list.color.pattern)
  if (interaction != 1){
      titre <- paste("PC", pcs[1], " vs PC", pcs[2], " - Factor ", factorName, sep="")
  }else {
    titre <- paste("PC", pcs[1], " vs PC", pcs[2], " - Interaction", sep="")
  }
  plot(asca[[ee]]$svd$t[, pcs[1]], asca[[ee]]$svd$t[, pcs[2]], 
       xlim=range(c(min(unlist(t.list.x)), max(unlist(t.list.x)))), 
       ylim=range(c(min(unlist(t.list.y)), max(unlist(t.list.y)))), 
       main=titre, 
       xlab=paste("PC", pcs[1], " (", formatC(asca[[ee]]$svd$var.explained[pcs[1]] * 100, digits=2, format="f"), "%)", sep=""), 
       ylab=paste("PC", pcs[2], " (", formatC(asca[[ee]]$svd$var.explained[pcs[2]] * 100, digits=2, format="f"), "%)", sep=""), 
       cex=1.5, lwd=3, col=colors()[color.per.variable], 
       pch=pattern.per.variable)
#  if (interaction != 1){
    legend(x="bottomright", legend=factorModalite, 
           cex=0.8, col=colors()[legend.colors.patterns[, 1]], pch=legend.colors.patterns[, 2])
#  }
#  else {
#    legend(x="bottomright", apply(asca[[ee]]$level.combinations$row.patterns, 1, paste, collapse=" "), 
#           cex=0.8, col=colors()[legend.colors.patterns[, 1]], pch=legend.colors.patterns[, 2])
#  }
    
  for (p in 1:length(t.list.x)) {
    points(t.list.x[[p]], t.list.y[[p]], col=colors()[list.color.pattern[[p]][1]], 
           pch=list.color.pattern[[p]][2])
  }
}
