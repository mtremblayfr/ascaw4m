ASCA.Calculate_w4m <- function (data, levels, equation.elements = "", scaling, only.means.matrix = FALSE, use.previous.asca = NULL) 
{
  ASCA.GetEquationElement <- function(asca, evaluation, previous.asca) {
    s <- list()
    s$factors.evaluated <- evaluation
    if (!is.null(previous.asca)) {
      s$level.combinations <- previous.asca[[paste(evaluation, 
                                                   collapse = "")]]$level.combinations
    }
    else {
      s$level.combinations <- ASCA.GetRowRepeats(asca$levels[, 
                                                             s$factors.evaluated, drop = FALSE])
    }
    s$means.matrix <- matrix(nrow = dim(asca$data)[1], ncol = dim(asca$data)[2])
    for (p in 1:dim(s$level.combinations$row.patterns)[1]) {
      mean.for.this.level.combination <- colMeans(asca$data[s$level.combinations$indices.per.pattern[[p]], 
                                                            , drop = FALSE])
      for (i in s$level.combinations$indices.per.pattern[[p]]) {
        s$means.matrix[i, ] <- mean.for.this.level.combination
      }
    }
    s
  }
  s <- list()
  dataAdjusted <- MetStaT.ScalePip(data, center = FALSE, scale = FALSE, 
                                   quietly = TRUE)
  s$ssq.mean <- sum(rep(dataAdjusted$center.vector/dataAdjusted$scale.vector, 
                        nrow(data))^2)
  s$ssq <- sum(data^2)
  s$data <- dataAdjusted$data
  if (!is.numeric(levels)) {
    stop("The supplied levels are not numeric.")
  }
  s$levels <- levels
  if (!only.means.matrix) {
    s$svd <- PCA.Calculate(s$data)
  }
  s$ee.names <- c()
  if (identical(equation.elements, "")) {
    equation.elements <- ASCA.GetPowerSet(c(1:dim(s$levels)[2]), 
                                          exclude.empty.set = TRUE)
  }
  if (is.character(equation.elements)) 
    equation.elements <- lapply(strsplit(strsplit(equation.elements, 
                                                  split = ",")[[1]], split = ""), as.numeric)
  for (ee in equation.elements) {
    for (f in ee) if (f > dim(levels)[2] || f < 1) {
      stop(paste("Factor ", f, " is beyond scope of study-design", 
                 sep = ""))
    }
  }
  if (dim(data)[1] != dim(levels)[1]) {
    stop(paste("Number of rows in data (", dim(data)[1], 
               ") and study design (", dim(levels)[1], ") do not match", 
               sep = ""))
  }
  order.to.evaluate.ee <- sort(as.numeric(unlist(lapply(equation.elements, 
                                                        paste, collapse = ""))), index.return = TRUE)$ix
  s$remainder <- s$data
  for (ee in order.to.evaluate.ee) {
    new.equation.element <- ASCA.GetEquationElement(s, equation.elements[[ee]], 
                                                    use.previous.asca)
    reductions <- ASCA.GetPowerSet(equation.elements[[ee]], 
                                   exclude.empty.set = TRUE, exclude.complete.set = TRUE)
    for (r in reductions) {
      new.equation.element$means.matrix <- new.equation.element$means.matrix - 
        s[[c(paste(r, collapse = ""))]]$means.matrix
    }
    new.equation.element$ssq <- sum(new.equation.element$means.matrix^2)
    if (!only.means.matrix) {
      s$remainder <- s$remainder - new.equation.element$means.matrix
      new.equation.element$reduced.matrix <- s$remainder
      new.equation.element$svd <- PCA.Calculate(new.equation.element$means.matrix)
    }
    ee.name <- paste(equation.elements[[ee]], collapse = "")
    s$ee.names <- c(s$ee.names, ee.name)
    s[[ee.name]] <- new.equation.element
  }
  s$ssq.remainder <- sum(s$remainder^2)
  if (!only.means.matrix) 
    asca.summary <- ASCA.GetSummary(s, quietly = TRUE)
  return(list(s, asca.summary))
}
 