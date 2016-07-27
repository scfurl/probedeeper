

GSEA.write.gct <- function (gct, filename) 
{
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct)[1], "\t", dim(gct)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "\t", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  names <- names(gct)
  cat("\t", names[1], file = f, append = TRUE, sep = "")
  for (j in 2:length(names)) {
    cat("\t", names[j], file = f, append = TRUE, sep = "")
  }
  cat("\n", file = f, append = TRUE, sep = "\t")
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] +  2)
  m[, 1] <- row.names(gct)
  m[, 2] <- row.names(gct)
  index <- 3
  for (i in 1:dim(gct)[2]) {
    m[, index] <- gct[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
  return(gct)
}

GSEA.ConsPlot <- function(V, col.names, main = " ", sub = " ", xlab=" ", ylab=" ") {
  
  # Plots a heatmap plot of a consensus matrix
  
  cols <- length(V[1,])
  B <- matrix(0, nrow=cols, ncol=cols)
  max.val <- max(V)
  min.val <- min(V)
  for (i in 1:cols) {
    for (j in 1:cols) {
      k <- cols - i + 1
      B[k, j] <-  max.val - V[i, j] + min.val
    }
  }
  
  
  
  #     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
  col.map <- rev(c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D"))
  
  #     max.size <- max(nchar(col.names))
  par(mar = c(5, 15, 15, 5))
  image(1:cols, 1:cols, t(B), col = col.map, axes=FALSE, main=main, sub=sub, xlab= xlab, ylab=ylab)
  
  for (i in 1:cols) {
    col.names[i]  <- substr(col.names[i], 1, 25)
  }
  col.names2 <- rev(col.names)
  
  size.col.char <- ifelse(cols < 15, 1, sqrt(15/cols))
  
  axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.col.char, font.axis=1, line=-1)
  axis(3, at=1:cols, labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=1, line=-1)
  
  return()
}




GSEA.HeatMapPlot2 <- function(V, row.names = "NA", col.names = "NA", main = " ", sub = " ", xlab=" ", ylab=" ", color.map = "default") {
  #
  # Plots a heatmap of a matrix
  
  n.rows <- length(V[,1])
  n.cols <- length(V[1,])
  
  if (color.map == "default") {
    color.map <- rev(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5))
  }
  
  heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
  heatm[1:n.rows,] <- V[seq(n.rows, 1, -1),]
  
  par(mar = c(7, 15, 5, 5))
  image(1:n.cols, 1:n.rows, t(heatm), col=color.map, axes=FALSE, main=main, sub = sub, xlab= xlab, ylab=ylab)
  
  if (length(row.names) > 1) {
    size.row.char <- ifelse(n.rows < 15, 1, sqrt(15/n.rows))
    size.col.char <- ifelse(n.cols < 15, 1, sqrt(10/n.cols))
    #            size.col.char <- ifelse(n.cols < 2.5, 1, sqrt(2.5/n.cols))
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 40)
    }
    row.names <- row.names[seq(n.rows, 1, -1)]
    axis(2, at=1:n.rows, labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=1, line=-1)
  }
  
  if (length(col.names) > 1) {
    axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
  }
  
  return()
}