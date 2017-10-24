PDViolin<-function(obj, gene, pc, main, intgroup, returnData=T, transform=T){
  # gene<-"TUBA1B"
  # obj<-dds
  # intgroup<-"group"
  # returnData = FALSE
  # transform = TRUE
  # pc<-0.5
  # replaced =FALSE
  # normalized=TRUE
  # main<-gene


  obj.class<-class(obj)
  if(obj.class=="DESeqDataSet"){
    stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) &
                                                           (gene >= 1 & gene <= nrow(dds)))))
    if (!all(intgroup %in% names(colData(dds))))
      stop("all variables in 'intgroup' must be columns of colData")
    stopifnot(returnData | all(sapply(intgroup, function(v) is(colData(dds)[[v]],
                                                               "factor"))))
    if (missing(pc)) {
      pc <- if (transform)
        0.5
      else 0
    }
    if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
      dds <- estimateSizeFactors(dds)
    }
    cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene,
                                                                      ]
    group <- if (length(intgroup) == 1) {
      colData(dds)[[intgroup]]
    }
    else if (length(intgroup) == 2) {
      lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]),
                                levels(colData(dds)[[intgroup[2]]]), function(x,
                                                                              y) paste(x, y, sep = " : "))))
      droplevels(factor(apply(as.data.frame(colData(dds)[,
                                                         intgroup, drop = FALSE]), 1, paste, collapse = " : "),
                        levels = lvls))
    }
    else {
      factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]),
                   1, paste, collapse = " : "))
    }
    data <- data.frame(count = cnts + pc, group = as.integer(group))
    logxy <- if (transform)
      "y"
    else ""
    if (missing(main)) {
      main <- if (is.numeric(gene)) {
        rownames(dds)[gene]
      }
      else {
        gene
      }
    }
    ylab <- ifelse(normalized, "normalized count", "count")
    g<-ggplot(dat)+
      geom_violin(aes(x=group, y=count, fill=group), scale="width", trim=F)+
      scale_fill_manual(values = ColObj@match$line)+
      geom_boxplot(aes(x=group, y=count, fill=group), width=0.5, fill=ColObj@match$fill)+
      geom_point(aes(x=group, y=count), color="black", size=2)+
      geom_point(aes(x=group, y=count), color=ColObj@full$fill, size=1)+
      scale_shape(solid = FALSE)+
      ggtitle(label=main)+
      theme_bw()
    print(g)
    if (returnData)
      return(data.frame(count = data$count, colData(dds)[intgroup]))
  }
  if(obj.class=="DESeqTransform"){
    stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) &
                                                           (gene >= 1 & gene <= nrow(dds)))))
    if (!all(intgroup %in% names(colData(dds))))
      stop("all variables in 'intgroup' must be columns of colData")
    stopifnot(returnData | all(sapply(intgroup, function(v) is(colData(dds)[[v]],
                                                               "factor"))))
    if (missing(pc)) {
      pc <- if (transform)
        0.5
      else 0
    }
    if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
      dds <- estimateSizeFactors(dds)
    }
    cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene,
                                                                      ]
    group <- if (length(intgroup) == 1) {
      colData(dds)[[intgroup]]
    }
    else if (length(intgroup) == 2) {
      lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]),
                                levels(colData(dds)[[intgroup[2]]]), function(x,
                                                                              y) paste(x, y, sep = " : "))))
      droplevels(factor(apply(as.data.frame(colData(dds)[,
                                                         intgroup, drop = FALSE]), 1, paste, collapse = " : "),
                        levels = lvls))
    }
    else {
      factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]),
                   1, paste, collapse = " : "))
    }
    data <- data.frame(count = cnts + pc, group = as.integer(group))
    logxy <- if (transform)
      "y"
    else ""
    if (missing(main)) {
      main <- if (is.numeric(gene)) {
        rownames(dds)[gene]
      }
      else {
        gene
      }
    }
    ylab <- ifelse(normalized, "normalized count", "count")
    g<-ggplot(dat)+
      geom_violin(aes(x=group, y=count, fill=group), scale="width", trim=F)+
      scale_fill_manual(values = ColObj@match$line)+
      geom_boxplot(aes(x=group, y=count, fill=group), width=0.5, fill=ColObj@match$fill)+
      geom_point(aes(x=group, y=count), color="black", size=2)+
      geom_point(aes(x=group, y=count), color=ColObj@full$fill, size=1)+
      scale_shape(solid = FALSE)+
      ggtitle(label=main)+
      theme_bw()
    print(g)
    if (returnData)
      return(data.frame(count = data$count, colData(dds)[intgroup]))
  }
}
