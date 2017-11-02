PDViolin<-function(obj, gene, pc, main, intgroup, returnData=T, transform=T, replaced=F, normalized=F, ColObj=NULL){
  # gene<-"TUBA1B"
  # obj<-obj
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
                                                           (gene >= 1 & gene <= nrow(obj)))))
    if (!all(intgroup %in% names(colData(obj))))
      stop("all variables in 'intgroup' must be columns of colData")
    stopifnot(returnData | all(sapply(intgroup, function(v) is(colData(obj)[[v]],
                                                               "factor"))))
    if(is.null(ColObj)){
      ColObj<-ColObjCreate(colData(obj)[[intgroup]], pie = F)
    }
    if (missing(pc)) {
      pc <- if (transform)
        0.5
      else 0
    }
    if (is.null(sizeFactors(obj)) & is.null(normalizationFactors(obj))) {
      obj <- estimateSizeFactors(obj)
    }
    cnts <- counts(obj, normalized = normalized, replaced = replaced)[gene,
                                                                      ]
    group <- if (length(intgroup) == 1) {
      colData(obj)[[intgroup]]
    }
    else if (length(intgroup) == 2) {
      lvls <- as.vector(t(outer(levels(colData(obj)[[intgroup[1]]]),
                                levels(colData(obj)[[intgroup[2]]]), function(x,
                                                                              y) paste(x, y, sep = " : "))))
      droplevels(factor(apply(as.data.frame(colData(obj)[,
                                                         intgroup, drop = FALSE]), 1, paste, collapse = " : "),
                        levels = lvls))
    }
    else {
      factor(apply(as.data.frame(colData(obj)[, intgroup, drop = FALSE]),
                   1, paste, collapse = " : "))
    }
    data <- data.frame(count = cnts + pc, group = as.integer(group))
    logxy <- if (transform)
      "y"
    else ""
    if (missing(main)) {
      main <- if (is.numeric(gene)) {
        rownames(obj)[gene]
      }
      else {
        gene
      }
    }
    dat<-data.frame(count = data$count, colData(obj)[intgroup])
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
      return(dat)
  }
  
  if(obj.class=="DESeqTransform"){
    stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) &
                                                           (gene >= 1 & gene <= nrow(obj)))))
    if (!all(intgroup %in% names(colData(obj))))
      stop("all variables in 'intgroup' must be columns of colData")
    stopifnot(returnData | all(sapply(intgroup, function(v) is(colData(obj)[[v]],
                                                               "factor"))))
    if(is.null(ColObj)){
      ColObj<-ColObjCreate(colData(obj)[[intgroup]], pie = F)
    }
    if (missing(pc)) {
      pc <- if (transform)
        0.5
      else 0
    }
    cnts <-assay(obj)[gene,
                                                                      ]
    group <- if (length(intgroup) == 1) {
      colData(obj)[[intgroup]]
    }
    else if (length(intgroup) == 2) {
      lvls <- as.vector(t(outer(levels(colData(obj)[[intgroup[1]]]),
                                levels(colData(obj)[[intgroup[2]]]), function(x,
                                                                              y) paste(x, y, sep = " : "))))
      droplevels(factor(apply(as.data.frame(colData(obj)[,
                                                         intgroup, drop = FALSE]), 1, paste, collapse = " : "),
                        levels = lvls))
    }
    else {
      factor(apply(as.data.frame(colData(obj)[, intgroup, drop = FALSE]),
                   1, paste, collapse = " : "))
    }
    data <- data.frame(count = cnts + pc, group = as.integer(group))
    logxy <- if (transform)
      "y"
    else ""
    if (missing(main)) {
      main <- if (is.numeric(gene)) {
        rownames(obj)[gene]
      }
      else {
        gene
      }
    }
    dat<-data.frame(count = data$count, colData(obj)[intgroup])
    ylab <- "normalized expression"
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
      return(dat)
  }
}
