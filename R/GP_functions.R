###rJAVA###

getJobFiles<-function(job.result){
  if (is.null(job.result)) {
    stop("job.result is NULL")
  }
  download.directory<-as.character(job.result.get.job.number(job.result))
  if (is.null(download.directory)) {
    stop("download directory cannot be NULL")
  }
  files.array <- .jcall(job.result, "[Ljava/io/File;", "downloadFiles", 
                        download.directory, overwrite=TRUE, check = FALSE)
  if (!is.null(e <- .jgetEx(clear = TRUE))) {
    e.message <- .jcall(e, "S", "toString")
    stop(e.message)
  }
  filenames.list <- list()
  for (i in 1:length(files.array)) {
    file <- files.array[[i]]
    filenames.list[i] <- .jcall(file, "S", "getPath")
  }
  workdir<-getwd()
  filenames.location<-lapply(filenames.list, function(X) paste(workdir, X, sep="/"))
  return(filenames.location)
}
