##XLConnect SUPPORT###
filename<-(paste(rootpath,"/Bioinformatics Resources/Templates/template.xls", sep=""))
style.wb <- loadWorkbook (filename)
cellstyle <- getCellStyle (style.wb, "GeneData")
rm(style.wb)
