#.ClassDefinitions.R

setClass("LimmaObj", slots=c(
  Contrasts="data.frame",
  DEGenes="list",
  AllGenes="list",
  Inputs="list"), package="probedeeper")

setClass("ColObj", slots=c(
  assign="data.frame",
  classvec="factor",
  full="list",
  match="list"), package="probedeeper")


setClass("PDObj", slots=c(
  eset="ExpressionSet",
  ColObj="ColObj",
  LimmaObj="LimmaObj"), package="probedeeper")

setClass("MultipleClassGSEAObject", representation(
  ObjectInfo="list",
  #Contains GSEAcomplist, GMTList, timestamp, odirectory
  Data="list",
  #Contains df of Comparisons grouped by GMTList
  PlotData="list", Stats="list"))
#Contains statistical data for Enrichment plots

setClass("MA.PP.datalist", representation(
  Suffix="character",
  Phenodata="data.frame",
  Data="list",
  Other="list"))

setClass("DEAnalyzerInput", representation(
  InputGeneSet="list",
  DEData="list"))

setClass("DEAnalyzerCalcs", representation(
  ObjectInfo="list",
  InputGeneSet="list",
  DEData="list",
  FoldChange="list",
  Pvalues="list",
  UpDnClassvec="list",
  Up="list",
  Down="list",
  GO.BP.UP="list",
  GO.BP.DN="list",
  GO.MF.UP="list",
  GO.MF.DN="list",
  PW.UP="list",
  PW.DN="list"))





