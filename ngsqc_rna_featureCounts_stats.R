#' Plot read alignment and assignment statistics based on featureCounts result file
#'
#' @param CountsFile featureCounts main counts file, data columns should have the .bam suffix
#' @param SummaryFile featureCounts summary file (.summary suffix)
#' @param DatasetName name of the dataset, used as a prefix in the file name
#' @param geneTypeCol column name that contains the gene classes, if left NULL the script will look for either type_of_gene or gene_type column
#' @param Version version of the plot, appended to the output file name
#'
#' @return plot which shows the percentage of reads that are associated with certain classes of RNAs
#' Warning:  this script ignores the multiassignment problem (one read assigned to multiple genes),
#' therefore the total number of reads obtained for all features may not match the total number of assigned reads
#' @export
#'
#' @examples
ngsqc_rna_featureCounts_stats = function(CountsFile, SummaryFile, DatasetName, geneTypeCol = NULL, Version=1) {

  library('openxlsx')

  #read featureCounts read counts
  fcCounts = read.table(CountsFile,header=T,sep="\t",check.names = F)
  tSampleCol = colnames(fcCounts)
  colnames(fcCounts) = gsub('.bam','',colnames(fcCounts))
  sampleCol = colnames(fcCounts)[grepl('.bam',tSampleCol)]

  #read featureCounts statistics
  fcStats = read.table(SummaryFile,header=T,sep="\t",stringsAsFactors = F,check.names = F)
  rownames(fcStats) = fcStats[,1]
  fcStats[,1] = NULL
  samples = gsub('.bam','',colnames(fcStats))
  fcStatsSel = reshape2::melt(data.frame(samples=samples,#TotalReads=colSums(fcStats),
                                         Unmapped = as.numeric(fcStats['Unassigned_Unmapped',]),
                                         Unassigned = as.numeric(fcStats['Unassigned_NoFeatures',]+fcStats['Unassigned_Ambiguity',])))
  #basic checks
  if(is.null(geneTypeCol)) {
    candidateGeneTypeCol = c('gene_type','type_of_gene')
    geneTypeCol = candidateGeneTypeCol[candidateGeneTypeCol %in% colnames(fcCounts)]
  }
  if (!geneTypeCol %in% colnames(fcCounts)) {
    stop(paste0('Counts file doesnt contain gene_type/type_of_gene column please provide valid  geneTypeCol parameter'))
  }

  colnames(fcStatsSel) = c('variable',geneTypeCol,'sum')
  geneIDCol = colnames(fcCounts)[1]
  GeneDataAnnotWide = reshape2::melt(fcCounts[,c(geneIDCol,geneTypeCol,sampleCol)],id.vars=c(geneIDCol,geneTypeCol))

  dfwc= summarySE(data=GeneDataAnnotWide, measurevar="value", groupvars=c("variable",geneTypeCol), na.rm=T, conf.interval=.95)
  dfwc = dfwc[,c('variable',geneTypeCol,'sum')]
  PlotData = rbind(dfwc,fcStatsSel)

  SelectedTypes = c('Unmapped','Unassigned','rRNA','snRNA','lncRNA','pri-miRNA','protein_coding')
  PlotData$type_of_gene = PlotData[,geneTypeCol]
  PlotData$type_of_gene[!PlotData$type_of_gene %in% SelectedTypes] = 'other'

  PlotData$type_of_gene = factor(PlotData$type_of_gene,
                                 levels=c('Unmapped','Unassigned','rRNA','snRNA','lncRNA','pri-miRNA','other','protein_coding'))

  colors = c("darkgray",'gray',"#FFFF33","#9E7FA5","#FF7F00","#E41A1C","#4DAF4A","#377EB8")
  p=ggplot() +geom_bar(data=PlotData,aes(x=variable,y=sum,fill=type_of_gene),stat = "identity",position="fill")  +
    scale_fill_manual(values=colors) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0)) +
    theme(legend.position="bottom") +
    scale_y_continuous(expand = c(0, 0)) +
    xlab('Sample') + ylab('Fraction of reads') +
    coord_flip()+
    scale_x_discrete(limits = rev(levels(PlotData$variable)))
  p
  ggsave(p,file=paste0(DatasetName,'_FeatureTypeReads_v',Version,'.png'), scale=2, width=4, height=3)

  #wb <- createWorkbook()
  #addWorksheet(wb = wb, sheetName = "ReadNumbers", gridLines = T);   writeDataTable(wb = wb, sheet = 1, x = dcast(PlotData,variable~type_of_gene))
  #saveWorkbook(wb, 'Table S1 - FeatureTypeReads_v1.xlsx', overwrite = TRUE)
}
