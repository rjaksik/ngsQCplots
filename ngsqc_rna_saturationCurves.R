#' Plot RNA-seq feature number saturation curves
#'
#' @description
#' This plot shows the total number of transcripts with a minimum of 10x coverage independently of the total number of reads (by sub sampling the counts table),
#' and the total number of reads obtained for each sample, which marks the end of each curve.
#'
#' @param CountsFile featureCounts main counts file, data columns should have the .bam suffix, all other annotation columns are unused
#' @param DatasetName name of the dataset, used as a prefix in the file name
#' @param Version version of the plot, appended to the output file name
#'
#' @return saturation plot
#' @export
#'
#' @examples
ngsqc_rna_saturationCurves = function(CountsFile, DatasetName, Version=1) {

  library(ggplot2)
  library(ggrepel)

  fcCounts = read.table(CountsFile,header=T,sep="\t",check.names = F)
  tSampleCol = colnames(fcCounts)
  colnames(fcCounts) = gsub('.bam','',colnames(fcCounts))
  sampleCol = colnames(fcCounts)[grepl('.bam',tSampleCol)]
  fcCounts_dta = fcCounts[,sampleCol]

  MaxReads = colSums(fcCounts_dta)
  MaxInterval = ceiling(max(MaxReads)/100000000)*100000000
  Intervals = seq(0,MaxInterval,by=1e07)
  Nsamp=length(sampleCol)

  PlotData = data.frame()
  for(k in 1:Nsamp) {
    sample = colnames(fcCounts_dta)[k]
    for (i in 2:(length(Intervals)-1)) {
      scale = Intervals[i]/MaxReads
      if (scale[k]>1) {
        scale[k] = 1
        xpt = MaxReads[k]
        lab=sample
      } else {
        xpt = Intervals[i]
        lab=NA
      }
      nRNA = sum(round(fcCounts_dta[,k]*scale[k]) >= 10)
      PlotData = rbind(PlotData,data.frame(sample,xpt,nRNA,lab))
      if (!is.na(lab)) break #this is the last point
    }
  }

  PlotData_labels = PlotData[!is.na(PlotData$lab),]

  colors = c("#4DAF4A","#377EB8","red","orange")
  customtheme = theme(panel.grid.minor = element_blank())+ theme(legend.position="none") + theme_bw() ;
  p <- ggplot(PlotData, aes(x=xpt/1000/1000,y= nRNA, group=sample)) +
    geom_line(color="#377EB8")+
    customtheme +
    geom_point()+
    scale_color_manual(values=colors) +
    xlab('Reads [mln]') + ylab('Transcripts at 10x') +
    geom_label_repel(data=PlotData_labels, aes(x=xpt/1000/1000, y= nRNA, group=sample,label=lab), min.segment.length = 0, size=2, color="orange",alpha=0.8)
  p
  ggsave(p,file=paste0(DatasetName,'_SaturationCurves_v',Version,'.png'), scale=2, width=5, height=3.5)
}
