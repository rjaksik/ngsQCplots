#' Create STAR alignment statistics plots
#'
#' @param FileDir directory which contains STAR log files (*.Log.final.out)
#' @param DatasetName name of the dataset, used as a prefix in the file name
#' @param Version version of the plot, appended to the output file name
#' @param TrimLabelPart string which will be removed from file names used to label features on the plot
#'
#' @return STAR alignment statistics
#' @export
#'
#' @examples
#' ngsqc_star_alignment('QC',SampleLabels,'IGCZ_T-AAL_miRge1')

ngsqc_star_alignment = function(FileDir,DatasetName,Version=1,TrimLabelPart=NULL) {

  library('tidyverse')
  library('ggplot2')

  total = 'Number of input reads'
  mapped = 'Uniquely mapped reads number'
  multi_map = c('Number of reads mapped to multiple loci','Number of reads mapped to too many loci')
  unmapped = c('Number of reads unmapped: too short','Number of reads unmapped: other','Number of reads unmapped: too many mismatches')

  files <- list.files(path = FileDir, pattern=".Log.final.out$")

  if (length(files)>0) {
    labs = gsub(".Log.final.out", "", files, perl=TRUE)
    if(!is.null(TrimLabelPart)) {
      labs = gsub(TrimLabelPart, "", labs, perl=TRUE)
    }
    alnstatsTable = data.frame()
    for (i in 1:length(files)) {

      #read all count files
      talnstatsTable =read.table(paste0(FileDir,'/',files[i]),header=F,comment.char = "",stringsAsFactors = F,sep="|",col.names=1:2,fill=T,strip.white =T)
      talnstatsTable=talnstatsTable[talnstatsTable$X1 %in% c(mapped,multi_map,unmapped),]
      colnames(talnstatsTable) = c('name','value')
      talnstatsTable$value = as.numeric(talnstatsTable$value)
      rownames(talnstatsTable) = talnstatsTable$name

      alnstatsTable  = rbind(alnstatsTable, data.frame('Sample' = labs[i],
                                                       'Uniquely mapped' = talnstatsTable[mapped,'value'],
                                                       'Multi mapping' = sum(talnstatsTable[multi_map,'value']),
                                                       'Unmapped' = sum(talnstatsTable[unmapped,'value']),check.names = FALSE ))
    }

    alnstatsTable_proc = alnstatsTable  %>%
      tibble() %>%
      pivot_longer(!Sample) %>%
      mutate(Sample = factor(Sample,levels = rev(labs)),
             name = factor(name,levels = rev(c('Uniquely mapped','Multi mapping','Unmapped'))))


    colors = c('#B1084C',"#7CB5EC","#377EB8")
    p=ggplot(data=alnstatsTable_proc,aes(x=Sample,y=value,fill=name)) + geom_bar(stat="identity",position="stack") +
      theme_bw() +
      coord_flip() +
      xlab('Sample') +
      ylab('Number of read pairs') +
      scale_fill_manual(values=colors,name="") +
      theme(legend.position="bottom")
    p
    ggsave(p,file=paste0(DatasetName,"_STAR_AlignmentStats_v",Version,".png"), scale=2, width=4, height=3, dpi=600)

  } else {
    stop('Missing *.Log.final.out files')
  }
}

