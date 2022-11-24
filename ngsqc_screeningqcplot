#' Create screening plots for result files generated using fastq_screen
#'
#' @param FileDir directory which contains fastq_screen result files (*.txt format)
#' @param DatasetName name of the dataset, used as a prefix in the file name
#' @param ScreenName name of the screen appended to the output plot name
#' @param Version version of the plot, appended to the output file name
#' @param TrimLabelPart string which will be removed from file names used to label features on the plot
#'
#' @return Contamination screen plot
#' @export
#'
#' @examples
#' screeningqcplot('QC',SampleLabels,'ContaminationScreen','IGCZ_T-AAL_miRge1')
ngsqc_screeningqcplot = function(FileDir,DatasetName,ScreenName,Version=1,TrimLabelPart=NULL) {

  library('tidyverse')
  library('ggplot2')

  files <- list.files(path = FileDir, pattern="txt$")

  if (length(files)>0) {
    labs = gsub("_screen.txt", "", files, perl=TRUE)
    if(!is.null(TrimLabelPart)) {
      labs = gsub(TrimLabelPart, "", labs, perl=TRUE)
    }
    screen = data.frame()
    baseMax = 0
    for (i in 1:length(files)) {
      dta = read.table(paste0(FileDir,'/',files[i]),sep="\t",comment.char = "",col.names=1:20,fill=T,stringsAsFactors = F)
      #dta = read_tsv(paste0(FileDir,'/',files[i]),skip = 1)

      colnames(dta) = dta[2,]
      tScreen = dta[-1:-2,!is.na(dta[2,])] %>%
        tibble() %>%
        filter(`#Reads_processed`!='') %>%
        mutate(GenomeSpecific = as.numeric(`%One_hit_one_genome`) + as.numeric(`%Multiple_hits_one_genome`),
               MultipleGenomes = as.numeric(`%One_hit_multiple_genomes`) + as.numeric(`%Multiple_hits_multiple_genomes`)) %>%
        select(Genome,GenomeSpecific,MultipleGenomes) %>%
        pivot_longer(!Genome) %>%
        mutate(Genome = gsub('_',':',Genome),
               Sample = labs[i])
      tScreen

      if (i==1) {
        screen = tScreen
      } else {
        screen = rbind(screen,tScreen)
      }
    }

    screen$Sample = factor(screen$Sample,levels = rev(labs))

    colors = c("#377EB8","#4DAF4A")
    p=ggplot(data=screen,aes(x=Sample,y=value,fill=name)) + geom_bar(stat="identity",position="stack") +
      theme_bw() +
      facet_wrap(~Genome,nrow=1) +
      coord_flip() +
      xlab('Sample') +
      ylab('Fraction of reads') +
      scale_fill_manual(values=colors,name="") +
      theme(legend.position="bottom")
    p
    ggsave(p,file=paste0(DatasetName,"_",ScreenName,"_v",Version,".png"), scale=2, width=6.7, height=3, dpi=600)

  } else {
    stop('Missing *.txt files')
  }
}
