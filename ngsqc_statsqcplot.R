#' Basic sequencing statistics and plots
#'
#' @description
#' Create genreal alignment and insert size statistics using *.stats files generated using samtools stats.
#'  
#' @param FileDir directory that contains *.stats files
#' @param SampleTable data.frame with 3 columns: (1) File names without extensions, (2) Sample names, (3-X) Sample labels, e.g. Tumor/Normal (multiple columns available)
#' @param FilePrefix prefix of the output file names
#' @param Version version of the analysis e.g. "1"
#' @param FileType output plots file type (PDF/PNG)
#' @param InsertSmooth insert size plot smoothing parameter
#' @param xlimit limit of the x axis on the insert size plot
#' @param X_Marker marker line on the x axis of the read count bar plot e.g. indicating expected number of reads
#' @param BarHeight height of the bar plot, to be adjusted based on the number of samples in the analysis
#'
#' @return 
#' This function creates a statistics table and two plots: 1) insert size distribution and 2) bar plots with obtained/filtered read counts.
#' 
#' @export
#'
#' @examples
#' statsqcplot('QC_dir',SampleTable,'File_prefix')
ngsqc_statsqcplot = function(FileDir,SampleTable,FilePrefix,Version=1,FileType="PNG",InsertSmooth=1,xlimit=NA,X_Marker=NA,BarHeight=2.5) {

  library('reshape2')
  library('ggplot2')
  library('RColorBrewer')

  files <- list.files(path = FileDir, pattern="stats$")

  if (length(files)==0) {
    print('WARNING: Missing *.stats files')
  } else {
    labs = gsub(".stats", "", files, perl=TRUE)
    statistics = inserts = data.frame()

    for (i in 1:length(files)) {
      dta = read.table(paste0(FileDir,'/',files[i]),sep="\t",comment.char = "#",col.names=1:50,fill=T,stringsAsFactors = F)

      SampleName = SampleTable[SampleTable[,1] == labs[i],2]

      if (length(SampleName)==0) {
        stop(paste0("Sample table doesn't contain information for ",labs[i]," file"))
      }

      #Basic statistics
      tStatistics = dta[dta$X1=="SN",2:3]
      tStatistics[,1] = gsub(':','',tStatistics[,1])
      colnames(tStatistics) = c('ID',SampleName)


      #Insert size distribution
      tInserts = dta[dta$X1=="IS",2:3]
      tInserts[,1] = as.character(tInserts[,1])
      tInserts[,2] = as.numeric(tInserts[,2])
      tInserts[,2] = tInserts[,2]/sum(tInserts[,2])
      colnames(tInserts) = c('Length','Fraction')
      tInserts$Length = as.numeric(tInserts$Length)
      tInserts = tInserts[tInserts$Length>0,]  ######### is it good????

      #add median insert size
      FractionCumsum = cumsum(tInserts$Fraction)
      insMedian = head(tInserts[FractionCumsum>0.5,],1)$Length
      taStatistics = data.frame('insert size median',insMedian)
      colnames(taStatistics) = colnames(tStatistics)
      tStatistics = rbind(tStatistics,taStatistics)


      maxID = max(as.numeric(tInserts$Length))
      Npt = ceiling(maxID/InsertSmooth)
      tInsertsSmooth=NULL
      if (InsertSmooth>1) {
        for(z in 1:Npt) {
          idx = ((z-1)*InsertSmooth+1) : (z*InsertSmooth)
          tInsertsSmooth = rbind(tInsertsSmooth,colMeans(tInserts[idx,]))
        }
        tInserts=data.frame(tInsertsSmooth)
      }
      tInserts$Sample = SampleName


      if (i==1) {
        statistics = tStatistics
        inserts = tInserts
      } else {
        statistics = merge(statistics,tStatistics,by='ID',sort = F)
        inserts = rbind(inserts,tInserts)
      }
    }

    Dinserts = merge(inserts,SampleTable[,1:3],by.x='Sample',by.y=colnames(SampleTable)[2])
    colnames(Dinserts)[5] = 'Group'
    Dinserts$Group = as.character(Dinserts$Group)
    Dinserts = Dinserts[!is.na(Dinserts$Length),]

    if(is.na(xlimit)) {
      XMAX=max(Dinserts$Length[Dinserts$Fraction>0.01*max(Dinserts$Fraction) & Dinserts$Length<max(Dinserts$Length)])
      XMAX = ceiling(XMAX/100)*100
    } else {
      XMAX = xlimit
    }

    YMAX = max(Dinserts$Fraction[Dinserts$Length<XMAX])
    nFactors = length(unique(Dinserts$Group))
    if (nFactors<=5) {
      colors = c("#4DAF4A","#377EB8","red","orange",'darkblue')
    } else if (nFactors<=9) {
      colors = brewer.pal(n = nFactors, name = "Set1")
    } else {
      colors = brewer.pal(nFactors,'Spectral')
    }
    p=ggplot(Dinserts)+geom_line(aes(x=Length,y=Fraction,group=Sample,color=Group),alpha=0.5)+ 
      theme_bw()+ 
      scale_color_manual(values=colors) +
      xlab('Insert size [bp]') + 
      ylab('Fraction of reads')+
      xlim(0,XMAX)+ylim(0,YMAX)

    if(FileType=="PNG") {
      ggsave(p,file=paste0(FilePrefix,"_InsertSize_v",Version,".png"), scale=2, width=4, height=2.5, dpi=600)
    } else if(FileType=="PDF") {
      ggsave(p,file=paste0(FilePrefix,"_InsertSize_v",Version,".pdf"), scale=2, width=4, height=2.5)
    }



    library(data.table)
    exp_statistics = transpose(statistics)
    colnames(exp_statistics)=exp_statistics[1,]
    exp_statistics = exp_statistics[-1,]
    exp_statistics = cbind(SampleID=colnames(statistics)[-1], exp_statistics)
    for (i in 2:dim(exp_statistics)[2]) {
      exp_statistics[,i] = as.numeric(exp_statistics[,i])
    }
    exp_statistics$`filtered sequences [%]` = exp_statistics$`filtered sequences`/exp_statistics$`raw total sequences` * 100
    exp_statistics$`coverage X` = exp_statistics$`bases mapped`/2937639396
    exp_statistics = exp_statistics[,c('SampleID','raw total sequences','filtered sequences','filtered sequences [%]','sequences','reads mapped','coverage X','reads mapped and paired','reads unmapped','reads properly paired','reads paired','reads MQ0','non-primary alignments','total length','bases mapped','bases mapped (cigar)','mismatches','error rate','average quality','insert size median','insert size average','insert size standard deviation','inward oriented pairs','outward oriented pairs','pairs with other orientation','pairs on different chromosomes','percentage of properly paired reads (%)')]


    #basic stats multi-sample plot
    ColNames = c('reads mapped','reads unmapped','filtered sequences')
    SimpColNames = c('mapped','unmapped','filtered')
    exp_statistics_sel = exp_statistics[,c('SampleID',ColNames)]
    colnames(exp_statistics_sel) = c('SampleID',SimpColNames)
    exp_statistics_md = reshape2::melt(exp_statistics_sel)
    colors = c('#A6CEE3','#1F78B4','#B2DF8A')

    exp_statistics_sel_prc = exp_statistics_sel
    exp_statistics_sel_prc[2:4] = exp_statistics_sel[2:4]/rowSums(exp_statistics_sel[2:4])*100
    exp_statistics_prc_md = reshape2::melt(exp_statistics_sel_prc)
    colnames(exp_statistics_prc_md)[3] = 'prc'

    exp_statistics_merge = merge(exp_statistics_md,exp_statistics_prc_md,by=c('SampleID','variable'))
    exp_statistics_merge$variable = factor(exp_statistics_merge$variable,rev(SimpColNames))
    exp_statistics_merge$prcLab = paste0(round(exp_statistics_merge$prc,1),'%')
    exp_statistics_merge$SampleID = factor(exp_statistics_merge$SampleID,levels=rev(SampleTable[,2]))



    p=ggplot(exp_statistics_merge,aes(x=SampleID,y=value,fill=variable))+geom_bar(stat="identity",position="stack")+theme_bw()+ ylab('Number of reads') + xlab('Sample')+
      coord_flip()+scale_fill_manual(values=colors) +geom_text(aes(label = prcLab), size = 3, position = position_stack(vjust = 0.5)) + scale_y_continuous(expand = c(0,0))
    if (!is.na(X_Marker)) {
      p=p+geom_hline(yintercept=X_Marker,linetype=2,color="red")
    }
    if(FileType=="PNG") {
      ggsave(p,file=paste0(FilePrefix,"_ReadStats_v",Version,".png"), scale=2, width=6, height=BarHeight, dpi=600)
    } else if(FileType=="PDF") {
      ggsave(p,file=paste0(FilePrefix,"_ReadStats_v",Version,".pdf"), scale=2, width=6, height=BarHeight)
    }


    library('openxlsx')
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "AlignmentStats", gridLines = T);   writeDataTable(wb = wb, sheet = 1, x = exp_statistics)
    saveWorkbook(wb,paste0(FilePrefix,"_AlignmentStats_v",Version,".xlsx"), overwrite = TRUE)

    Ncol = ncol(exp_statistics)
    colnames(exp_statistics)[2:Ncol] <- paste0("BasicStats.", colnames(exp_statistics)[2:Ncol])
    return(exp_statistics)
  }
}
