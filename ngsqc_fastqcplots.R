#' Visualize QC statistics for all samples based on FastQC results processed with  fastqc_process.py custom script
#' example: fastqcplots('QC',SampleTable,'IGCZ_T-AAL_miRge1')
#' NOTE: NovaSeq qualities are based only one values: 2, 12, 23 and 37
#'
#' @param FileDir 
#' @param SampleTable data.frame with 3 columns: (1) File names without extensions, (2) Sample names, (3-X) Sample labels, e.g. Tumor/Normal (multiple columns available)
#' @param FilePrefix 
#' @param HeatmapHeight 
#' @param AlphaCol 
#' @param Version 
#'
#' @return
#' @export
ngsqc_fastqcplots = function(FileDir, SampleTable, FilePrefix, HeatmapHeight=8, AlphaCol = 0.8, Version=1) {
  
  library('scales')
  library('RColorBrewer')
  library('pheatmap')
  library('grid')
  library('openxlsx')
  library('reshape2')
  library('ggplot2')
  
  #convert range in string format e.g. 1-3 to range middle
  StrRangeToNum = function(sCols) {
    resCol = rep(0,length(sCols))
    for (i in 1:length(sCols)) {
      if (is.numeric(sCols[i])) {
        tmpstr=sCols[i]
      } else {
        tmpstr = unlist(strsplit(sCols[i],'-'))
      }
      if (length(tmpstr)>1) {
        resCol[i] = (as.numeric(tmpstr[1]) + as.numeric(tmpstr[2]))/2
      } else {
        resCol[i] = sCols[i]
      }
    }
    return(resCol)
  }
  
  #Create heatmap with stats
  drawQCheatmap = function(stats,xlabel,plotName,fileName,plotwidth,FactorsTable,SampNames) {
    ColorMatrix = matrix(0,nrow=length(stats),ncol=length(stats[[1]][,1]))
    colnames(ColorMatrix) = stats[[1]][,1]
    rownames(ColorMatrix) = SampNames
    
    for (i in 1:length(stats)) {
      tmp = stats[[i]][,2]
      names(tmp) = as.character(stats[[i]][,1])
      ColorMatrix[i,] = tmp[colnames(ColorMatrix)]
    }
    ColorMatrix[is.na(ColorMatrix)] = 0
    
    png(fileName,width = plotwidth, height = HeatmapHeight, units = "in", res=600)
    setHook("grid.newpage", function() pushViewport(viewport(x=1,y=0.95,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
    print(pheatmap(ColorMatrix,cluster_col=F,annotation_row=FactorsTable))
    setHook("grid.newpage", NULL, "replace")
    grid.text(xlabel, y=-0.02, gp=gpar(fontsize=16))
    grid.text("Sample", x=-0.07, rot=90, gp=gpar(fontsize=16))
    grid.text(plotName, y=1.02, gp=gpar(fontsize=16))
    dev.off()
  }
  
  
  readFlagstatStats = function(FileName) {
    flagstat <- read.table(FileName,stringsAsFactors=F,sep="#",header=F,nrow=12)
    stats = data.frame(stringsAsFactors = F)
    for (i in 1:dim(flagstat)[1]) {
      tStr = unlist(strsplit(flagstat[i,1],'+',fixed=T))
      Name =  unlist(strsplit(tStr[2],'(',fixed=T))[1]
      Name = trimws(sub(" 0 ","",Name))
      Name = sub("in total","total QC-passed reads",Name)
      stats = rbind(stats,data.frame(Name=Name,Value = as.numeric(trimws(tStr[1]))))
    }
    colnames(stats)[2] = basename(FileName)
    return(stats)
  }
  
  
  
  # Generate file names -> sample names translator
  SampNames = SampleTable[,2]
  names(SampNames) = SampleTable[,1]
  Nfiles = length(SampNames)
  
  #extract factor table
  FactorsTable = data.frame(SampleTable[,3:dim(SampleTable)[2]])
  rownames(FactorsTable) = SampleTable[,2]
  colnames(FactorsTable) = colnames(SampleTable)[3:dim(SampleTable)[2]]
  
  # Generate color vector
  factors = as.character(SampleTable[,3])
  Nfactor = length(unique(factors))
  if (Nfactor==1) {
    cols='#377EB8'
  } else if (Nfactor<=9)  {
    cols = brewer.pal(Nfactor,'Set1')
  } else {
    cols = brewer.pal(Nfactor,'Spectral')
  }
  labsDF = data.frame(row.names=SampNames, factor=factors, cols=cols[as.factor(factors)],stringsAsFactors = F)
  ulabsDF = unique(labsDF)
  
  ################ Alignment statistics ################
  
  files <- list.files(path = FileDir, pattern="flagstat$")
  
  if (length(files)>0) {
    labs = gsub(".flagstat", "", files, perl=TRUE)
    flagstats <- data.frame()
    for (i in 1:length(files)) {
      tStats <- readFlagstatStats(paste0(FileDir,'/',files[i]))
      if (i==1) {
        flagstats = tStats
      } else {
        flagstats = merge(flagstats,tStats,sort = F)
      }
    }
    rownames(flagstats) = flagstats$Name
    flagstats$Name=NULL
    flagstats = data.frame(t(flagstats),check.names =F)
    flagstats$`duplicates [%]` = flagstats$`duplicates`/flagstats$`total QC-passed reads` *100
    flagstats$`mapped [%]` = flagstats$`mapped`/flagstats$`total QC-passed reads` *100
    rownames(flagstats) = labs
    flagstats = cbind(data.frame(SampNames[rownames(flagstats)]),flagstats)
    colnames(flagstats)[1] = 'SampleID'
    
    
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "AlignmentStats", gridLines = T);   writeDataTable(wb = wb, sheet = 1, x = flagstats)
    saveWorkbook(wb,paste0(FilePrefix,'_AlignmentStats_v',Version,'.xlsx'), overwrite = TRUE)
    
    
    PlotData = flagstats[,c('SampleID','total QC-passed reads','duplicates')]
    PlotData$`unique reads` = PlotData$`total QC-passed reads` - PlotData$duplicates
    PlotData$`total QC-passed reads`=NULL
    PlotData = reshape2::melt(PlotData,id.vars="SampleID")
    p=ggplot(PlotData, aes(x = SampleID, y = value, fill = variable)) + geom_bar(stat = "identity") + theme_bw()  + scale_fill_brewer(palette="Set1") +xlab('Sample') +ylab('Number of mapped reads') +coord_flip()
    ggsave(p,file=paste0(FilePrefix,'_AlignmentStats_v',Version,'.png'), scale=2, width=4, height=0.2*length(files)+0.2, dpi=600)
  } else {
    print('WARNING: Missing flagstat files')
  }
  
  ################ Exon Coverage ################
  
  #from: http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
  
  files <- list.files(path = FileDir,pattern="tcov$")
  
  if (length(files)>0) {
    fnames = gsub("_ExomeCap.tcov", "", files, perl=TRUE)
    labs = as.character(SampNames[fnames])
    
    # Create lists to hold coverage and cumulative coverage for each alignment,
    # and read the data into these lists.
    cov <- list()
    cov_cumul <- list()
    baseMax = 0
    covMax = 0
    for (i in 1:length(files)) {
      cov[[i]] <- read.table(paste0(FileDir,"/",files[i]))
      cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
      tmax = max(cov[[i]][cov_cumul[[i]]>0.01,2])
      if (tmax>baseMax) {
        baseMax = tmax
      }
      tmax = max(cov_cumul[[i]])
      if (tmax>covMax) {
        covMax = tmax
      }
    }
    
    CovWidth = ceiling(baseMax/50)*50
    if (CovWidth>1000) CovWidth=1000
    CovHeight = ceiling(covMax/0.1)*0.1
    
    # Save the graph to a file
    png(paste0(FilePrefix,"_ExonCoverage_v",Version,".png"),width = 6, height = 6, units = "in", res=600)
    
    # Create plot area, but do not plot anything. Add gridlines and axis labels.
    plot(cov[[1]][2:(CovWidth+1), 2], cov_cumul[[1]][1:CovWidth], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,CovHeight), main="Exon Coverage")
    if (CovWidth<500) {
      abline(v = 20, col = "gray60")
      abline(v = 50, col = "gray60")
      abline(v = 80, col = "gray60")
      axis(1, at=c(20,50,80), labels=c(20,50,80))
    }
    abline(v = 100, col = "gray60")
    abline(h = 0, col = "gray60")
    abline(h = 0.50, col = "gray60")
    abline(h = 0.90, col = "gray60")
    axis(2, at=c(0.90), labels=c(0.90))
    axis(2, at=c(0.50), labels=c(0.50))
    axis(1, at=c(100), labels=c(100))
    
    # Actually plot the data for each of the alignments (stored in the lists).
    for (i in 1:length(cov)) points(cov[[i]][2:(CovWidth+1), 2], cov_cumul[[i]][1:CovWidth], type='l', lwd=2, col=alpha(labsDF[labs[i],'cols'], AlphaCol))
    if(dim(ulabsDF)[1]>1) legend("topright", legend=ulabsDF$factor, col=ulabsDF$cols, lty=1, lwd=4)
    
    dev.off()
  } else {
    print('WARNING: Missing tcov files')
  }
  
  ################ Per base sequence quality ################
  #read the data
  files <- list.files(path = FileDir, pattern="pbsq$")
  
  if (length(files)>0) {
    fnames = gsub("_fastqc.pbsq", "", files, perl=TRUE)
    labs = as.character(SampNames[fnames])
    
    # Create lists to hold coverage and cumulative coverage for each alignment,
    # and read the data into these lists.
    stats <- list()
    for (i in 1:length(files)) {
      stats[[i]] <- read.table(paste0(FileDir,'/',files[i]),stringsAsFactors=F)
      stats[[i]][,1]=StrRangeToNum(stats[[i]][,1])
    }
    
    # Save the graph to a file
    png(paste0(FilePrefix,"_pbsq_v",Version,".png"),width = 8.5, height = 6, units = "in", res=600)
    
    # Create plot area, but do not plot anything. Add gridlines and axis labels.
    plot(stats[[1]][,1],stats[[1]][,2], type='l', xlab="Base", ylab="Phred quality", ylim=c(0,42), main="Per base sequence quality",col=alpha(labsDF[labs[1],'cols'], AlphaCol))
    abline(h = 28, col = "green",lty=2)
    abline(h = 20, col = "red",lty=2)
    
    # Actually plot the data for each of the alignments (stored in the lists).
    for (i in 2:length(stats)) points(stats[[i]][,1],stats[[i]][,2], type='l', lwd=2, col=alpha(labsDF[labs[i],'cols'], AlphaCol))
    
    if(dim(ulabsDF)[1]>1) legend("bottomleft", legend=ulabsDF$factor, col=ulabsDF$cols, lty=1, lwd=4)
    
    dev.off()
    
    drawQCheatmap(stats,'Read length','Per base sequence quality',paste0(FilePrefix,"_heatmap_pbsq_v",Version,".png"),12,FactorsTable,labs)
    
    
  } else {
    print('WARNING: Missing pbsq files')
  }
  
  
  ################ Per sequence GC content ################
  #read the data
  files <- list.files(path = FileDir, pattern="psgc$")
  
  if (length(files)>0) {
    fnames = gsub("_fastqc.psgc", "", files, perl=TRUE)
    labs = as.character(SampNames[fnames])
    
    # Create lists to hold coverage and cumulative coverage for each alignment,
    # and read the data into these lists.
    stats <- list()
    for (i in 1:length(files)) {
      stats[[i]] <- read.table(paste0(FileDir,'/',files[i]),stringsAsFactors=F)
      stats[[i]][,2]=stats[[i]][,2]/sum(stats[[i]][,2]) * 100
    }
    
    # Save the graph to a file
    png(paste0(FilePrefix,"_psgc_v",Version,".png"), width = 8.5, height = 6, units = "in", res=600)
    
    # Create plot area, but do not plot anything. Add gridlines and axis labels.
    plot(stats[[1]][,1],stats[[1]][,2], type='l', xlab="Mean GC content [%]", ylab="Read percentage", ylim=c(0,9), main="GC distribution over all reads",
         col=alpha(labsDF[labs[1],'cols'], AlphaCol))
    
    # Actually plot the data for each of the alignments (stored in the lists).
    for (i in 2:length(stats)) points(stats[[i]][,1],stats[[i]][,2], type='l', lwd=2, col=alpha(labsDF[labs[i],'cols'], AlphaCol))
    if(dim(ulabsDF)[1]>1) legend("topleft", legend=ulabsDF$factor, col=ulabsDF$cols, lty=1, lwd=4)
    
    dev.off()
    
    drawQCheatmap(stats,'Mean GC content [%]','GC distribution over all reads',paste0(FilePrefix,"_heatmap_psgc_v",Version,".png"),16,FactorsTable,labs)
    
    
  } else {
    print('WARNING: Missing psgc files')
  }
  
  
  ################ Per sequence quality scores ################
  #read the data
  files <- list.files(path = FileDir, pattern="psqs$")
  
  if (length(files)>0) {
    fnames = gsub("_fastqc.psqs", "", files, perl=TRUE)
    labs = as.character(SampNames[fnames])
    
    stats <- list()
    gmax=0
    for (i in 1:length(files)) {
      stats[[i]] <- read.table(paste0(FileDir,'/',files[i]),stringsAsFactors=F)
      stats[[i]][,2]=stats[[i]][,2]/sum(stats[[i]][,2]) * 100
      tmax = max(stats[[i]][,2])
      if (tmax>gmax) {
        gmax = tmax
      }
    }
    
    MaxReadPrc = ceiling(gmax/10)*10
    
    # Save the graph to a file
    png(paste0(FilePrefix,"_psqs_v",Version,".png"), width = 8.5, height = 6, units = "in", res=600)
    
    # Create plot area, but do not plot anything. Add gridlines and axis labels.
    plot(stats[[1]][,1],stats[[1]][,2], type='l', ylim=c(0,MaxReadPrc), xlab="Mean phred quality", ylab="Read percentage", main="Per sequence quality scores",col=alpha(labsDF[labs[1],'cols'], AlphaCol))
    abline(v = 28, col = "green")
    abline(v = 20, col = "red")
    
    # Actually plot the data for each of the alignments (stored in the lists).
    for (i in 2:length(stats)) points(stats[[i]][,1],stats[[i]][,2], type='l', lwd=2, col=alpha(labsDF[labs[i],'cols'], AlphaCol))
    if(dim(ulabsDF)[1]>1) legend("topleft", legend=ulabsDF$factor, col=ulabsDF$cols, lty=1, lwd=4)
    dev.off()
    
    drawQCheatmap(stats,'Mean phred quality','Per sequence quality scores',paste0(FilePrefix,"_heatmap_psqs_v",Version,".png"),14,FactorsTable,labs)
  } else {
    print('WARNING: Missing psqs files')
  }
  
  
  ################ Sequence length distribution ################
  #read the data
  files <- list.files(path = FileDir, pattern="sldst$")
  
  if (length(files)>0) {
    fnames = gsub("_fastqc.sldst", "", files, perl=TRUE)
    labs = as.character(SampNames[fnames])
    
    stats <- list()
    gmax=0
    for (i in 1:length(files)) {
      stats[[i]] <- read.table(paste0(FileDir,'/',files[i]),stringsAsFactors=F)
      stats[[i]][,1]=StrRangeToNum(stats[[i]][,1])
      stats[[i]][,2]=stats[[i]][,2]/sum(stats[[i]][,2]) * 100
      tmax = max(stats[[i]][,2])
      if (tmax>gmax) {
        gmax = tmax
      }
    }
    
    MaxReadPrc = ceiling(gmax/10)*10
    
    # Save the graph to a file
    png(paste0(FilePrefix,"_sldst_v",Version,".png"), width = 8.5, height = 6, units = "in", res=600)
    
    # Create plot area, but do not plot anything. Add gridlines and axis labels.
    plot(stats[[1]][,1],stats[[1]][,2], type='l', ylim=c(0,MaxReadPrc), xlab="Read length", ylab="Read percentage", main="Sequence length distribution",col=alpha(labsDF[labs[1],'cols'], AlphaCol))
    
    # Actually plot the data for each of the alignments (stored in the lists).
    for (i in 2:length(stats)) points(stats[[i]][,1],stats[[i]][,2], type='l', lwd=2, col=alpha(labsDF[labs[i],'cols'], AlphaCol))
    if(dim(ulabsDF)[1]>1) legend("topleft", legend=ulabsDF$factor, col=ulabsDF$cols, lty=1, lwd=4)
    dev.off()
    
    drawQCheatmap(stats,'Read length','Sequence length distribution',paste0(FilePrefix,"_heatmap_sldst_v",Version,".png"),12,FactorsTable,labs)
  } else {
    print('WARNING: Missing sldst files')
  }
  
  ################ Adapter contamination ################
  #read the data
  files <- list.files(path = FileDir, pattern="adaptct$",full.names = TRUE)
  
  if (length(files)>0) {
    fnames = gsub("_fastqc.adaptct", "", files, perl=TRUE)
    labs = as.character(SampNames[fnames])
    
    read_adaptct = function(file) {
      sampleid = gsub('_fastqc.adaptct','',basename(file))
      read_tsv(file,col_types='cddddd') %>% 
        dplyr::rename(Position=`#Position`) %>%
        mutate(PositionNum = StrRangeToNum(Position),
               sampleid = sampleid)
    }
    
    data <- files %>%
      map(read_adaptct) %>% 
      reduce(rbind) 
      
    PlotData = data  %>%
      select(-Position) %>%
      pivot_longer(!c('PositionNum','sampleid'),values_to = "AdapterContamination",names_to = "AdapterType") %>%
      #mutate(AdapterMax = do.call(pmax, select(., matches("Adapter|Sequence")))) %>% 
      #mutate(AdapterMax = rowSums(select(., matches("Adapter|Sequence"))))  #rowSums is not good since the same read can be counted more than once
      #select(sampleid, PositionNum, AdapterMax) %>%
      mutate(PositionNum = as.numeric(PositionNum )) %>%
      left_join(SampleTable,by='sampleid') 
    
    
    LabelData <- PlotData %>% 
      group_by(sampleid, AdapterType) %>%
      filter(PositionNum == max(PositionNum)) 
    
    
    p = ggplot(PlotData,aes(x=PositionNum, y=AdapterContamination, group=interaction(sampleid,AdapterType), color=Platform)) + 
      geom_line(alpha=0.4) +theme_bw() +
      geom_label_repel(data=LabelData, aes(x=PositionNum, y= AdapterContamination, group=interaction(sampleid,AdapterType), label=sampleid), min.segment.length = 0, size=2, color="black",alpha=0.8)+
      xlab('Read position') +ylab('% of sequences') +
      scale_color_manual(values=mnmcolors) + 
      facet_wrap(~AdapterType,nrow=1)+
      theme(legend.position = "bottom")
    ggsave(p,file=paste0(FilePrefix,'_FastQC-Adaptct_v',Version,'.png'), scale=2, width=4, height=4, dpi=600)
    
    
  } else {
    print('WARNING: Missing adaptct files')
  }  
  
}
