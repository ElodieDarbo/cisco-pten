load_data <- function(corr.bins,frequencies){
    reps <- lapply(1:length(corr.bins),function(i,corr.bins,frequencies){
      message("Loading replicate ", i," ...")
      tmp.bins <- fread(corr.bins[i])
      tmp.freq <- fread(frequencies[i])
      frequencies <- data.table(chr1=tmp.bins$V1[tmp.freq$V1], 
                                locus1=tmp.bins$V2[tmp.freq$V1], 
                                chr2=tmp.bins$V1[tmp.freq$V2], 
                                locus2=tmp.bins$V2[tmp.freq$V2], 
                                frequencies=tmp.freq$V3)
      frequencies <- frequencies[locus1!=locus2]
      frequencies$int1 <-paste(frequencies$chr1,frequencies$locus1,sep='_')
      frequencies$int2 <-paste(frequencies$chr2,frequencies$locus2,sep='_')
      return(frequencies)
    },corr.bins,frequencies)
    
    message("Merging ...")
    merged <- rbindlist(reps)
    merged <- merged[,list(frequencies=sum(frequencies)),by=c("chr1","locus1","chr2","locus2","int1","int2")]
    all.data <- reps
    names(all.data) <- paste("rep",1:length(corr.bins),sep="_")
    all.data$merged <- merged
    return(all.data)
}

compute_binomial <- function(binned_df_filtered,smooth.dist=FALSE){
  message("Statistical analysis started ...")
  #all read pairs used in binomial
  numberOfReadPairs <- sum(binned_df_filtered$frequencies)
  
  #calculate coverage 
  all_bins <- unique(c(unique(binned_df_filtered$int1), unique(binned_df_filtered$int2)))
  all_bins <- sort(all_bins)
  
  covA <- binned_df_filtered[,sum(frequencies),by=int1]	
  covB <- binned_df_filtered[,sum(frequencies),by=int2]
  
  covA <- setkey(covA,key='int1')
  setnames(covB, 1,'int1')
  covB <- setkey(covB,key='int1')
  
  cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
  cov$V1.x[is.na(cov$V1.x)]=0
  cov$V1.y[is.na(cov$V1.y)]=0
  cov$coverage=cov$V1.x+cov$V1.y
  coverage=cov$coverage
  names(coverage)=cov$int1
  sumcov <- sum(coverage)
  relative_coverage <- coverage/sumcov
  names(relative_coverage)=names(coverage)
  binned_df_filtered$coverage_source <- relative_coverage[binned_df_filtered$int1]
  binned_df_filtered$coverage_target <- relative_coverage[binned_df_filtered$int2]
  
  numberOfAllInteractions <- length(all_bins)^2
  upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2
  
  chromos <- unique(binned_df_filtered$chr1)
  chrlens <- c()
  chrlens[chromos] <- length(unique(c(unique(binned_df_filtered$locus1),unique(binned_df_filtered$locus2))))
  
  cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
  
  diagonalProb <- sum(relative_coverage^2)
  
  binned_df_filtered$distance <- binned_df_filtered$locus2 - binned_df_filtered$locus1
  
  if (smooth.dist){
    binned_df_filtered$probabilityDistCorrection <- pgamma((binned_df_filtered$distance)/10000,shape = 1,scale = 1,lower.tail = F)
  }
  else {
    binned_df_filtered$probabilityDistCorrection <- 0
  }
  
  binned_df_filtered$probabilityOfInteraction <- binned_df_filtered$coverage_source*binned_df_filtered$coverage_target*2+binned_df_filtered$probabilityDistCorrection
  
  binned_df_filtered$predicted <- binned_df_filtered$probabilityOfInteraction * numberOfReadPairs
  
  message("Computing Binomial ...")
  
  binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
  {
    binom.test(as.numeric(x[["frequencies"]])-1, numberOfReadPairs, as.numeric(x[["probabilityOfInteraction"]]), alternative = "greater")$p.value
  }	
  )
  
  binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
  
  binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber)
  
  binned_df_filtered$sig <- -log10(binned_df_filtered$qvalue)
  binned_df_filtered[sig==Inf]$sig <- max(binned_df_filtered[sig!=Inf]$sig)
  return(binned_df_filtered)
}

plot_triangle <- function(binned_df_filtered,plot.name){
  pageCreate(width = 12, height = 10, default.units = "cm")
  
  plotText(label = plot.name, fontsize = 12,
           x = 1, y = 2, just = "left", default.units = "cm")
  
  params_c <- pgParams(chrom = "chr10", chromstart = 87652000, chromend = 88644000, 
                       assembly = "hg38",
                       x = 1, width = 10, default.units = "cm")
  hic_gm <- plotHicTriangle(data = binned_df_filtered[,list(locus1,locus2,sig)], params = params_c,
                            zrange = c(0, 10), resolution = 2000,
                            y = 8, height = 7, just = c("left", "bottom"))
  
  annoHeatmapLegend(plot = hic_gm, fontsize = 7, 
                  x = 10, y = 3, width = 0.25, height = 2, 
                  just = c("right", "top"), default.units = "cm")
  
  annoGenomeLabel(plot = hic_gm, params = params_c, 
                  scale = "Kb", fontsize = 7,
                  y = 8.25)
  pageGuideHide()
}

export.plot <- function (file.prefix="PlotExport",
                         export.formats="pdf", # supported: postscript, jpg, png, bmp, pdf
                         width=11, # in inches
                         height=8, # in inches
                         horizontal=T,
                         ... ## Additional parameters are passed to the export method
) {
  
  ppi <- 72
  file.ext <- c(
    postscript = "ps",
    pdf = "pdf",
    ps = "ps",
    eps = "eps",
    jpeg="jpg",
    jpg="jpg",
    bmp="bmp",
    png="png",
    svg="svg",
    tiff="tiff")
  for (f in export.formats) {
    from.dev <- dev.cur();
    
    file.name <- paste(file.prefix,file.ext[f], sep=".")
    
    if ((f == "postscript") || (f == "ps")) {
      postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal, ...)
    } else if (f == "eps") {
      postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal,onefile=F, ...)
    } else if (f == "pdf") {
      pdf(file.name, paper="special",width=width,height=height, ...)
    } else if ((f == "jpg") || (f == "jpeg")) {
      jpeg(file.name,width=(width*ppi),height=(height*ppi),quality=100, ...)
    } else if (f == "png") {
      png(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "bmp") {
      bitmap(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "svg") {
      svg(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "tiff") {
      #tiff(filename = "Rplot%03d.tiff", width = 480, height = 480, units = "px", pointsize = 12, compression = c("none", "rle", "lzw", "jpeg", "zip"), bg = "white", res = NA,  ..., type = c("cairo", "Xlib", "quartz"), antialias)
      tiff(file.name,width=width*ppi,height=height*ppi, compression = 'none', ...)
    }
    else {
      print(paste("Error: format ", f, " is not supported", sep=""))
      return()
    }
    to.dev <- dev.cur()
    dev.set(which=from.dev)
    dev.copy(which=to.dev)
    dev.set(which=to.dev)
    dev.off()
    dev.set(which=from.dev) ## This is required because dev.off() returns to the first, not the last, device
  }
}

merge.interactions <- function(x,g, print.debug=F){
  if (counter %in% seq(1,30000,500)){
    print(counter)
  }
  tmp <- as_data_frame(induced_subgraph(g,x), what = "vertices")
  ints <- data.table(crossing(int1 = as.numeric(tmp$name), int2 = as.numeric(tmp$name)))
  ints <- ints[int1!=int2 & int1<int2,list(start1=int1,end1=int1+2000, start2=int2, end2=int2+2000)]
  ints$distance <- ints$start2 - ints$end1
  ints$type <- "raw"
  anchors <- reduce(GRanges("chr10",IRanges(as.numeric(tmp$name),as.numeric(tmp$name)+2000)),min.gapwidth=2000)
  if (length(anchors)>1){
    int.merge <- data.table(crossing(int1 = start(anchors), int2 = start(anchors)))
    int.merge <- int.merge[int1!=int2 & int1<int2]
    int.merge <- int.merge[,list(start1=int1,end1=end(anchors)[match(int1,start(anchors))],
                                 start2=int2,end2=end(anchors)[match(int2,start(anchors))])]
    int.merge$distance <- int.merge$start2 - int.merge$end1
    int.merge$type <- "merged"
    ints <- rbind(ints,int.merge)
  }
  ints$width1 <- ints$end1 - ints$start1
  ints$width2 <- ints$end2 - ints$start2
  ints$clique <- counter
  if (counter %in% seq(1,30000,500) & print.debug){
    print(ints)
  }
  counter <<- counter + 1
  return(list(raw=ints[type=="raw"],merged=ints[type=="merged"]))
}

analyze.replicates <- function(binned_df_binomial,main.dir,experience,plot=T){
  binned_df_binomial <- lapply(binned_df_binomial,function(x){
    x$ID <- paste(x$chr1,x$locus1,x$locus2,sep="_")
    return(x)
  })
  
  if (length(binned_df_binomial)<4){
    message("Identify 2 replicates")
    common <- sum(binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID)
    rep1 <- sum(!binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID)
    rep2 <- sum(!binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_1[sig>2]$ID)
    
    commonIDs <- binned_df_binomial$rep_1[sig>2][binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID]$ID
    rep1IDs <- binned_df_binomial$rep_1[sig>2][!binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID]$ID
    rep2IDs <- binned_df_binomial$rep_2[sig>2][!binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_1[sig>2]$ID]$ID
    
    if (plot){
      #specify values to use in venn diagram
      fit <- euler(c('rep1' = rep1, 'rep2' = rep2, 'rep1&rep2' = common))
      
      #create venn diagram with custom colors
      print(plot(fit, fill=c('coral2', 'steelblue'),quantities=T))
      export.plot(file.path(main.dir,experience,paste0(experience,"_rep_venn")),width=4,height=4)
    }
    
    comp <- merge(binned_df_binomial$rep_1[,list(ID,chr1,locus1,locus2,distance,rep1=sig)],binned_df_binomial$rep_2[,list(ID,chr1,locus1,locus2,distance,rep2=sig)],by=c("ID","distance","chr1","locus1","locus2"))
    comp$mean.sig <- apply(comp[,list(rep1,rep2)],1,mean)
    
    if (plot){
      g <- ggscatter(comp, x="rep1", y="rep2", size = 0.5, alpha=0.5, add="reg.line") +
        geom_smooth(method="lm",linewidth=0.5,fill=NA,color="blue") +
        theme_light() +
        theme(legend.position = "right") +
        stat_cor(label.y = 300) +
        theme(aspect.ratio = 1)
      print(g)
      export.plot(file.path(main.dir,experience,paste0(experience,"_rep_correlation")),width=7,height=7,export.formats = c("pdf","png"))
    }
    
  } else if (length(binned_df_binomial)==4) {
    message("Identify 3 replicates")
    common <- sum(binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID & binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID)
    rep1 <- sum(!binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID)
    rep2 <- sum(!binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_1[sig>2]$ID & !binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID)
    rep3 <- sum(!binned_df_binomial$rep_3[sig>2]$ID%in%binned_df_binomial$rep_1[sig>2]$ID & !binned_df_binomial$rep_3[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID)
    
    commonIDs <- binned_df_binomial$rep_1[sig>2][binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID & binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID]$ID
    rep1IDs <- binned_df_binomial$rep_1[sig>2][!binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID]$ID
    rep2IDs <- binned_df_binomial$rep_2[sig>2][!binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_1[sig>2]$ID & !binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID]$ID
    rep3IDs <- binned_df_binomial$rep_3[sig>2][!binned_df_binomial$rep_3[sig>2]$ID%in%binned_df_binomial$rep_1[sig>2]$ID & !binned_df_binomial$rep_3[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID]$ID
    
    rep12 <- sum(binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%commonIDs)
    rep13 <- sum(binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%commonIDs)
    rep23 <- sum(binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID & !binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_1[sig>2]$ID & !binned_df_binomial$rep_2[sig>2]$ID%in%commonIDs)
    
    rep12IDs <- binned_df_binomial$rep_1[sig>2][binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%commonIDs]$ID
    rep13IDs <- binned_df_binomial$rep_1[sig>2][binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%binned_df_binomial$rep_2[sig>2]$ID & !binned_df_binomial$rep_1[sig>2]$ID%in%commonIDs]$ID
    rep23IDs <- binned_df_binomial$rep_2[sig>2][binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_3[sig>2]$ID & !binned_df_binomial$rep_2[sig>2]$ID%in%binned_df_binomial$rep_1[sig>2]$ID & !binned_df_binomial$rep_2[sig>2]$ID%in%commonIDs]$ID
    
    if (plot){
      #specify values to use in venn diagram
      fit <- euler(c('rep1' = rep1, 'rep2' = rep2, 'rep3' = rep3, 'rep1&rep2' = rep12, 'rep1&rep3' = rep13, 'rep2&rep3' = rep23 , 'rep1&rep2&rep3' = common))
      
      #create venn diagram with custom colors
      print(plot(fit, fill=c('coral2', 'steelblue','gold2'),quantities=T))
      export.plot(file.path(main.dir,experience,paste0(experience,"_rep_venn")),width=4,height=4)
    }
    
    comp <- merge(binned_df_binomial$rep_1[,list(ID,chr1,locus1,locus2,distance,rep1=sig)],binned_df_binomial$rep_2[,list(ID,chr1,locus1,locus2,distance,rep2=sig)],by=c("ID","distance","chr1","locus1","locus2"))
    comp <- merge(comp,binned_df_binomial$rep_3[,list(ID,chr1,locus1,locus2,distance,rep3=sig)],by=c("ID","distance","chr1","locus1","locus2"))
    comp$mean.sig <- apply(comp[,list(rep1,rep2,rep3)],1,mean)
    comp.facet <- rbind(comp[,list(ID,comparison="rep1&rep2",score1=rep1,score2=rep2)],
                        comp[,list(ID,comparison="rep1&rep3",score1=rep1,score2=rep3)],
                        comp[,list(ID,comparison="rep2&rep3",score1=rep2,score2=rep3)])
    if (plot){
      g <- ggscatter(comp.facet, x="score1", y="score2", size = 0.5, alpha=0.5, add="reg.line") +
        geom_smooth(method="lm",linewidth=0.5,fill=NA,color="blue") +
        theme_light() +
        facet_wrap(~comparison, ncol=3) +
        theme(legend.position = "right") +
        stat_cor(label.y = 300) +
        theme(aspect.ratio = 1)
      print(g)
      export.plot(file.path(main.dir,experience,paste0(experience,"_rep_correlation")),width=12,height=4,export.formats = c("pdf","png"))
      
    }
  }
  return(list(comp=comp,commonIDs=commonIDs))
}

