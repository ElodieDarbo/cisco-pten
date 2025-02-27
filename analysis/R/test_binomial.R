library(data.table)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(plotgardener)
library(eulerr)
library(igraph)
library(GenomicInteractions)
library(tidyr)

source("functions_for_binomial.R")

main.dir <- "tiled-C"
# the name of the folder containing the hic-pro results
experience <- "GM12878"
data.dir <- file.path(main.dir,experience,"hic_results/matrix")

if (file.exists(file.path(main.dir,experience,"binomial_all.RData"))){
  load(file.path(main.dir,experience,"binomial_all.RData"))
} else {
  replicates <- list.dirs(data.dir,recursive = F)
  replicates <- replicates[!grepl("merge|ipynb",replicates)]
  
  corr.bins <- file.path(replicates,paste0("raw/2000/",sub("-","",basename(replicates)),"_2000_abs.bed"))
  frequencies <- file.path(replicates,paste0("raw/2000/",sub("-","",basename(replicates)),"_2000.matrix"))
  
  
  # Load all replicates, keep separated and compute merging
  binned_df <- load_data(corr.bins,frequencies)
  binned_df_binomial <- lapply(binned_df,compute_binomial)
  save(binned_df_binomial,file=file.path(main.dir,experience,"binomial_all.RData"))
}


# plot triangles heatmaps
for (r in names(binned_df_binomial)){
  print(r)
  print(plot_triangle(binned_df_binomial[[r]],paste0(experience," replicates: ",r)))
  export.plot(file.path(main.dir,experience,paste0(experience,"_",r,"_garden")),width=5,height=4)
}

### Compare replicates
comp.results <- analyze.replicates(binned_df_binomial,main.dir,experience,plot=T)

## Pour K562 et Francis, j'ai enlever le réplicat 2, voir plot venn et correlation
if (experience %in% c("K562", "Francis")){
  binned_df_binomial_2 <- binned_df_binomial[c("rep_1","rep_3")] 
  names(binned_df_binomial_2)[2] <- "rep_2"
  comp.results <- analyze.replicates(binned_df_binomial_2,main.dir,experience,plot=F)
} else if(experience %in% c("HCC-1937")){ # may be replicate 1 is not to keep for HCC-1937
  binned_df_binomial_2 <- binned_df_binomial[c("rep_1","rep_2")] 
  comp.results <- analyze.replicates(binned_df_binomial_2,main.dir,experience,plot=F)
}
##

comp <- comp.results$comp
commonIDs <- comp.results$commonIDs

#ici ça serait peut-être pas mal de filtrer sur la distance
bed.pe <- data.table(comp[ID%in%commonIDs & distance>0,list(chr1,start1=locus1,end1=locus1+2000,chr2=chr1,start2=locus2,end2=locus2+2000,mean.sig)])
write.table(bed.pe,file.path(main.dir,experience,paste0(experience,"_reproductible_ints.bedpe")),sep="\t",quote=F, col.names=F, row.names = F)


