# Charger les bibliothèques nécessaires
library(GenomicRanges)

# Charger les packages
library(readr)
library(dplyr)
library(data.table)
library(MASS)

setwd("../data")


# Lire le fichier BED avec read.table en spécifiant les types de colonnes
bed_data <- fread("merged_file.PTEN_TAD.txt")

bed_data <- bed_data[,list(score=max(score),name=name[which.max(score)]),by=c("chr","start","end","tissue")]

# Convertir en GRanges avec métadonnées
atac_gr <- GRanges(seqnames = bed_data$chr,
                   ranges = IRanges(start = bed_data$start, end = bed_data$end),
                   name = bed_data$name,
                   score = bed_data$score,
                   tissue=bed_data$tissue)

atac_gr <- split(atac_gr,atac_gr$tissue)
atac_gr <- lapply(atac_gr,function(gr){
  red <- reduce(gr,with.revmap=T)
  red$tissue <- unique(gr$tissue)
  red$score <- lapply(red$revmap,function(x,gr) { 
    max(gr$score[x])
    },gr)
  return(red)
})

atac_gr <- unlist(as(atac_gr,"GRangesList"))

atac_reduce <- reduce(atac_gr, with.revmap=T)
atac_reduce$lr <- width(atac_reduce)
atac_reduce$nbr <- unlist(lapply(atac_reduce$revmap,length))

plot(start(atac_reduce),atac_reduce$nbr,pch=16, cex=0.5, col="blue")
lines(start(atac_reduce),atac_reduce$nbr,pch=16, cex=0.5, col="blue")

atac_reduce$nbr <- unlist(lapply(atac_reduce$revmap,function(x,atac_gr){
  nb <- length(unique(atac_gr[x]$tissue))
},atac_gr))

params.nb <- fitdistr(atac_reduce$nbr - 1,densfun = "negative binomial")
r.nb <- rnegbin(n=length(atac_reduce),mu = params.nb$estimate["mu"], theta = params.nb$estimate["size"])

params.pois <- fitdistr(atac_reduce$nbr - 1,densfun = "Poisson")
r.pois <- rpois(n=length(atac_reduce),lambda =  params.pois$estimate["lambda"])

params.geom <- fitdistr(atac_reduce$nbr - 1,densfun = "geometric")
r.geom <- rgeom(n=length(atac_reduce),prob =  params.geom$estimate["prob"])

atac_reduce$pVal.pois <- 1 - ppois(atac_reduce$nbr - 1,lambda =  params.pois$estimate["lambda"])
atac_reduce$pVal.nb <- 1 - pnbinom(atac_reduce$nbr - 1,size = params.nb$estimate["size"], mu = params.nb$estimate["mu"])
atac_reduce$pVal.geom <- 1 - pgeom(atac_reduce$nbr - 1,prob = params.geom$estimate["prob"])

# Compute AIC
aic_nb <- AIC(params.nb)
aic_geo <- AIC(params.geom)
aic_pois <- AIC(params.pois)

# Print AIC values
cat("AIC for Negative Binomial:", aic_nb, "\n")
cat("AIC for Geometric:", aic_geo, "\n")

# Perform goodness-of-fit tests (e.g., Chi-Square Test)
observed <- table(cut(atac_reduce$nbr-1, breaks = 56)) # Observed frequencies
names(observed) <- 1:56
expected_nb <- dnbinom(as.numeric(names(observed)), size = params.nb$estimate["size"], mu = params.nb$estimate["mu"]) * sum(observed)
expected_geo <- dgeom(as.numeric(names(observed)), prob = params.geom$estimate["prob"]) * sum(observed)
expected_pois <- dpois(as.numeric(names(observed)),lambda =  params.pois$estimate["lambda"]) * sum(observed)

# Chi-Square Test
chisq_nb <- sum((observed - expected_nb)^2 / expected_nb)
chisq_geo <- sum((observed - expected_geo)^2 / expected_geo)
chisq_pois <- sum((observed - expected_pois)^2 / expected_pois)

cat("Chi-Square for Negative Binomial:", chisq_nb, "\n")
cat("Chi-Square for Geometric:", chisq_geo, "\n")

par(mfrow=c(2,2))
hist(atac_reduce$nbr,breaks=55, main="number cell type per ATAC")
abline(v=min(atac_reduce$nbr[atac_reduce$pVal.nb<0.05]),col="red",lty=2)
text(x=min(atac_reduce$nbr[atac_reduce$pVal.nb<0.05])+1, y= 120, col="red", label=sum(atac_reduce$pVal.nb<0.05))
abline(v=min(atac_reduce$nbr[atac_reduce$pVal.pois<0.05]),col="green",lty=2)
text(x=min(atac_reduce$nbr[atac_reduce$pVal.pois<0.05])+1, y= 120, col="green", label=sum(atac_reduce$pVal.pois<0.05))
abline(v=min(atac_reduce$nbr[atac_reduce$pVal.geom<0.05]),col="blue",lty=2)
text(x=min(atac_reduce$nbr[atac_reduce$pVal.geom<0.05])+1, y= 120, col="blue", label=sum(atac_reduce$pVal.geom<0.05))
hist(r.nb,breaks=55,main="Negative binomial fitting",border = "red")
text(x = 30, y=100,label=paste("AIC:",round(aic_nb), ", ChiSq:",round(chisq_nb)))
hist(r.pois,breaks=55,main="Poisson fitting",border = "green")
text(x = 10, y=50,label=paste("AIC:",round(aic_pois), ", ChiSq:",round(chisq_pois)))
hist(r.geom,breaks=55,main="Geometric fitting",border = "blue")
text(x = 20, y=50,label=paste("AIC:",round(aic_geo), ", ChiSq:",round(chisq_geo)))
par(mfrow=c(1,1))

write.table(as.data.frame(atac_reduce)[,c("seqnames","start","end","nbr","pVal.nb")],"NegBinom_all_ATAC.bed", sep="\t", col.names=F,row.names=F, quote=F)

