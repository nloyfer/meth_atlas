## To run:
## Rscript process_array.R {idat directory} {output csv path} {path to ref_sample.RData}
## Note: Requires minfi package.
## To install minfi:
## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite("minfi")

args<-commandArgs(trailingOnly=TRUE)
if (!length(args)==3){
	stop('Usage: Rscript process_array.R {idat directory} {output csv path} {path to ref_sample.RData}')
}


library(minfi)
get_betas <- function(grSet,RGChannelSet,p_threshold=0.01,min_beads=3,sex_remove=T){
  # filter by P-value
  betas <- getBeta(grSet)
  pval <- detectionP(RGChannelSet)[rownames(betas),]
  betas[pval>p_threshold] <- NA
  # filter sex chromosomes
  if (sex_remove){
    chrom <- as.character(seqnames(granges(grSet)))
    chrom.sex <- chrom=='chrX' | chrom=='chrY'
    betas <- betas[!chrom.sex,]
  }
  #filter by bead number
  nbeads <- getNBeads(RGChannelSet)
  t1 <- getProbeInfo(RGChannelSet,'I')
  t2 <- getProbeInfo(RGChannelSet,'II')
  t1.nbeads.u <- nbeads[rownames(nbeads) %in% t1[,2],]; rownames(t1.nbeads.u) <- t1[match(rownames(t1.nbeads.u),t1[,2]),1]
  t1.nbeads.m <- nbeads[rownames(nbeads) %in% t1[,3],]; rownames(t1.nbeads.m) <- t1[match(rownames(t1.nbeads.m),t1[,3]),1]
  t1.nbeads.m <- t1.nbeads.m[rownames(t1.nbeads.u),]
  t1.nbeads.min <- pmin(t1.nbeads.u,t1.nbeads.m)
  t2.nbeads <- nbeads[rownames(nbeads) %in% t2[,2],]; rownames(t2.nbeads) <- t2[match(rownames(t2.nbeads),t2[,2]),1]
  nbeads.array <- rbind(t1.nbeads.min,t2.nbeads)
  cgs <- intersect(rownames(nbeads.array),rownames(betas))
  nbeads.array <- nbeads.array[cgs,]
  betas <- betas[cgs,]
  betas[nbeads.array<min_beads] <- NA
  return(betas)
}

process_array <- function(array_dir, out_path, ref_sample_path=args[3])
{
	array_path <- array_dir
	
	print(paste0("1. process_array script - READING ARRAY - ",Sys.time()))
	rgSet <- read.metharray.exp(array_path,extended=T)
	
	# Get p-values	
	print(paste0("2. process_array script - P_VALUES CALCULATION - ",Sys.time()))
	pVal <- detectionP(rgSet, type = "m+u")

	# illumina normalization
	## join with reference sample, normalize, then remove reference sample
	print(paste0("3. process_array script - ILLUMINA - ",Sys.time()))
	load(ref_sample_path)
	ref.sample@colData <- rgSet@colData[1,]; rownames(ref.sample@colData)[1] <- 'REF'
	rgSet <- combineArrays(ref.sample,rgSet)
	MethSet <- preprocessIllumina(rgSet)
	
	# Context removal
	print(paste0("4. process_array script - CONTEXT_REMOVAL - ",Sys.time()))
	MethDrop <- dropMethylationLoci(MethSet, dropRS = TRUE, dropCH = TRUE)

	# Get RatioSet
	print(paste0("5. process_array script - CONVERT TO RATIOSET - ",Sys.time()))
	RSet <- ratioConvert(MethDrop, what = "beta", keepCN = FALSE,offset=100)

	# Get GenomicRatioSet
	print(paste0("6. process_array script - CONVERT TO GENOMICRATIOSET - ",Sys.time()))
    grSet <- mapToGenome(RSet)
	
	# Remove SNPs
	print(paste0("7. process_array script - REMOVE SNPs - ",Sys.time()))
    grSet <- dropLociWithSnps(grSet, snps=c("SBE","CpG"), maf=0.01)
	
	# Get filtered beta values
	print(paste0("8. process_array script - BETA VALUES CALCULATION - ",Sys.time()))
	betas<-get_betas(grSet,rgSet)
	
	#remove reference sample
	betas <- betas[,2:ncol(rgSet)]

	print(paste0("13. process_array script - WRITING ARRAY TO CSV - ",Sys.time()))
	if (length(sampleNames(rgSet)) == 2){
		betas <- as.data.frame(betas)
		colnames(betas)<-sampleNames(rgSet)[2]
	}
	write.csv(round(betas,5),out_path,quote=FALSE)
	print(paste0("14. process_array script - FINISHED - ",Sys.time()))	
}

process_array(args[1], args[2])
