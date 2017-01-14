globalVariables("mslib")
## Intern Classes for BaitMet:
#data(mslib, package="erah", envir=environment())
suppressForeignCheck("mslib")

setClass(Class = "MetaData", representation = representation(Instrumental = "data.frame", Phenotype = "data.frame", DataDirectory="character"))

setClass(Class = "Statistics", representation = representation(Univariate="data.frame", Multivariate="data.frame"))	

setClass(Class="MSResultsParameters", representation = representation(Alignment = "list", Identification = "list"))

setClass(Class="Data", representation = representation(FeatureList = "list", FactorList = "list", Parameters = "list"))

setClass(Class = "Results", representation = representation(Parameters="MSResultsParameters", Alignment = "data.frame", Identification="data.frame", Statistics="Statistics"))

setClass(Class="MetaboSet",representation= representation(Info = "character", Data="Data", MetaData="MetaData", Results = "Results"))
	
setClass(Class = "BaitMetSoftParameters", representation = representation(algorithm="character", ri.error = "numeric", min.peak.width = "numeric", min.peak.height = "numeric", noise.threshold = "numeric", avoid.processing.mz = "vector", compression.coef = "numeric", analysis.time="vector"))

#setClass(Class = "Baitmet_DB", representation = representation(name="character", version="character", info="character", database="list"))


## Main Software functions:

#normalize <- osd::normalize
#paste.sp <- erah::paste.sp
#is.even <- erah::is.even

normalize <- function (x) 
{
    x[is.na(x)] <- 0
    if (is.matrix(x) == T) 
        norm.x <- sweep(x, 2, apply(x, 2, function(k) max(k, 
            na.rm = T)), "/")
    if (is.matrix(x) == F) 
        norm.x <- x/max(x, na.rm = T)
    norm.x[is.na(norm.x)] <- 0
    norm.x
}

is.even <- function (x) { x%%2 == 0 } 

paste.sp <- function (x, y) {
    paste(x, y, sep = ",")
    }


setBaitPar <- function(ri.error=0.05, min.peak.width, min.peak.height=500, noise.threshold=500, avoid.processing.mz=c(1:69,73:75,147:149), compression.coef=2, analysis.time=0)
{
	softPar <- new("BaitMetSoftParameters",algorithm="BaitMet", ri.error = ri.error, min.peak.width = min.peak.width/60, min.peak.height = min.peak.height, noise.threshold = noise.threshold, avoid.processing.mz = avoid.processing.mz, compression.coef = compression.coef, analysis.time=analysis.time)
	softPar
}


setChrmMethod <- function(method=c("alk.var5","alk.mdn35","fame.var5","fame.mdn35"), rt, ri, name="ChromMethod #1")
{

	method <- match.arg(method, c("alk.var5","alk.mdn35","fame.var5","fame.mdn35"))
		
	if(length(rt)!=length(ri)) stop("The RT and RI vector lengths are different")
	
	list(name=name, method=method, ref.rt=rt, ref.ri=ri)	
}

decBaitMet <- function(Experiment, BaitParameters, ms.library=mslib, chrom.method, samples.to.process=NULL)
{
	
	home.library=NULL
	softParameters <- BaitParameters
	
	Number.of.Samples <- nrow(Experiment@MetaData@Instrumental)
	if(is.null(samples.to.process)) samples.to.process <- 1:Number.of.Samples
	stopifnot(samples.to.process>=1, max(samples.to.process)<=Number.of.Samples, length(samples.to.process)<=Number.of.Samples)
	
	soft.par <- list(processed.with="BaitMet", scans.per.second=0, ri.error = softParameters@ri.error, min.peak.width = softParameters@min.peak.width, min.peak.height = softParameters@min.peak.height, noise.threshold = softParameters@noise.threshold, avoid.processing.mz = softParameters@avoid.processing.mz, compression.coef = softParameters@compression.coef, analysis.time = softParameters@analysis.time)
	Experiment@Data@Parameters <- soft.par
	
	ri.db <- sapply(1:length(ms.library@database), function(i) ms.library@database[[i]]$RI.VAR5.ALK)
	ref.mat <- sapply(1:length(ms.library@database), function(i)erah:::convertMSPspectra.dot(ms.library@database[[i]]$Spectra,2000))
		
	k <- 1
	for(index in samples.to.process)
	{
		cat("\n Searching metabolites in",as.character(Experiment@MetaData@Instrumental$filename[index]),"... Processing", k,"/",length(samples.to.process),"\n")  
		Experiment <- processSample.baitmet(Experiment, index, chrom.method, ri.db, ref.mat)
		k <- k + 1
	}
	cat("\n Grouping metabolites across samples... \n")
	Experiment <- alignBaitMet(Experiment, chrom.method)
	cat("\n Computing spectral match scores... \n")
	Experiment <- identifyBaitMet(Experiment, id.database=ms.library)
	cat("\n Done! \n")
	Experiment	
}

alignBaitMet <- function(Experiment, chrom.method)
{
	#Experiment <- ex
	al.par <- list(alignment.algorithm="BaitMet", min.spectra.cor=0, max.time.dist=0, mz.range=1:600, chrom.method=chrom.method)
	Experiment@Results@Parameters@Alignment <- al.par
	
	#Experiment@Results@Parameters@Alignment
	#str(c(u, chrom.method))
	
	#min.spectra.cor <- Experiment@Results@Parameters@Alignment$min.spectra.cor
	#max.time.dist <- Experiment@Results@Parameters@Alignment$max.time.dist
	#mz.range <- Experiment@Results@Parameters@Alignment$mz.range
	#maxMZ <- max(mz.range)
	Experiment@Data@FactorList <- align.across.samples(Experiment@Data@FactorList)
	Experiment@Results@Alignment <- create.factorlist.table(Experiment)
	
	Experiment
}


identifyBaitMet<- function(Experiment, id.database)
{
	Experiment <- identifyComp_baitmet(Experiment, id.database)
	Experiment	
}











