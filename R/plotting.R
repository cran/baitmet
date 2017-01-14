plotMZ <- function(Experiment, AlignId, mz.ind=1, per.class=T, aligned=FALSE, xlim=NULL)
{	
	if(!(any(unlist(lapply(Experiment@Data@FeatureList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("MZ have to be quantified first")
	
	if(nrow(Experiment@MetaData@Phenotype)==0) 
		if(per.class==T) {per.class=F; warning("The experiment does not contain phenotypic metadata. Profiles are shown per sample.")}
	
	empty.samples <- which(lapply(Experiment@Data@FeatureList,nrow)==0)
	if(length(empty.samples)!=0)
	{
		Experiment@Data@FeatureList <- Experiment@Data@FeatureList[-empty.samples]
		Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	}
	
	alignId <- lapply(Experiment@Data@FeatureList,function(x){x$AlignID})
	N.groups <- max(unique(as.numeric(as.vector(unlist(alignId)))))	
	N.samples <- length(Experiment@Data@FeatureList)
	
	samples.name <- names(Experiment@Data@FeatureList)		
	for(i in 1:N.samples) samples.name[i] <- strsplit(as.character(samples.name[i]), split="\\.")[[1]][1]
	
	profile.list <- lapply(Experiment@Data@FeatureList,function(x) {
			outp <- as.character(x[which(as.numeric(as.vector(x$AlignID))==AlignId),"Profile"])[mz.ind]
			outMZ <- x[which(as.numeric(as.vector(x$AlignID))==AlignId),"mz"][mz.ind]
			time <- NA
			int <- NA
			if(length(outp)!=0)
			{
				output <- erah:::sparse.to.vector(outp)
				time <- output$time
				int <- output$int*as.numeric(as.character(x[which(as.numeric(as.vector(x$AlignID))==AlignId),"Int"]))[mz.ind]
			}
			list(time=time,int=int, mz=outMZ)
			})
	
	profile.len <- max(unlist(lapply(unlist(profile.list, recursive=F),length)))
	
	profile.time = matrix(NA, nrow=profile.len,ncol=N.samples)
	profile.int = matrix(NA, nrow=profile.len,ncol=N.samples)

	for(i in 1:N.samples) {profile.time[1:length(profile.list[[i]]$time),i] <- profile.list[[i]]$time; profile.int[1:length(profile.list[[i]]$int),i] <- profile.list[[i]]$int}
	
	na.samples.i <- which(apply(apply(profile.int,2,is.na),2,all)==T)
	na.samples.t <- which(apply(apply(profile.time,2,is.na),2,all)==T)
	na.samples <- unique(na.samples.i,na.samples.t)
	
	if(length(na.samples)!=0)
	{ 
		samples.name <- samples.name[-na.samples]
		profile.int <- profile.int[,-na.samples]
		profile.time <- profile.time[,-na.samples]
	}
	
	if(aligned==TRUE) {
		vector.time <- diag(profile.time[apply(profile.int,2,which.max),])
		time.mean <- mean(vector.time)	
		profile.time <- sweep(profile.time,2,(vector.time-time.mean),"-") 
	}
	
	compound.name <- as.character(Experiment@Results@Identification[which(Experiment@Results@Identification$AlignID==AlignId),"Name.1"])
	compound.mz <- as.numeric(as.vector(profile.list[[1]]$mz))
#Experiment@Data@FeatureList[which(x$AlignID==AlignId),"mz"][mz.ind]
	if(per.class==F)
	{
		
		matplot(profile.time, profile.int, type="l", lty=1, col=(1:length(samples.name)), main=paste(compound.name, "\n m/z:", compound.mz), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.name, pch=19, col=(1:length(samples.name)), title="Samples")
		par(font=1)
	}else{
		pn <- Experiment@MetaData@Phenotype
		indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
		class.names <- pn[indx,"class"]
		
		samples.class.type <- levels(pn$class)
			
		matplot(profile.time,profile.int, type="l", lty=1, col=class.names, main=paste(compound.name, "\n m/z:", compound.mz), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
		par(font=1)

	}
	
}

