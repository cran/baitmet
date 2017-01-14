align.across.samples <- function(factors.list)
{
	#factors.list <- Smp.list
	#factors.list <- ex@Data@FactorList
	
	empty.samples <- which(lapply(factors.list,nrow)==0)
	factors.list.original <- NULL
	if(length(empty.samples)!=0)
	{
		factors.list.original <- factors.list
		factors.list <- factors.list[-empty.samples] 
	}

	
	if(!(any(unlist(lapply(factors.list,function(x) {is.null(x$AlignID)}))==FALSE)))
	{	
		factors.list <- lapply(factors.list, function(x){
			outp <- cbind(x,matrix(0,nrow=length(x$ID)))
			colnames(outp)[ncol(outp)] <- "AlignID"
			outp
			})
	}else{
		factors.list <- lapply(factors.list, function(x){
			x$AlignID <- rep(0,nrow=length(x$ID))
			x
			})
	}
	
	#str(factors.list[[1]],1)
	
	DBID.list <- lapply(factors.list, function(x) as.numeric(as.vector(x$DB.Id)))
	existing.DBID <-unique(unlist(DBID.list))
	
	AlId.vect <- 1:length(existing.DBID)
	
	for(i in 1:length(factors.list))
	{
		SubInds <- sapply(DBID.list[[i]], function(x) which(existing.DBID==x))
		factors.list[[i]]$AlignID <- AlId.vect[SubInds]
	}
			
	if(!is.null(factors.list.original))
	{
		factors.list.original[-empty.samples] <- factors.list
		factors.list <- factors.list.original
	}		
		
	factors.list
}



create.factorlist.table <- function(object)
{
	empty.samples <- which(lapply(object@Data@FactorList,nrow)==0)
	if(length(empty.samples)!=0) object@Data@FactorList <- object@Data@FactorList[-empty.samples]
	
	factors.list <- object@Data@FactorList
	
	#Spectra table
		
	alignId <- lapply(factors.list,function(x){x$AlignID})
	N.groups <- unique(unlist(alignId))
	if(length(which(N.groups==0))!=0) N.groups <- N.groups[-which(N.groups==0)]
	N.samples <- length(object@Data@FactorList)
	samples.name <- names(object@Data@FactorList)		
	max.mz <- max(object@Results@Parameters@Alignment$mz.range)
	
	spectra.list <- apply(as.matrix(N.groups),1,function(Ng){
		#cat(Ng, " \n")
		align.iterator <- as.vector(which(lapply(alignId,function(x){
			if(length(intersect(Ng,x))!=0) {
				return(T)
				}else{
					return(F)
				}
			})==T))
			group.spectra <- lapply(as.matrix(align.iterator), function(x){ 
				list(spectra=erah:::convertMSPspectra(factors.list[[x]][which(alignId[[x]]==Ng),"Spectra"],max.mz),time=factors.list[[x]][which(alignId[[x]]==Ng),"RT"]) 
				})
			if(length(align.iterator)!=0)
			{	
				common.spectra <- normalize(rowSums(do.call(cbind,lapply(group.spectra,function(x){normalize(x$spectra)})))) 
				common.time <- mean(do.call(cbind,lapply(group.spectra,function(x){x$time})))
				return(list(spectra=common.spectra,time=common.time, nsamples=length(align.iterator)))
			}else{return(NULL)}
	})
	
	#list.delete <- unlist(lapply(spectra.list, is.null))
	#if(any(list.delete)) spectra.list <- spectra.list[-which(list.delete==T)]
	#N.groups <- length(spectra.list)	

	spectra.matrix <- do.call(cbind,lapply(spectra.list,function(x){x$spectra}))
	#time.vector <- unlist(lapply(spectra.list,function(x){x$time}))
	#foundin.vector <- unlist(lapply(spectra.list,function(x){x$nsamples}))	
	

	for(i in 1:N.samples) samples.name[i] <- strsplit(as.character(samples.name[i]), split="\\.")[[1]][1]
	
	align.matrix <- t(apply(as.matrix(N.groups),1,function(i){
		local.index <- as.vector(unlist(lapply(object@Data@FactorList,function(x) {
			outp <- as.character(x[which(x$AlignID==i),"Peak Height"])
			if(length(outp)==0) outp <- 0
			as.numeric(as.vector(outp))
			})))
	}))

	
	foundIn.vector <- apply(align.matrix,1,function(x) length(which(x!=0)))
	
	time.align.matrix <- t(apply(as.matrix(N.groups),1,function(i){
		local.index <- as.vector(unlist(lapply(object@Data@FactorList,function(x) {
			outp <- as.character(x[which(x$AlignID==i),"RT"])
			if(length(outp)==0) outp <- 0
			outp
			})))
	}))
	
	class(time.align.matrix) <- "numeric"
	rt.mean <- rowSums(time.align.matrix)/apply(time.align.matrix,1,function(x) length(which(x!=0)))
	
	align.List <- matrix(0,nrow=nrow(align.matrix),ncol=(5+N.samples))
	colnames(align.List) <- c("AlignID","Factor","tmean","FoundIn","Spectra",as.character(samples.name))

	align.List[,6:(5+N.samples)] <- align.matrix
	
	# for(i in 1:N.groups) 
	# {
	#	local.pos <- which(object@Results@Identification$AlignID==i)
	#	align.List[i,"AlignID"] <- object@Results@Identification[local.pos,"AlignID"]
	#	align.List[i,"Name"] <- as.character(object@Results@Identification[local.pos,"Name"])
	#  	align.List[i,"tmean"] <- object@Results@Identification[local.pos,"tmean"]
	# 	#alignList[i,c(1:3)] <- object@Results@Identification[which(object@Results@Identification$AlignID==i),c("AlignID","Name","tmean")]
	# }
	
	align.List[,"AlignID"] <- N.groups
	align.List[,"tmean"] <- round(rt.mean, digits=4)
	align.List[,"Factor"] <- apply(align.List,1,function(x) paste("Factor #", x["AlignID"] , sep="")) 
	align.List[,"FoundIn"] <- foundIn.vector
	align.List[,"Spectra"] <- apply(as.matrix(1:ncol(spectra.matrix)),1,function(j){
			spectra <- spectra.matrix[,j]
			spectra.index <- which(spectra!=0) 
			spectra.pos <- spectra.index
			spectra.int <- round(spectra[spectra.index]*1000)
			spectra.text <- paste(sweep(as.matrix(spectra.pos),1,as.matrix(spectra.int),"paste.sp"), collapse=" ")
			spectra.text
		})	
	
	align.List <- as.data.frame(align.List[order(as.vector(as.numeric(align.List[,"tmean"]))),], row.names=1:nrow(align.List)) 
	align.List[,"tmean"] <- as.numeric(as.vector(align.List[,"tmean"]))
	align.List[,"FoundIn"] <- as.numeric(as.vector(align.List[,"FoundIn"]))
	align.List[,"AlignID"] <- as.integer(as.vector(align.List[,"AlignID"]))
	align.List[,"Spectra"] <- as.character(align.List[,"Spectra"])
	align.List[,"Factor"] <- as.character(align.List[,"Factor"])

	for(i in 6:ncol(align.List)) align.List[,i] <- as.numeric(as.vector(align.List[,i]))
	
	align.List
}









