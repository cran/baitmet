
processSample.baitmet <- function(Experiment, index, Chrm.Method, Rindexes, RefMat)
{
	if(Experiment@MetaData@DataDirectory=="") {filename <- as.character(Experiment@MetaData@Instrumental$filename[index])
		}else{filename <- paste(Experiment@MetaData@DataDirectory,"/",Experiment@MetaData@Instrumental$filename[index], sep="")}
	
	sampleObject <- NULL
	sampleObject <- erah:::load.file(filename)
			
	Experiment@Data@Parameters$scans.per.second <- sampleObject@scans.per.second
	sampleObject@avoid.processing.mz <- Experiment@Data@Parameters$avoid.processing.mz
	sampleObject@min.peak.width <- Experiment@Data@Parameters$min.peak.width*Experiment@Data@Parameters$scans.per.second*60
	sampleObject@min.peak.height <- Experiment@Data@Parameters$min.peak.height
	sampleObject@noise.threshold <- Experiment@Data@Parameters$noise.threshold
	sampleObject@compression.coef <- Experiment@Data@Parameters$compression.coef

	Ref.Loc.mat <- RefMat
	Ref.Loc.mat[sampleObject@avoid.processing.mz,] <- 0				
	Ref.Loc.mat <- Ref.Loc.mat[sampleObject@min.mz:sampleObject@max.mz,]
	
	#Fer tambe per la REFMAT, treure les masses que no li toquin.
	# if(min(sampleObject@avoid.processing.mz)<ssampleObject@min.mz) stop(paste("Error: Minimum Mz adquired is higher than Avoid Processing Mz. Please change Avoid Processing Mz parameter in SetAlPar() function, whith at least from Mz number",sampleRD@min.mz))
	# avoid.mz <- sampleRD@avoid.processing.mz - (sampleRD@min.mz - 1)
	# sampleRD@data[,avoid.mz] <- 0
	# sampleRD


	# Pre-processing:
	sampleObject <- erah:::avoid.processing(sampleObject)
	k.filt <- round(sampleObject@min.peak.width)
	sampleObject@data <- erah:::pre.process(sampleObject@data, k.filt)
	#sampleObject@data <- apply(sampleObject@data,2,function(x) erah:::removeBaseline(x,k.filt*10))
	#if(erah:::is.even(k.filt)) k.filt <- k.filt + 1
	#sampleObject@data <- erah:::soft.filter(sampleObject@data,3,k.filt)
		
	factor.list <- try(process.samples.target(sampleObject, Chrm.Method, Rindexes, Ref.Loc.mat, Experiment@Data@Parameters$ri.error, round(sampleObject@min.peak.width*4)), silent=F)
	if(class(factor.list)=="try-error") {factor.list <- as.data.frame(NULL); warning("Unable to extract factors from ", Experiment@MetaData@Instrumental$filename[index], ". Data may be corrupted.", sep="")}
	
	Experiment@Data@FactorList[[index]] <- factor.list		
	Experiment
}



process.samples.target <- function(sampleRD, Chrm.Method, RI.vect, REF.mat, ri.error, span.len.local, DB.id.vect=NULL)	
{	
	DB.id.vect=NULL
	
	RI.curve <- smooth.spline(Chrm.Method$ref.ri, Chrm.Method$ref.rt, df=length(Chrm.Method$ref.ri), tol=1)
		
	Chrm.Data <- sampleRD@data		
	C.matrix <- matrix(0,nrow=nrow(Chrm.Data), ncol=0)
	S.matrix <- matrix(0,nrow=ncol(Chrm.Data), ncol=0)
	ID.vector <- vector()
	MF.vector <- vector()
	
	if(is.null(DB.id.vect)) DB.id.vect <- 1:length(RI.vect)
		
	## Sigma Scans:
	
	sigma.scans <-  sampleRD@min.peak.width
	ksp <- vector()
	sewq <- seq(0.1,1000,0.01)
	for( i in 1:length(sewq))
	{
		krn <- normalize(dnorm(1:nrow(sampleRD@data),trunc(nrow(sampleRD@data)/2),sewq[i]))
		k.span <- trunc((length(which(normalize(krn)>0.01))/2)+0.5)	 	
		ksp[i] <- k.span
		if(k.span - sigma.scans >0) break
	}
	sigma.model <- sewq[(i-1)]		
		
	pb <- txtProgressBar(min=1,max=ncol(REF.mat), width=50, style=3)
	for(Cmp in 1:ncol(REF.mat))
	{
		setTxtProgressBar(pb, Cmp)
		Cmp.RI <- RI.vect[Cmp]
		Put.RT <- predict(RI.curve,Cmp.RI)$y
		Target.S <- REF.mat[,Cmp]
		
		if(Put.RT<sampleRD@start.time/60) next
		if(Cmp.RI==0) next
		
		Low.scan <- trunc((predict(RI.curve,Cmp.RI-Cmp.RI*ri.error)$y - sampleRD@start.time/60)*sampleRD@scans.per.second*60)
		Upp.scan <- trunc((predict(RI.curve,Cmp.RI+Cmp.RI*ri.error)$y - sampleRD@start.time/60)*sampleRD@scans.per.second*60)
		window.span <- Low.scan:Upp.scan
		if(min(window.span)<1) window.span <- window.span[-which(window.span<1)]
		if(max(window.span)>nrow(Chrm.Data)) window.span <- window.span[-which(window.span>ncol(Chrm.Data))]
		Cmp.Matrix <- Chrm.Data[window.span,] #* hamming(length(window.span))
						
		sim.vector <- suppressWarnings(cor(t(Cmp.Matrix),Target.S))
		sim.vector[!is.finite(sim.vector)] <- 0
		MaxVal <- suppressWarnings(max(sim.vector, na.rm=T))
		if(!is.finite(MaxVal)) next
		if(MaxVal<=0.5) next
		
		Put.scans <- window.span[1] + which.max(sim.vector)	- 1		
		window.span <- trunc((Put.scans-span.len.local):(Put.scans+span.len.local))
		if(min(window.span)<1) window.span <- window.span[-which(window.span<1)]
		if(max(window.span)>nrow(Chrm.Data)) window.span <- window.span[-which(window.span>ncol(Chrm.Data))]
		Cmp.Matrix <- Chrm.Data[window.span,]
	
			## Profile Deconvolution:			
			C.mod <- try(erah:::getC.rq(Cmp.Matrix*hanning(nrow(Cmp.Matrix)), Target.S), silent=T)
			if(class(C.mod)=="try-error") next
			#C.mod <- try(erah:::getC.tP(Cmp.Matrix*hanning(nrow(Cmp.Matrix)), Target.S), silent=T)			
			#C.mod <- apply(Cmp.Matrix*hanning(nrow(Cmp.Matrix)), 1, function(x) as.numeric(nnls(as.matrix(normalize(Target.S)),x)$x))
			if(sum(C.mod)==0) next

			CMod.OSD <- dnorm(1:length(C.mod), which.max(C.mod), sigma.model)			

			S.mod <- erah:::getS.OSD(CMod.OSD,Cmp.Matrix)  #, ref.response=Target.S)			
			if(sum(S.mod)==0) next
					
			## _New inserction
			#C.mod <- try(erah:::getC.rq(Cmp.Matrix, S.mod), silent=T)
			#if(class(C.mod)=="try-error") next
			
			if(which.max(C.mod)==1) next
			if(which.max(C.mod)==length(C.mod)) next
			
			#C.mod <- chrom.isoreg(C.mod)	
			MF.score <- suppressWarnings(cor(Target.S, S.mod))
			
		ID.vector <- c(ID.vector,DB.id.vect[Cmp])	
		MF.vector <- c(MF.vector, MF.score)
		S.matrix <- cbind(S.matrix,as.matrix(S.mod))
		C.matrix.aux <- matrix(0,nrow=nrow(Chrm.Data), ncol=1)
		C.matrix.aux[window.span,] <- C.mod
		C.matrix <- cbind(C.matrix, C.matrix.aux)
	}

	## START: Per evitar overlap:
		# CorGen <- suppressWarnings(fastCor(C.matrix))
		
		# C.out.inds <- vector()
		# is.with <- list()
		# for(j in 1:ncol(C.matrix))
		# {				
			# is.with[[j]] <- ID.vector[j]
			# if(j %in% C.out.inds) next
			# #if(j %in% C.ok.inds) next
			# conf <- which(CorGen[,j]>0.95)
			# if(length(conf)==1) next
			
			# delete.i <- conf[-which.max(MF.vector[conf])]
			# C.out.inds <- c(C.out.inds,delete.i)
			
			# is.with[[j]] <- c(is.with[[j]],ID.vector[delete.i])
		# }
		# C.out.inds <- unique(C.out.inds)
		#window.features <- list(id.vector=is.with[-C.out.inds], profile=C.matrix[,-C.out.inds], spectra=S.matrix[,-C.out.inds])

	## FI: Per evitar overlap:

	
	#matplot(sampleRD@data, type="l", lty=1, col="gray")
	#matplot(C.matrix[-(1:299),-C.out.inds], type="l", lty=1, add=T)
	
	## CREATING TABLE:
		
	window.features <- list(id.vector=ID.vector, profile=C.matrix, spectra=S.matrix)
	
	## Peak Position
	iter.ind <- seq(1,ncol(window.features$profile))
	models.maximas <- unlist(apply(window.features$profile,2,which.max))
	#models.maximas <- models.maximas + (trunc(from.s) - 1)
	
	## Quantification
	#models.areas <- sapply(iter.ind,function(j){ sum(as.matrix(window.features$profile)[,j]%*%t(as.matrix(window.features$spectra)[,j])) })
	
	colSum.C <- colSums(window.features$profile)
	models.areas <- sapply(iter.ind,function(j){ sum(colSum.C[j]*t(as.matrix(window.features$spectra)[,j])) })	
		
	## Peak Height
	models.peakheight <-  apply(as.matrix(window.features$profile),2, function(x) max(x, na.rm=T))
		
	## Spectra:
	models.spectra <- apply(as.matrix(window.features$spectra),2,function(spectra){
			spectra.index <- which(normalize(spectra)>0.005) 
			spectra.pos <- spectra.index + (sampleRD@min.mz - 1)
			spectra.int <- round(spectra[spectra.index]*1000)
			spectra.text <- paste(sweep(as.matrix(spectra.pos),1,as.matrix(spectra.int),"paste.sp"), collapse=" ")
			spectra.text
		})	

	## Profile:
	models.profile <- apply(as.matrix(window.features$profile),2,function(profile){
			profile.index <- which(normalize(profile)>0.0001) 
			profile.pos <- ((profile.index)/sampleRD@scans.per.second/60) + (sampleRD@start.time/60)
			profile.int <- profile[profile.index]/max(profile)
			profile.text <- paste(sweep(as.matrix(profile.pos),1,as.matrix(profile.int),"paste.sp"), collapse=" ")
			profile.text
		})	
	
	## Cmp.DB:
	
	id.db <- sapply(window.features$id.vector,function(dbID){ paste(dbID, collapse=",") })	
	
	feature.list <- data.frame()
	
	feature.prelist <- apply(as.matrix(iter.ind),1,function(window.number){
		items <- length(models.maximas[window.number])
		subfeat.list <- matrix(0,ncol=7,nrow=items)		
		subfeat.list[,2] <- models.maximas[window.number]
		subfeat.list[,3] <- id.db[window.number]
		subfeat.list[,4] <- round(models.areas[window.number], digits=0)
		subfeat.list[,5] <- round(models.peakheight[window.number], digits=0)
		subfeat.list[,6] <- models.spectra[window.number]
		subfeat.list[,7] <- models.profile[window.number]
		subfeat.list
		})

	feature.list <- as.data.frame(t(feature.prelist))
	feature.list[,1] <- seq(1,nrow(feature.list))
	feature.list[,2] <- round((as.numeric(as.vector((feature.list[,2])))/sampleRD@scans.per.second)/60 + (sampleRD@start.time/60), digits=4)	
	colnames(feature.list) <- c("ID","RT","DB.Id","Area","Peak Height","Spectra","Profile")	

	feature.list
	
}

