mzList <- function(Experiment, by.area = TRUE)
{
	QGenTab <- 0
	RGenTab <- 0
	for(i in 1:length(Experiment@Data@FeatureList))
	{
		#AlV <- paste(Experiment@Data@FeatureList[[i]]$AlignID, Experiment@Data@FeatureList[[i]]$mz, sep=".")
		AlVG <- t(apply(Experiment@Data@FeatureList[[i]][,c("AlignID","mz")], 1, function(x) c(as.vector(as.matrix(x[1])),paste(as.matrix(x), collapse="."))))
		colnames(AlVG) <- c("AlignID", "AlV")
		
		if(by.area==TRUE) QTab <- Experiment@Data@FeatureList[[i]][,c(which(colnames(Experiment@Data@FeatureList[[i]])=="AlignID"),which(colnames(Experiment@Data@FeatureList[[i]])=="mz"),which(colnames(Experiment@Data@FeatureList[[i]])=="RT"),which(colnames(Experiment@Data@FeatureList[[i]])=="Area"))]
		if(by.area==FALSE) QTab <- Experiment@Data@FeatureList[[i]][,c(which(colnames(Experiment@Data@FeatureList[[i]])=="AlignID"),which(colnames(Experiment@Data@FeatureList[[i]])=="mz"),which(colnames(Experiment@Data@FeatureList[[i]])=="RT"),which(colnames(Experiment@Data@FeatureList[[i]])=="Int"))]

		LocTab <- as.data.frame(t(apply(QTab, 1, function(x){ as.vector(c(paste(as.matrix(x[1:2]), collapse="."),as.matrix(x[3]),as.matrix(x[4])))})))
		if(by.area==TRUE) colnames(LocTab) <- c("AlV","RT","Area")
		if(by.area==FALSE) colnames(LocTab) <- c("AlV","RT","Int")
		
		if(by.area==TRUE) QGenTab <- LocTab[,c("AlV","Area")]
		if(by.area==FALSE) QGenTab <- LocTab[,c("AlV","Int")]
		RGenTab <- LocTab[,c("AlV","RT")]

		if(i>1){QGen <- merge(QGen, QGenTab, by="AlV")
			colnames(QGen) <- c("AlV", names(Experiment@Data@FactorList)[1:i])
			RGen <- merge(RGen, RGenTab, by="AlV")
			colnames(RGen) <- c("AlV", names(Experiment@Data@FactorList)[1:i])	
			}
		if(i==1) {QGen <- QGenTab
			RGen <- RGenTab
			}
	}	
	
	AlV <- as.vector(QGen[,"AlV"])
	RTmean <- round(apply(RGen[,-1],1,function(x) mean(as.numeric(as.vector(x)))),3)
	MetaTab <- as.matrix(t(sapply(AlV, function(x) as.numeric(as.vector(strsplit(x, "\\.")[[1]])))))
	rownames(MetaTab) <- NULL
	FinMZtab <- (cbind(MetaTab,RTmean,QGen[,-1]))
	colnames(FinMZtab)[1:3] <- c("AlignID", "MZ", "RT")		
	FinMZtab 
}

quantSM <- function(Experiment, ms.library, AlignID=NULL, pre.process=TRUE, fit.gaussian=FALSE)
{
	Number.of.Samples <- nrow(Experiment@MetaData@Instrumental)	
	k <- 1
	FeatList <- list()
	for(index in 1:Number.of.Samples)
	{
		cat("\n Quantifying selective masses in",as.character(Experiment@MetaData@Instrumental$filename[index]),"... Processing", k,"/",Number.of.Samples,"\n")  
		FeatList[[k]] <- quantifySelMasses(Experiment, index, ms.library, AlignID=AlignID, pre.process=pre.process, fit.gaussian=fit.gaussian)
		k <- k + 1
	}
	names(FeatList) <- names(Experiment@Data@FactorList) 
	Experiment@Data@FeatureList <- FeatList
	Experiment
}


quantifySelMasses <- function(Experiment, index, ms.library, AlignID, pre.process, fit.gaussian)
{
	if(Experiment@MetaData@DataDirectory=="") {filename <- as.character(Experiment@MetaData@Instrumental$filename[index])
		}else{filename <- paste(Experiment@MetaData@DataDirectory,"/",Experiment@MetaData@Instrumental$filename[index], sep="")}
		
	sampleObject <- NULL
	sampleObject <- erah:::load.file(filename)
			
	#Experiment@Data@Parameters$scans.per.second <- sampleObject@scans.per.second
	sampleObject@avoid.processing.mz <- Experiment@Data@Parameters$avoid.processing.mz
	sampleObject@min.peak.width <- Experiment@Data@Parameters$min.peak.width*Experiment@Data@Parameters$scans.per.second*60
	sampleObject@min.peak.height <- Experiment@Data@Parameters$min.peak.height
	sampleObject@noise.threshold <- Experiment@Data@Parameters$noise.threshold
	sampleObject@compression.coef <- Experiment@Data@Parameters$compression.coef

	#Ref.Loc.mat <- RefMat
	#Ref.Loc.mat[sampleObject@avoid.processing.mz,] <- 0				
	#Ref.Loc.mat <- Ref.Loc.mat[sampleObject@min.mz:sampleObject@max.mz,]
	if(is.null(AlignID)) AlignID <- unique(as.numeric(as.vector(Experiment@Results@Identification$AlignID)))
	
	if(pre.process) 
	{
		sampleObject <- erah:::avoid.processing(sampleObject)
		k.filt <- round(sampleObject@min.peak.width)
		sampleObject@data <- erah:::pre.process(sampleObject@data, k.filt)
	}
	
	## FeatureList
	
	FeatList <- matrix(0, ncol=6, nrow=0)
	FeatList <- as.data.frame(FeatList)
	
	## Quantification:
	id.List <- Experiment@Results@Identification
	SpanLocal <- round(sampleObject@min.peak.width*3)
	AlignID.vect <- as.numeric(as.vector(Experiment@Results@Identification$AlignID))
	maxMZ <- max(Experiment@Results@Parameters@Identification$compare.only.mz)
	for(Cmp.AlID in AlignID)
	{
		ListI <- which(AlignID.vect==Cmp.AlID)
		if(length(ListI)==0) stop("The following AlignID number: ", Cmp.AlID, " does not exist in the experiment, please, make sure that all the AlignID selected exist by executing idList(Experiment)")
		LocalCMP.DBId <- as.numeric(as.vector(id.List[ListI,]$DB.Id.1))
		
		FL.LocInd <- which(as.numeric(as.vector(Experiment@Data@FactorList[[index]]$AlignID))==Cmp.AlID)
		LocalCMP.RT <- as.numeric(as.vector(id.List[ListI,]$tmean))
		if(length(FL.LocInd)!=0) LocalCMP.RT <- as.numeric(as.vector(Experiment@Data@FactorList[[index]]$RT))[FL.LocInd]
	
		if(is.null(ms.library@database[[LocalCMP.DBId]]$SelMZ))
		{
			LocalSpect <- erah:::convertMSPspectra(id.List[ListI,]$Spectra, maxMZ)
			selMass <- order(LocalSpect, decreasing=T)[1:5]
			selMass[LocalSpect[selMass]==0] <- 0
			
		}else{
			if(all(is.na(ms.library@database[[LocalCMP.DBId]]$SelMZ)))
			{
				LocalSpect <- erah:::convertMSPspectra(id.List[ListI,]$Spectra, maxMZ)
				selMass <- order(LocalSpect, decreasing=T)[1:5]
				selMass[LocalSpect[selMass]==0] <- 0
			}else{
				selMass <- ms.library@database[[LocalCMP.DBId]]$SelMZ
			}
		}
		selMass <- selMass[-selMass<=0]

		LocalScan <- trunc((LocalCMP.RT-sampleObject@start.time/60)*sampleObject@scans.per.second*60)  
		
		
		Low.scan <- LocalScan -	SpanLocal
		Upp.scan <- LocalScan +	SpanLocal
		window.span <- Low.scan:Upp.scan
		if(min(window.span)<1) window.span <- window.span[-which(window.span<1)]
		if(max(window.span)>nrow(sampleObject@data)) window.span <- window.span[-which(window.span>ncol(sampleObject@data))]
		Low.scan <- min(window.span)
		Cmp.Matrix <- sampleObject@data[window.span,] #* hamming(length(window.span))
		
		selMass.Q <- c(selMass - (sampleObject@min.mz - 1))
		
		if(fit.gaussian==TRUE)
		{
			Cmp.MatrixFitted <- apply(Cmp.Matrix[,selMass.Q], 2, fitGaussian)
			SelAreas <- trunc(colSums(Cmp.MatrixFitted))
			SelInt <- trunc(apply(as.matrix(Cmp.MatrixFitted),2,max))
		}else{
			SelAreas <- trunc(colSums(Cmp.Matrix[,selMass.Q]))
			SelInt <- trunc(apply(as.matrix(Cmp.Matrix[,selMass.Q]),2,max))
		}
		models.profile <- apply(as.matrix(Cmp.Matrix[,selMass.Q]),2,function(profile){
			profile.index <- which(normalize(profile)>0.0001) 
			profile.pos <- round(((profile.index + (Low.scan - 1))/sampleObject@scans.per.second/60) + (sampleObject@start.time/60),5)
			profile.int <- profile[profile.index]/max(profile)
			profile.text <- paste(sweep(as.matrix(profile.pos),1,as.matrix(profile.int),"paste.sp"), collapse=" ")
			profile.text
		})	
		
		Mass.Scan <- apply(as.matrix(Cmp.Matrix[,selMass.Q]),2,which.max) + min(window.span)
		Mass.RT <- round(Mass.Scan/sampleObject@scans.per.second/60 + sampleObject@start.time/60,4)
		Mass.AI <- rep(Cmp.AlID, length(selMass))
		FeatList <- rbind(FeatList,cbind(Mass.AI, selMass, Mass.RT, SelAreas, SelInt, models.profile)) 	
	}
	colnames(FeatList) <- c("AlignID", "mz", "RT", "Area", "Int", "Profile")
	FeatList
}