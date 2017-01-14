identifyComp_baitmet <- function(Experiment, id.database)
{
	
	#	Experiment <- exta
	#	id.database <- golm.kegg
		####
		
		compare.only.mz <- Experiment@Results@Parameters@Alignment$mz.range

		id.par <- list(database.name = id.database@name, compare.only.mz = compare.only.mz, n.putative=1)
		Experiment@Results@Parameters@Identification <- id.par
	
	avoid.processing.mz <- Experiment@Data@Parameters$avoid.processing.mz
	maxMZ <- max(compare.only.mz)
	Experiment@Results@Identification <- idenBaitmet(Experiment, maxMZ, compare.only.mz, avoid.processing.mz, id.database)
	Experiment

}

idenBaitmet <- function(Experiment, maxMZ, compare.only.mz, avoid.processing.mz, id.database)
{	

	
	factors.list <- Experiment@Data@FactorList
	empty.samples <- which(lapply(factors.list,nrow)==0)
	if(length(empty.samples)!=0) factors.list <- factors.list[-empty.samples]
		
	
	spectra.msp <- Experiment@Results@Alignment$Spectra	
	spectra.list <- lapply(spectra.msp,function(x) erah:::convertMSPspectra(x,maxMZ))	
	spectra.mat <- do.call(cbind, spectra.list)
	
	time.vector <- round(as.numeric(as.vector(Experiment@Results@Alignment$tmean)),3)
	alignID.vector <- as.numeric(as.vector(Experiment@Results@Alignment$AlignID))
	foundIn.vector <- as.numeric(as.vector(Experiment@Results@Alignment$FoundIn))
	identification.metadata <- cbind(matrix(round(time.vector, digits=4),ncol=1), matrix(alignID.vector, ncol=1), matrix(foundIn.vector,ncol=1))
		
	vect.j <- vector()
	#for(j in 1:max(as.numeric(as.vector(Experiment@Results@Alignment$AlignID))))
	for(j in alignID.vector)	
	{
		for(i in 1:length(factors.list))
		{				
			ID.loc <- which(as.numeric(as.vector(factors.list[[i]]$AlignID))==alignID.vector[j])
			if(length(ID.loc)!=0) vect.j[j] <- as.numeric(as.vector(factors.list[[i]]$DB.Id[ID.loc]))
		}		
	}
	
	id.name <- unlist(lapply(vect.j, function(x) id.database@database[[x]]$Name))
	id.form <- unlist(lapply(vect.j, function(x) id.database@database[[x]]$Formula))
	id.kegg <- unlist(lapply(vect.j, function(x) id.database@database[[x]]$KEGG))
	id.cas <- unlist(lapply(vect.j, function(x) id.database@database[[x]]$CAS))

	ref.list <- lapply(vect.j,function(x) id.database@database[[x]]$Spectra)	
	ref.list.s <- lapply(ref.list, function(x) erah:::convertMSPspectra.dot(x,maxMZ))
	ref.mat <- do.call(cbind, ref.list.s)
	ref.mat[Experiment@Data@Parameters$avoid.processing.mz,] <- 0
	
	id.match <- cor(spectra.mat,ref.mat)
	
	identification.results <- cbind(id.name, round(diag(id.match)*100,2), vect.j, id.form, id.kegg, id.cas)
	colnames(identification.results) <- c("Name.1","MatchFactor.1","DB.Id.1","Formula.1","KEGG.1","CAS.1")

	identification.list <- cbind(spectra.msp,identification.metadata,identification.results)
	colnames(identification.list) <- c("Spectra","tmean","AlignID","FoundIn",colnames(identification.results))
	identification.list <- identification.list[,c(1,3,2,5,4,6,7,8,9,10)]
	
	#cat("Done! \n")
	id.list <- as.data.frame(identification.list, row.names=1:nrow(identification.list))
	id.list[,"Spectra"] <- as.character(id.list[,"Spectra"])
	id.list
}



