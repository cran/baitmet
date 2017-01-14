computeRI <- function(Experiment, ms.library=mslib, IS.alignid=NULL)
{	

	id.database <- ms.library
	
	if(!is.null(IS.alignid))
	{
		
		factors.list <- Experiment@Data@FactorList
		empty.samples <- which(lapply(factors.list,nrow)==0)
		factors.list.original <- NULL
		if(length(empty.samples)!=0)
		{
			factors.list.original <- factors.list
			factors.list <- factors.list[-empty.samples] 
		}
		if(length(factors.list)==1) stop("Only one sample has been processed. No alignment needed")
		N.samples <- length(factors.list)
		
		IS.genId <- which(Experiment@Results@Identification$AlignID %in% IS.alignid)
		IS.ri <- sapply(as.numeric(as.vector(Experiment@Results@Identification$DB.Id[IS.genId])), function(i) id.database@database[[i]]$RI.VAR5.ALK)
		
		#load("FAMErt.rda")
		#dim(FAME.rt)
		#err.v <- FAME.rt*0
		
		for(Smp in 1:N.samples)
		{
		
			IS.inds <- as.numeric(as.vector(sapply(IS.alignid, function(y) which(y==factors.list[[Smp]]$AlignID))))
			
			if(any(is.na(IS.inds)))
			{
				IS.ri <- IS.ri[-which(is.na(IS.inds))]
				IS.inds <-  IS.inds[-which(is.na(IS.inds))]
			}
				
			IS.rt <- factors.list[[Smp]]$RT[IS.inds]
			
			#err.v[,Smp] <- IS.rt - FAME.rt[,Smp]
			
			chrom.method <- setChrmMethod(method="alk.var5", rt=IS.rt, ri=IS.ri, name="IS Std")
			RtIS <- smooth.spline(chrom.method$ref.rt , chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1))
			rt.loc.o <- as.numeric(as.vector(factors.list[[Smp]]$RT))
			emp.ri.new <- predict(RtIS, rt.loc.o)$y
			factors.list[[Smp]]$Emp.ri <- emp.ri.new		
		}
		#boxplot(t(abs(err.v)), main="RT error")
	
	}
	
	if(is.null(IS.alignid))
	{
			chrom.method <- Experiment@Results@Parameters@Alignment$chrom.method

			DataDBID <- as.numeric(as.vector(Experiment@Results@Identification$DB.Id.1))
			RefRI <- sapply(DataDBID, function(x) id.database@database[[x]]$RI.VAR5.ALK)
	
			RIm <- matrix(rep(RefRI,length(RefRI)), ncol=length(RefRI))
			HDist <- 100*as.matrix(dist(RefRI), upper=T, diag=F)/t(RIm)
			RIov.vect <- unlist(apply(HDist,1,function(x) length(which(x<Experiment@Data@Parameters$ri.error*2))))
			RIov.vect <- as.vector(as.numeric(which(RIov.vect>1)))
			if(length(RIov.vect)==0) RIov.vect <- 0
	
			SPov.vect <- 0
			if(sum(RIov.vect)!=0)
			{
				RIov.SPECT <- sapply(RIov.vect,function(x) erah:::convertMSPspectra(Experiment@Results@Identification$Spectra[[x]],600))		
				CorSpect <- suppressWarnings(fastCor(RIov.SPECT))
				SPov.vect <- unlist(apply(CorSpect,1,function(x) length(which(x>0.75))))
				SPov.vect <- as.vector(as.numeric(which(SPov.vect>1)))
				if(length(SPov.vect)==0) SPov.vect <- 0
			}
			#which(RIov.vect %in% SPov.vect)
			#RIov.vect[which(RIov.vect %in% SPov.vect)]
			MFvect <- as.numeric(as.vector(Experiment@Results@Identification$MatchFactor.1))
	
			#if(sum(SPov.vect) + sum(RIov.vect)!=0)
			#{
			#	LeaveOut.vect <- as.numeric(as.vector(Experiment@Results@Identification$AlignID[RIov.vect[which(RIov.vect %in% SPov.vect)]]))
			#	MFvect[LeaveOut.vect] <- 0
			#}
			
			select.IDinds <- which(MFvect>90)
			#select.IDinds
			select.inds <- as.numeric(as.vector(Experiment@Results@Identification$AlignID[select.IDinds]))
	
			if(length(select.inds)<15) stop("Not enough number of natural standards to perform elastic statistics")
			## Fer servir una curva estatica en aquest cas!
			#############
			############		
		#}
		
	
		factors.list <- Experiment@Data@FactorList
		empty.samples <- which(lapply(factors.list,nrow)==0)
		factors.list.original <- NULL
		if(length(empty.samples)!=0)
		{
			factors.list.original <- factors.list
			factors.list <- factors.list[-empty.samples] 
		}
		if(length(factors.list)==1) stop("Only one sample has been processed. No alignment needed")
		N.samples <- length(factors.list)
		
		Rt.curve <- smooth.spline(chrom.method$ref.rt, chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1))
		
		RIerr <- matrix(0, nrow=length(select.inds), ncol=N.samples)
		
		##Per validar
		Fame.pred.rt <- matrix(0, nrow=length(chrom.method$ref.rt), ncol=N.samples)
		ElastMod <- matrix(0, nrow=2, ncol=N.samples)
		ErBai <- matrix(0, nrow=34, ncol=N.samples)
		ErFam <- matrix(0, nrow=34, ncol=N.samples)
		ErBai.Fame <- matrix(0, nrow=9, ncol=N.samples)
		
		for(Smp in 1:N.samples)
		{
			rt.loc.o <- as.numeric(as.vector(factors.list[[Smp]]$RT))
			dbid.loc.o <- as.numeric(as.vector(factors.list[[Smp]]$DB.Id))
			
			alid.loc <- as.numeric(as.vector(factors.list[[Smp]]$AlignID))
			sel.alid <- which(alid.loc %in%select.inds)
			
			if(chrom.method$method=="alk.var5") ri.emp.o <- as.numeric(as.vector(sapply(dbid.loc.o, function(x) id.database@database[[x]]$RI.VAR5.ALK)))
			if(chrom.method$method!="alk.var5") stop("Chromatographic method not supported. Please, select alk.var.5 as a chromatographic method in function setChrmMethod() ")
			
			rt.loc <- rt.loc.o[sel.alid]
			ri.emp <- ri.emp.o[sel.alid]
			
		## Per LS:
			
		#LM.Obj <- lm(rt.loc~predict(RI.curve,ri.emp)$y)
		# LinCor.ls <- as.numeric(coefficients(LM.Obj)[2])
		# LinOff.ls <- as.numeric(coefficients(LM.Obj)[1])
	
		### Per OPT:	
			LinCor <- 1
			LinOff <- 0
			error.vect <- vector()
			RI.curve <- smooth.spline(chrom.method$ref.ri, chrom.method$ref.rt, df=(length(chrom.method$ref.ri)-1))
			pred.rt <- predict(RI.curve,ri.emp)$y
			#opt.ERROR.old <- sum((pred.rt - rt.loc)^2)
			opt.ERROR.old <- median(100*(abs(ri.emp - predict(Rt.curve, rt.loc)$y))/ri.emp)
	
			local.Rt.curve <- Rt.curve
			for(i in 1:50)
			{
				if(!is.even(i))
				{
					opt.f1 <- function(yk) {
						localFunct.Rt.curve <- smooth.spline(chrom.method$ref.rt*yk+LinOff, chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1))
						RI.error <-  100*(abs(ri.emp - predict(localFunct.Rt.curve, rt.loc)$y))/ri.emp
						median(RI.error^1)
						}
										
					opt.res.pos <- optimize(opt.f1, interval=c(-200,200), tol=1e-200) 
					opt.ERROR <- opt.res.pos$objective
					if(opt.ERROR<opt.ERROR.old)
					{ 
						LinCor <- opt.res.pos$minimum
						opt.ERROR.old <- opt.ERROR
						local.Rt.curve <- smooth.spline(chrom.method$ref.rt*LinCor, chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1))
	
					}
				}else{
					opt.f1 <- function(yk) {
						localFunct.Rt.curve <- smooth.spline(chrom.method$ref.rt*LinCor+yk, chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1))
	
						RI.error <-  100*(abs(ri.emp - predict(local.Rt.curve, rt.loc)$y))/ri.emp
						median(RI.error^1)				
						}
										
					opt.res.pos <- optimize(opt.f1, interval=c(-1000,1000), tol=1e-100) 
					opt.ERROR <- opt.res.pos$objective
					if(opt.ERROR<opt.ERROR.old)
					{ 
						LinOff <- opt.res.pos$minimum
						opt.ERROR.old <- opt.ERROR
						local.Rt.curve <- smooth.spline(chrom.method$ref.rt+LinOff, chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1))
	
					}		
				}
			error.vect[i] <- opt.ERROR
			if(i>5) if(opt.ERROR>=error.vect[(i-1)]) break
			}
												
			if(is.null(factors.list[[Smp]]$Emp.RI))
			{
				factors.list[[Smp]] <- c(factors.list[[Smp]],1)	
				names(factors.list[[Smp]])[length(factors.list[[Smp]])] <- "Emp.RI"
			}
			
			RtBaiTMet <- smooth.spline((chrom.method$ref.rt*LinCor) + LinOff , chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1))
			emp.ri.new <- predict(RtBaiTMet, rt.loc.o)$y
			factors.list[[Smp]]$Emp.ri <- emp.ri.new
	
		## For Validation
		## End Validation
		}	
			
	}
	
	N.groups <- as.numeric(as.vector(Experiment@Results@Identification$AlignID))
		
		ri.align.matrix <- t(apply(as.matrix(N.groups),1,function(i){
			local.index <- as.vector(unlist(lapply(factors.list,function(x) {
				outp <- as.character(x$Emp.ri[which(x$AlignID==i)])
				if(length(outp)==0) outp <- 0
				outp
				})))
		}))
		class(ri.align.matrix) <- "numeric"
		ri.mean.adj <- rowSums(ri.align.matrix)/apply(ri.align.matrix,1,function(x) length(which(x!=0)))
		
		Gen.DBId <- as.numeric(as.vector(Experiment@Results@Identification$DB.Id.1))
		
		if(chrom.method$method=="alk.var5") ri.emp.db <- as.numeric(as.vector(sapply(Gen.DBId, function(x) id.database@database[[x]]$RI.VAR5.ALK)))
		if(chrom.method$method!="alk.var5") stop("Chromatographic method not supported. Please, select alk.var.5 as a chromatographic method in function setChrmMethod() ")
		
		#Rt.curve <- smooth.spline(chrom.method$ref.rt, chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1), tol=1)
		RI.error <-  100*(abs(ri.emp.db - ri.mean.adj))/ri.emp.db
	
		Experiment@Results@Identification$RI.error.1 <- round(RI.error,2)
	
		if(!is.null(factors.list.original))
		{
			factors.list.original[-empty.samples] <- factors.list
			factors.list <- factors.list.original
		}
		Experiment@Data@FactorList <- lapply(factors.list, function(x) {ouT <- as.data.frame(x)
			colnames(ouT)[which(colnames(ouT)=="Peak.Height")] <- "Peak Height"
			ouT
			})
	
	Experiment
}

