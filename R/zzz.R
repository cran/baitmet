
.onAttach <- function(libname, pkgname) {
	        
	data <- try(XML::xmlParse("http://cran.r-project.org/package=baitmet", isHTML=T), silent=T)
	currVersion <- as.character(utils::packageVersion("baitmet"))
	
	msg <- paste("Welcome to Baitmet. This is an early release of Baitmet (V",currVersion,"). For bugs, problems and issues, please do not hesitate in contacting xdomingo@scripps.edu and describing your problem.", sep="") 
    #packageStartupMessage(msg)   
	
	if(class(data)[1]!="try-error")
	{
		xml_data <- XML::xmlToList(data)
		cranVersion <- xml_data$body$table[[1]][[3]]
		#vers_state <- (cranVersion!=currVersion)
		vers_state <- utils::compareVersion(cranVersion, currVersion)
		if(vers_state==-1) vers_state <- 0
		wMsg <- paste("The current version of Baitmet (", currVersion, ") is outdated. There is a new version of Baitmet (", cranVersion, "), available at CRAN. To update your version execute the follwing: \n install.packages('baitmet')", sep="")
		if(vers_state) warning(wMsg)
	}else{
		wMsg <- paste("The current available version of Baitmet in CRAN cannot be checked. Please, make sure that your current version of Baitmet (", currVersion, ") is the same as in CRAN (http://cran.r-project.org/package=baitmet). To update your version execute the follwing: \n install.packages('baitmet')", sep="")
		warning(wMsg)
	}       	        
	 	
}    
