subSetLib <- function(database, indexes)
{
	database.sub <- database
	database.sub@database <- database@database[indexes]
	database.sub	
}

