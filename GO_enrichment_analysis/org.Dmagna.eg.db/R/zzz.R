datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Dmagna.eg <- function() showQCData("org.Dmagna.eg", datacache)
org.Dmagna.eg_dbconn <- function() dbconn(datacache)
org.Dmagna.eg_dbfile <- function() dbfile(datacache)
org.Dmagna.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Dmagna.eg_dbInfo <- function() dbInfo(datacache)

org.Dmagna.egORGANISM <- "Daphnia magna"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Dmagna.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Dmagna.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Dmagna.eg_dbconn())
}

