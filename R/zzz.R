.onLoad <- function(libname, pkgname) {
    library.dynam("rarpack", pkgname, libname);
}

.onUnload <- function(libpath) {
    library.dynam.unload("rarpack", libpath);
}
