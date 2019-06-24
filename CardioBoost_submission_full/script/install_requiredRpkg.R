require(devtools)
pkg_list=read.table("./requiredRpkgs.txt",stringsAsFactors=FALSE,sep="_")
for (i in 1:nrow(pkg_list)){
install_version(package = pkg_list[i,1],version=pkg_list[i,2],repos = "http://cran.us.r-project.org")
}
