install.packages("devtools")
load("EmpiricalStudy/SessionInfo.RData")
for (i in 1:length(pckglist)) {
  devtools::install_version(package=pckglist[[i]]$Package, 
                            version = pckglist[[i]]$Version)
}
if (grepl("linux", sessionInfo()$platform)) install.packages("snow")
install.packages("dynConfiR_0.0.3.tar.gz", repos = NULL, type = "source")