severity<-read.table("cardiovascular_severe.dat",header=T)
mortality<-read.table("cardiovascular.dat",header=T)
RPP_filtered<-read.table("rpp_73.dat",header = T)
usethis::use_data(severity,overwrite = T)
usethis::use_data(mortality,overwrite = T)
usethis::use_data(RPP_filtered,overwrite = T)
