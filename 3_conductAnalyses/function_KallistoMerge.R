countmerge = function(mypath){
    filenames=list.files(path=mypath,pattern = "*.tsv", recursive=T, full.names=TRUE)
    datalist = lapply(filenames, function(x){read.table(file=x,header=T,colClasses=c(NA,"NULL","NULL",NA,"NULL"))})
    Reduce(function(x,y) {merge(x,y,by="target_id")}, datalist) 
}



tpmmerge = function(mypath){
    filenames=list.files(path=mypath,pattern = "*.tsv", recursive=T, full.names=TRUE)
    datalist = lapply(filenames, function(x){read.table(file=x,header=T,colClasses=c(NA,"NULL","NULL","NULL",NA))})
    Reduce(function(x,y) {merge(x,y,by="target_id")}, datalist) 
}