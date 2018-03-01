
library(rmr2)
library(rhdfs)

#setwd('/home/student1/m2johar/Desktop/bme777') # initialize with the path to your working directory.

databme777<-read.csv("new_diabetic_data_V2.csv") # study the read.csv command parameters and read the dataset csv file.

#which(complete.cases(databme777))
#databme777[!complete.cases(databme777),]
#head(databme777)

hdfs.init()
databme777.values <- to.dfs(databme777)

# write your own map() and reduce() functions based on your assigned queries.
d <- 800:999
databme777.map.fn <- function(k,v) {
  p <- which((as.numeric(v[,2]) <= 7) & v[,8] %in% d)
keyval(v[p,], v[p,c(4,8)])
}
databme777.reduce.fn <- function(k,v) {
keyval( k,(unlist(v)))
}


databme777.map.fns <- function(k,v) {
  p <- which((as.numeric(v[,2]) > 7) & v[,8] %in% d)
keyval(v[p,], v[p,c(4,8)])
}
databme777.reduce.fns <- function(k,v) {
keyval( k,(unlist(v)))
}

# study mapreduce function and pass appropriate inputs and ouputs.

dataex <- mapreduce(input=databme777.values ,
                   map =databme777.map.fn ,
                   reduce = databme777.reduce.fn)

dataexs <- mapreduce(input=databme777.values ,
                   map =databme777.map.fns ,
                   reduce = databme777.reduce.fns)



new_var<-from.dfs(dataex)
new_vars<-from.dfs(dataexs)

k <- rbind(as.data.frame(new_var[1]), as.data.frame(new_vars[1]))
j <- as.data.frame(k[!duplicated(as.data.frame(k)),])

# write appropriate code to format the data matrix you want to write to a csv file.
#k <- unlist(new_var)
#l <- length(k)/2
#y<-matrix(k,nrow=l,ncol=2,byrow=TRUE)

write.csv(j,'BME777_BigData_Extract1.csv', row.names = FALSE)
