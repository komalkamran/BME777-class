
library(rmr2)
library(rhdfs)

setwd( ) # initialize with the path to your working directory.

databme777<-read.csv() # study the read.csv command parameters and read the dataset csv file.

#head(databme777)

hdfs.init()
databme777.values <- to.dfs(databme777)

# write your own map() and reduce() functions based on your assigned queries.

databme777.map.fn <- function(k,v) {
keyval( k, v)
}
databme777.reduce.fn <- function(k,v) {
keyval( k,v)
}

# study mapreduce function and pass appropriate inputs and ouputs.

dataex <- mapreduce(input= ,
                   map = ,
                   reduce = )

new_var<-from.dfs(dataex)

# write appropriate code to format the data matrix you want to write to a csv file.

#y<-matrix(k,nrow=l,ncol=2,byrow=TRUE)

#write.csv(y,'BME777_BigData_Extract.csv')
