
library(rmr2)
library(rhdfs)

#setwd('/home/student1/k2kamran/Desktop/BME77/Project') # initialize with the path to your working directory.

databme777<-read.csv("new_diabetic_data_V2.csv") # study the read.csv command parameters and read the dataset csv file.



hdfs.init()
databme777.values <- to.dfs(databme777)

# write your own map() and reduce() functions based on your assigned queries.


d <- c(460:519, 786) #range of values

databme777.map.fn <- function(k,v) {
  p <- which((v[,4] == 1) & (v[,8] %in% d)) #Query: 
keyval(v[p,], v[p,c(4,8)])
}
databme777.reduce.fn <- function(k,v) {
keyval( k,(unlist(v)))
}

d <- c(460:519, 786) #range of values

databme777.map.fns <- function(k,v) {
  p <- which((v[,4] != 1) & (v[,8] %in% d)) #Query: 
keyval(v[p,], v[p,c(4,8)])
}
databme777.reduce.fns <- function(k,v) {
keyval( k,(unlist(v)))
}
#777.map.fns <- function(k,v) {
 #p <- which((v[,4] != 1) & (v[,8] %in% d))#Query: 
#keyval(v[p,], v[p,c(4,8)])

#}
#databme777.reduce.fns <- function(k,v) {
#keyval( k,(unlist(v)))







#}

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

#k<- unlist(new_var[2])
#l <- length(k)/2

#o <- unlist(new_vars[2])
#ll <- length(o)/2


# write appropriate code to format the data matrix you want to write to a csv file.
#k <- unlist(new_var)
#l <- length(k)/2

#y<-matrix(k,nrow=l,ncol=10,byrow=TRUE)
#q<-matrix(o,nrow=ll,ncol=10,byrow=TRUE)

write.csv(j,'BME777_BigData_Extract4.csv', row.names = FALSE)
#write.csv(q,'BME777_BigData_Extract5.csv', row.names = FALSE)


