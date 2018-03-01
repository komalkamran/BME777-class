bme<- read.csv("diabetic_data_V2.csv", na.strings = c("?", "E909"))
z<-bme[!complete.cases(bme),]
z<- bme[which(bme$diag_1=='?'),]
x<-y[, c(2,8)]
str(y)
levels(y$age)
which.min(y$age)
dim(z)
x
mean(as.numeric(as.character(bme$diag_1)), na.rm = T)
is.numeric(as.numeric(as.character(bme$age)))
str(bme)
bme$diag_1 <- round(as.numeric(as.character(bme$diag_1)))
hist(bme$diag_1)
ggplot(data=bme, aes(x=bme$diag_1))+
  geom_histogram()
median(bme$diag_1)
z
dim(z)
x<-z[, c(2,8)]
x$age<- factor(x$age)
levels(x$age)
which(is.na(x))
a <- bme[as.numeric(bme$age)==10,]
b <-median(a$diag_1, na.rm=T)
which(as.numeric(x$age)==10)
x
bme[c(1006,49516,57058,60314,98396),8] <- b
?write.csv

new_diabetic_data_V2 <- write.csv(bme, file = "new_diabetic_data_V2.csv", row.names = F)
mean(bme$diag_1)

810 %in% d

new_data <- databme777
d <- c(800:999)
which(round(new_data$diag_1) %in% d)
str(new_data)
?between
levels(bme$age)
data_more_70 <- new_data[(as.numeric(new_data$age)>7 & round(new_data$diag_1) %in% d) ,]
data_less_70 <- new_data[(as.numeric(new_data$age)<=7 & round(new_data$diag_1) %in% d) ,]
dim(new_data[round(new_data$diag_1) %in% d ,])
idk <- new_data[new_data$diag_1 %in% d,]
dim(data_more_70)
dim(data_less_70)

data_less_70$state <- 1
data_more_70$state <- 2

final <- rbind(data_less_70, data_more_70)
ggplot(data=final, aes(x=final$state, y = final$diag_1, color= final$age))+
  geom_boxplot()


which(any(idk$diag_1 < 800 | idk$diag_1 > 999))
