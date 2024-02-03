tab <- read.table("rawCounts.txt")
head(tab)
tab2 <- tab[,-c(1,2)]
head(tab2)
dim(tab2)
tab3 <- tab2[complete.cases(tab2),]
head(tab3)
tab4 <- round(tab3)
head(tab4)
write.csv(tab4,file="rounded_matrix.csv")
getwd()
