#install.packages("DirichletReg")
library(DirichletReg)

Npop   <- read.csv(file='Npop.csv',header=FALSE)
Npop$Y <- DR_data(Npop)

matplot(t(Npop[,1:19]),type="l")

fit    <- DirichReg(Y~1,Npop)
alpha  <- as.matrix(fit$coefficients)

marq25 <- qbeta(0.25,alpha,sum(alpha)-alpha,ncp=0,lower.tail=TRUE,log.p=FALSE)
marq50 <- qbeta(0.50,alpha,sum(alpha)-alpha,ncp=0,lower.tail=TRUE,log.p=FALSE)
marq75 <- qbeta(0.75,alpha,sum(alpha)-alpha,ncp=0,lower.tail=TRUE,log.p=FALSE)
points(seq(1,19,1),marq25,pch=25,bg="black",col="black")
points(seq(1,19,1),marq50,pch=15,bg="black",col="black")
points(seq(1,19,1),marq75,pch=24,bg="black",col="black")

write.table(fit$coefficients,file='alpha.csv',row.names=FALSE,col.names=FALSE)