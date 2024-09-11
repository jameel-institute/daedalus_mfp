# suppressWarnings(suppressMessages({require(voi,quietly=T)
#   require(scales,quietly=T)
#   require(earth,quietly=T)
#   require(pracma,quietly=T)
#   require(heavy,quietly=T)}))

combos       <- list(c('av_age','av_lfx'),
                     c('Pwfp','Ppri','Psec','Pter',
                       'WFHpri','WFHsec','WFHter'),
                     c('comm','travelA3','schoolA2','workp'),
                     c('Gpri','Gsec','Gter'),
                     c('Tres','t_tit','trate'),
                     c('sdl','sdb'),
                     c('Hmax'),
                     c('t_vax','arate','puptake'))
diseases     <- list('Influenza 2009','Influenza 1957','Influenza 1918',
                     'Covid Omicron','Covid Wildtype','Covid Delta','SARS');
strategies   <- list('No Closures','School Closures',
                     'Economic Closures','Elimination');
informations <- array(0,
                dim=c(length(combos),length(diseases),length(strategies)))
#pairs     <- combos
#for (i in 1:(length(combos)-1)){
#for (j in (i+1):length(combos)){
#pairs <- append(pairs,list(c(combos[[i]],combos[[j]])))
#}}

for (j in 1:length(diseases)){
for (k in 1:length(strategies)){
fname1 <- paste('LLMIC_',diseases[j],'_',strategies[k],'.csv',sep="")
T1     <- read.csv(fname1)
fname2 <- paste('UMIC_',diseases[j],'_',strategies[k],'.csv',sep="")
T2     <- read.csv(fname2)
fname3 <- paste('HIC_',diseases[j],'_',strategies[k],'.csv',sep="")
T3     <- read.csv(fname3)
T      <- rbind(T1,T2,T3)

SEC  <- T$SEC
gdp  <- 365*rowSums(T[, 331:375])
y    <- 100*SEC/gdp
#y    <- SEC
inds <- !is.nan(y)
y    <- y[inds]
y    <- data.frame(y)

x      <- T[which(inds),]
Npop   <- x[, 4:24]
bands  <- seq(from=2,to=102,by=5)
av_age <- (as.matrix(Npop)%*%bands)/rowSums(Npop)
x      <- cbind(x,av_age=av_age)
# Pret   <- rowSums(Npop[,14:21])/rowSums(Npop)
# x      <- cbind(x,Pret=Pret)
# Pwrk   <- rowSums(Npop[,5:13])/rowSums(Npop)
# x      <- cbind(x,Pwrk=Pwrk)
# Psch   <- rowSums(Npop[,2:4])/rowSums(Npop)
# x      <- cbind(x,Psch=Psch)
# Ppre   <- Npop[,1]/rowSums(Npop)
# x      <- cbind(x,Ppre=Ppre)
Npop2  <- cbind(Npop[,1:17],rowSums(Npop[,18:21]))
la     <- x[,476:493]
av_lfx <- diag(as.matrix(Npop2)%*%t(la))/rowSums(Npop)
x      <- cbind(x,av_lfx=av_lfx)
###
NNs    <- x[, 25:69]
Pwfp   <- rowSums(NNs)/rowSums(Npop)
x      <- cbind(x,Pwfp=Pwfp)
Pwfa   <- rowSums(NNs)/rowSums(Npop[,5:13])
x      <- cbind(x,Pwfa=Pwfa)
Ppri   <- rowSums(NNs[,1:5])/rowSums(NNs)
x      <- cbind(x,Ppri=Ppri)
Psec   <- rowSums(NNs[,6:25])/rowSums(NNs)
x      <- cbind(x,Psec=Psec)
Pter   <- rowSums(NNs[,26:45])/rowSums(NNs)
x      <- cbind(x,Pter=Pter)
###
wfhu   <- x[, 421:465]
WFHpri <- diag(as.matrix(NNs[,1:5])%*%t(wfhu[,1:5]))/rowSums(NNs[,1:5])
x      <- cbind(x,WFHpri=WFHpri)
WFHsec <- diag(as.matrix(NNs[,6:25])%*%t(wfhu[,6:25]))/rowSums(NNs[,6:25])
x      <- cbind(x,WFHsec=WFHsec)
WFHter <- diag(as.matrix(NNs[,26:45])%*%t(wfhu[,26:45]))/rowSums(NNs[,26:45])
x      <- cbind(x,WFHter=WFHter)
###
obj    <- x[, 331:375]
Gpri   <- rowSums(obj[,1:5])/rowSums(NNs[,1:5])
x      <- cbind(x,Gpri=Gpri)
Gsec   <- rowSums(obj[,6:25])/rowSums(NNs[,6:25])
x      <- cbind(x,Gsec=Gsec)
Gter   <- rowSums(obj[,26:45])/rowSums(NNs[,26:45])
x      <- cbind(x,Gter=Gter)
#
#
#

informations[,j,k] <- voi::evppivar(y$y,x,par=combos,verbose = F)$evppi/var(y$y)

}}

infoplot <- rbind(informations[,1,],informations[,2,],informations[,3,],
                  informations[,4,],informations[,5,],informations[,6,],
                  informations[,7,])
heatmap(infoplot,Colv=NA,Rowv=NA,scale="column")
# labels <- vector("list",length=length(combos))
# for (i in 1:length(combos)){labels[i] <- paste(combos[[i]], collapse = "")}
# plot(informations,1:length(informations),type="p",pch=16,col="blue",
#      xlim=c(0,1.2),xlab='EVPPI',ylab='P2 Measure(s)')
# text(x=informations+0.2,y=1:length(informations),
#      labels=unlist(labels),srt = 0)
infoplot