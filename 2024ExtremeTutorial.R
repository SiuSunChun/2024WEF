##install.packages('terra') ###for your first run you need to install this packages
library(terra)

###Loading data
HD<-read.csv('https://raw.githubusercontent.com/SiuSunChun/LearnRwithSun/main/Data/A0_Historical.csv')
FF1<-read.csv('https://raw.githubusercontent.com/SiuSunChun/LearnRwithSun/main/Data/A1_RCP85_1980_2000.csv')
FF2<-read.csv('https://raw.githubusercontent.com/SiuSunChun/LearnRwithSun/main/Data/A1_RCP85_2020_2040.csv')
FF3<-read.csv('https://raw.githubusercontent.com/SiuSunChun/LearnRwithSun/main/Data/A1_RCP85_2060_2080.csv')
UK1<-readRDS(gzcon(url('https://github.com/SiuSunChun/LearnRwithSun/raw/main/Data/UKMap2.rds')))
UK2<-readRDS(gzcon(url('https://github.com/SiuSunChun/LearnRwithSun/raw/main/Data/UKMap4.rds')))

###UK Catchment Maps
plot(UK1,border=grey(0.6))
plot(UK2,add=T,border=grey(0.6))
axis(1)
axis(2)
box()
text(coordinates(UK1)[,1],coordinates(UK1)[,2],UK1$geo_region,
pos=c(1,1,1, 3,1,1, 1,1,3,
3,1,1, 1,3,3, 1,1,1,
1,1,1,1,1),
offset=0.2,
cex=0.7)

##Avaliable Years
HD[,1]

##Catchment Names
names(HD)[-1]

###Select Catchment
#SeC<-'Severn'
SeC<-'Thames'

###Historical annual precipitation maximum time series 
#options(repr.plot.width = 10, repr.plot.height = 5)
plot(HD$Year,HD[,SeC],type='h',xlab='Year',ylab='Annual Precipitation Maximum (mm)')
title(SeC,line=0.5)

###UKCP18 time series for 3 time periods 
###RCP85_1980_2000
###RCP85_2020_2040
###RCP85_2060_2080
#options(repr.plot.width = 10, repr.plot.height = 5)
plot(c(FF1$Year,FF2$Year,FF3$Year),c(FF1[,SeC],FF2[,SeC],FF3[,SeC]),
type='h',xlab='Year',ylab='Annual Precipitation Maximum (mm)')
title(SeC,line=1.5)
title('UKCP18 RCP85',line=0.5)

####Extreme analysis###
##install.packages('ismev') ###for your first run you need to install this packages
library('ismev')

#####Generialised Extreme Value (GEV) quantile functions#######
qgevs<-function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
    if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 
        1) 
        stop("`p' must contain probabilities in (0,1)")
    if (min(scale) < 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if (!lower.tail) 
        p <- 1 - p
    if (shape == 0) 
        return(loc - scale * log(-log(p)))
    else return(loc + scale * ((-log(p))^(-shape) - 1)/shape)
}

#####Generialised Extreme Value (GEV) probability functions#######
pgevs<-function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
    if (min(scale) <= 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    q <- (q - loc)/scale
    if (shape == 0) 
        p <- exp(-exp(-q))
    else p <- exp(-pmax(1 + shape * q, 0)^(-1/shape))
    if (!lower.tail) 
        p <- 1 - p
    p
}

###Extreme analysis plots
#options(repr.plot.width = 8, repr.plot.height = 8)
gev.diag(gev.fit(HD[,SeC],show=F))
title(SeC)

###Fitted model
F1<-gev.fit(HD[,SeC],show=F)

HH<-hist(HD[,SeC],plot=F)
XX0<-seq(min(HH$breaks),max(HH$breaks),length=100)
YY0<-gev.dens(F1$mle, XX0)

#options(repr.plot.width = 6, repr.plot.height = 6)
plot(HH,ylim=range(HH$density,YY0),freq=F,main=SeC,xlab='Annual Precipitation Maximum (mm)')
lines(XX0,YY0)

###UKCP18
###Generate extreme model data fits
PP1<-gev.fit(FF1[,SeC],show=F)
PP2<-gev.fit(FF2[,SeC],show=F)
PP3<-gev.fit(FF3[,SeC],show=F)

###Generate the fitted extreme curves
XX<-seq(min(c(FF1[,SeC],FF2[,SeC],FF3[,SeC]))*0.8,max(c(FF1[,SeC],FF2[,SeC],FF3[,SeC])),length=100)
YY1<-gev.dens(PP1$mle, XX)
YY2<-gev.dens(PP2$mle, XX)
YY3<-gev.dens(PP3$mle, XX)

####plot density plots
#options(repr.plot.width = 6, repr.plot.height = 6)
par(mar=c(4,4.5,3,1.5))
plot(XX,YY1, ylim =range(YY1,YY2,YY3),main='',xlab='Precipitation extreme (mm)',col='brown',lwd=2,ty='l',ylab='Density')
lines(XX,YY2,col='orange',lwd=2)
lines(XX,YY3,col='red',lwd=2)
title(SeC,line=0.5)
legend('topright',c('1980-2000','2020-2040','2060-2080'),lwd=2,col=c('brown','orange','red'),cex=0.8)

####The 1 in 100 year event for the annual precipitation maximum at the Thames catchment is 46.86mm in 1980-2000
Q1<-qgevs(1-1/100,PP1$mle[1],PP1$mle[2],PP1$mle[3])
print(paste0('The ',SeC,' results'))
print(paste0('The 1 in 100 year event is ',round(Q1,2),'mm for the UKCP18 simulation in 1980-2000'))
####It will be only 1 in 35 year event for the annual precipitation maximum at the Thames catchment in 2060-2080 under the RCP8.5 scenarios
print(paste('1 in',round(1/(1-pgevs(Q1,PP3$mle[1],PP3$mle[2],PP3$mle[3])),0),'year in  2060-2080under the RCP8.5 scenarios'))

