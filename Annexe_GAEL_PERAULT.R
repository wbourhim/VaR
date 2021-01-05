#################################################################################
################################## GAEL PERAULT   ##################
#################################################################################


#################################################################################
#########   ANNEXE ACTION EUROFINS SCIENTIFICS    #######
#################################################################################


rm(list=ls(all=TRUE))

library(QuantTools)
library(xts)
library(forecast)
library(moments)
library(BatchGetSymbols)
library(rugarch)
library(ghyp)
first.date <- "2010-01-04" #date debut
last.date <- "2019-12-31"#date fin
freq.data <- 'daily'#frequence
type.return<-'log'#type de rendement
tickers <- 'ERF.PA' #symbole de l'action sur yahoo/finance
tab <- BatchGetSymbols(tickers = tickers,
                       first.date = first.date,
                       last.date = last.date,
                       freq.data = freq.data,
                       type.return=type.return,
                       cache.folder = file.path(tempdir(),
                                                'BGS_Cache') ) # cache in tempdir()
pt<-tab$df.tickers$price.adjusted

dates<-tab$df.tickers$ref.date[-1]
rendement=tab$df.tickers$ret.adjusted.prices[-1]
N<-length(rendement)
rt<-xts(x=rendement,order.by=dates)
length(rt)
rte=rendement[1:1510]#rt sur ensemble estimation de 2010 `a 2015 inclu
datesrte<-dates[1:1510]
datesrtt = dates[1511:N]
rtt=rendement[1511:N]#rt sur l'ensemble de test de 2016 `a 2019 inclu
alpha=0.95
Ne=length(rte)
Nt=length(rtt)

#############################################
####### VAR normal, Historique, Cornish Fisher 
############################################


backTestVaR <- function(x, p = alpha) {
  normal.VaR = as.numeric(VaR(x, p=p, method="gaussian"))
  historical.VaR = as.numeric(VaR(x, p=p, method="historical"))
  modified.VaR = as.numeric(VaR(x, p=p, method="modified"))
  ans = c(normal.VaR, historical.VaR, modified.VaR)
  names(ans) = c("Normal", "HS", "Modified")
  return(ans)
}

# rolling 1-step ahead estimates of VaR
VaR.results = rollapply(as.zoo(rt), width=Ne, 
                        FUN = backTestVaR, p=alpha, by.column = FALSE,
                        align = "right")
#VaR.results = lag(VaR.results, k=-1)


violations.mat = matrix(0, 3, 5)
rownames(violations.mat) = c("Normal", "HS", "Modified")
colnames(violations.mat) = c("En1", "n1", "1-alpha", "Percent", "VR")
violations.mat[, "En1"] = (1-alpha)*Nt
violations.mat[, "1-alpha"] = 1 - alpha

# Show Normal VaR violations
normalVaR.violations = as.numeric(as.zoo(rt[index(VaR.results)])) < VaR.results[, "Normal"]
violation.dates = index(normalVaR.violations[which(normalVaR.violations)])

for(i in colnames(VaR.results)) {
  VaR.violations = as.numeric(as.zoo(rt[index(VaR.results)])) < VaR.results[, i]
  violations.mat[i, "n1"] = sum(VaR.violations)
  violations.mat[i, "Percent"] = sum(VaR.violations)/Nt
  violations.mat[i, "VR"] = violations.mat[i, "n1"]/violations.mat[i, "En1"]
}
violations.mat

resultats<-data.frame(matrix(NA,ncol=5,nrow=3))
colnames(resultats)<-c("expected.exceed","actual.exceed","Kupiecpv","Christoffersenpv","Violation_rate")
rownames(resultats)<-c("Normale","HS","CF")

# normale
VaR.test1 = VaRTest(1-alpha,actual=coredata(rt[index(VaR.results)]), VaR=coredata(VaR.results[,"Normal"]))
resultats[1,1]=VaR.test1$expected.exceed
resultats[1,2]=VaR.test1$actual.exceed
resultats[1,3]=VaR.test1$uc.LRp
resultats[1,4]=VaR.test1$cc.LRp
resultats[1,5]=(VaR.test1$actual.exceed/length(VaR.results))

# historique
VaR.test2 = VaRTest(1-alpha,actual=coredata(rt[index(VaR.results)]), VaR=coredata(VaR.results[,"HS"]))
resultats[2,1]=VaR.test2$expected.exceed
resultats[2,2]=VaR.test2$actual.exceed
resultats[2,3]=VaR.test2$uc.LRp
resultats[2,4]=VaR.test2$cc.LRp
resultats[2,5]=VaR.test2$actual.exceed/length(VaR.results)


# modifie
VaR.test3 = VaRTest(1-alpha, actual=coredata(rt[index(VaR.results)]), VaR=coredata(VaR.results[,"Modified"]))

resultats[3,1]=VaR.test3$expected.exceed
resultats[3,2]=VaR.test3$actual.exceed
resultats[3,3]=VaR.test3$uc.LRp
resultats[3,4]=VaR.test3$cc.LRp
resultats[3,5]=VaR.test3$actual.exceed/length(VaR.results)



###################################################
#### iGARCH Norm AND BACKTESTING
###################################################


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
spec_ig = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="norm",fixed.pars = list(gamma1 = 0 ))
fit = ugarchfit(spec =spec_ig, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit)

roll_ig=ugarchroll(spec_ig, data=rt,n.ahead=1,forecast.length=length(rtt),refit.every=10,
                   refit.window="moving",solver = "hybrid",
                   cluster=cl,calculate.VaR=TRUE,VaR.alpha=0.05,keep.coef = TRUE)




roll_VaR_ig <- zoo(roll_ig@forecast$VaR[, 1])
tail(roll_VaR_ig)

report(roll_ig,type="VaR",VaR.alpha=0.05,conf.level=0.95)

# calcul ES
spec_ig = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="norm",fixed.pars = list(gamma1 = 0 ))
fit = ugarchfit(spec =spec_ig, data = rt,out.sample=length(rtt),solver="hybrid")
spec_a = spec_ig
setfixed(spec_a)<-as.list(coef(fit))
filt = ugarchfilter(spec_a, data=rtt)
f = function(x) qdist("norm",p=x,mu=0,sigma=1,skew=coef(fit)["skew"],shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
VaR_ig=fitted(filt)+sigma(filt)*qdist("norm",p=0.05,mu=0,sigma=1,skew = coef(fit)["skew"],shape=coef(fit)["shape"])
index(VaR_ig) = dates[1511:N]
# EXPECTED SHORTFALL FOR ROLLING WINDOW VAR METHOD
print(ESTest(0.05, rtt, ES, roll_VaR_ig, boot = TRUE))
# EXPECTED SHORTFALL FOR ONE ESTIMATION VAR METHOD
print(ESTest(0.05, rtt, ES, VaR_ig, boot = TRUE))
stopCluster(cl)



###################################################
####### ARCH M  NORM AND BACKTESTING
###################################################


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
spec_arch = ugarchspec(mean.model=list(armaOrder=c(2,2)),distribution.model="norm")
fit = ugarchfit(spec =spec_arch, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit)
roll_arch=ugarchroll(spec_arch, data=rt,n.ahead=1,forecast.length=length(rtt),refit.every=10,
                    refit.window="moving",solver = "hybrid",solver.control=list(tol=1e-6, trace=1),
                    cluster=cl,fit.control=list(scale=1),calculate.VaR=TRUE,VaR.alpha=0.05,keep.coef = TRUE)




roll_VaR_arch <- zoo(roll_arch@forecast$VaR[, 1])

tail(roll_VaR_arch)

report(roll_arch,type="VaR",VaR.alpha=0.05,conf.level=0.95)

# calcul ES
spec_arch = ugarchspec(mean.model=list(armaOrder=c(2,2)),distribution.model="norm")
fit = ugarchfit(spec =spec_arch, data = rt,out.sample=length(rtt),solver="hybrid")
spec_a = spec_arch
setfixed(spec_a)<-as.list(coef(fit))
filt = ugarchfilter(spec_a, data=rtt)
f = function(x) qdist("sstd",p=x,mu=0,sigma=1,skew=coef(fit)["skew"],shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
VaR_arch=fitted(filt)+sigma(filt)*qdist("norm",p=0.05,mu=0,sigma=1,skew = coef(fit)["skew"],shape=coef(fit)["shape"])

# EXPECTED SHORTFALL FOR ROLLING WINDOW VAR METHOD
print(ESTest(0.05, rtt, ES, roll_VaR_arch, boot = TRUE))
# EXPECTED SHORTFALL FOR ONE ESTIMATION VAR METHOD
print(ESTest(0.05, rtt, ES, VaR_arch, boot = TRUE))

stopCluster(cl)

###################################################
### GARCH norm AND BACKTESTING
###################################################


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
spec_g = ugarchspec(distribution.model="norm")
fit = ugarchfit(spec =spec_g, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit)
roll_g=ugarchroll(spec_g, data=rt,n.ahead=1,forecast.length=length(rtt),refit.every=10,
                  refit.window="moving",solver = "hybrid",
                  cluster=cl,calculate.VaR=TRUE,VaR.alpha=0.05,keep.coef = TRUE)



roll_VaR_g <- zoo(roll_g@forecast$VaR[,1])
tail(roll_VaR_g)

report(roll_g,type="VaR",VaR.alpha=0.05,conf.level=0.95)

# calcul ES

spec_a = spec_g
setfixed(spec_a)<-as.list(coef(fit))
filt = ugarchfilter(spec_a, data=rtt)
f = function(x) qdist("norm",p=x,mu=0,sigma=1,skew=coef(fit)["skew"],shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
VaR_g=fitted(filt)+sigma(filt)*qdist("norm",p=0.05,mu=0,sigma=1,skew = coef(fit)["skew"],shape=coef(fit)["shape"])


# EXPECTED SHORTFALL FOR ROLLING WINDOW VAR METHOD
print(ESTest(0.05, rtt, ES, roll_VaR_g, boot = TRUE))
# EXPECTED SHORTFALL FOR ONE ESTIMATION VAR METHOD
print(ESTest(0.05, rtt, ES, VaR_g, boot = TRUE))

stopCluster(cl)

###################################################
#### iGARCH std AND BACKTESTING
###################################################


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
spec_std = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="std",fixed.pars = list(gamma1 = 0 ))
fit = ugarchfit(spec =spec_std, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit)

roll_std=ugarchroll(spec_std, data=rt,n.ahead=1,forecast.length=length(rtt),refit.every=10,
                   refit.window="moving",solver = "hybrid",
                   cluster=cl,calculate.VaR=TRUE,VaR.alpha=0.05,keep.coef = TRUE)




roll_VaR_std <- zoo(roll_std@forecast$VaR[, 1])
tail(roll_VaR_std)

report(roll_std,type="VaR",VaR.alpha=0.05,conf.level=0.95)

# calcul ES
spec_ig = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="norm")
fit = ugarchfit(spec =spec_std, data = rt,out.sample=length(rtt),solver="hybrid")
spec_a = spec_std
setfixed(spec_a)<-as.list(coef(fit))
filt = ugarchfilter(spec_a, data=rtt)
f = function(x) qdist("std",p=x,mu=0,sigma=1,skew=coef(fit)["skew"],shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
VaR_std=fitted(filt)+sigma(filt)*qdist("std",p=0.05,mu=0,sigma=1,skew = coef(fit)["skew"],shape=coef(fit)["shape"])
index(VaR_std) = dates[1511:N]
# EXPECTED SHORTFALL FOR ROLLING WINDOW VAR METHOD
print(ESTest(0.05, rtt, ES, roll_VaR_std, boot = TRUE))
# EXPECTED SHORTFALL FOR ONE ESTIMATION VAR METHOD
print(ESTest(0.05, rtt, ES, VaR_std, boot = TRUE))
stopCluster(cl)

###################################################
#### iGARCH std AND BACKTESTING
###################################################


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
spec_nig = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                      mean.model=list(armaOrder=c(2,2)),distribution.model="nig",fixed.pars = list(gamma1 = 0 ))
fit = ugarchfit(spec =spec_nig, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit)

roll_nig=ugarchroll(spec_nig, data=rt,n.ahead=1,forecast.length=length(rtt),refit.every=10,
                    refit.window="moving",solver = "hybrid",
                    cluster=cl,calculate.VaR=TRUE,VaR.alpha=0.05,keep.coef = TRUE)




roll_VaR_nig <- zoo(roll_nig@forecast$VaR[, 1])
tail(roll_VaR_nig)

report(roll_nig,type="VaR",VaR.alpha=0.05,conf.level=0.95)

# calcul ES
spec_nig = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="nig")
fit = ugarchfit(spec =spec_nig, data = rt,out.sample=length(rtt),solver="hybrid")
spec_a = spec_nig
setfixed(spec_a)<-as.list(coef(fit))
filt = ugarchfilter(spec_a, data=rtt)
f = function(x) qdist("nig",p=x,mu=0,sigma=1,skew=coef(fit)["skew"],shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
VaR_nig=fitted(filt)+sigma(filt)*qdist("nig",p=0.05,mu=0,sigma=1,skew = coef(fit)["skew"],shape=coef(fit)["shape"])
index(VaR_nig) = dates[1511:N]
# EXPECTED SHORTFALL FOR ROLLING WINDOW VAR METHOD
print(ESTest(0.05, rtt, ES, roll_VaR_nig, boot = TRUE))
# EXPECTED SHORTFALL FOR ONE ESTIMATION VAR METHOD
print(ESTest(0.05, rtt, ES, VaR_nig, boot = TRUE))
stopCluster(cl)



#############################
###########Ploting of All Var
############################

plot(dates[1511:N],rtt,type='c',xlab="Dates",ylab="Return/VaR",main = "95% 1 day VaR Backtesting model iGARCH norm")
lines(dates[1511:N],roll_VaR_ig,type='l',col="red")
legend("topright", inset=.05, c("ERF.PA return","VaR (iGARCH - norm)"),col = c("black","red"), lty = c(1,1))

plot(dates[1511:N],rtt,type='b',xlab="Dates",ylab="Return/VaR",main = "95% 1 day VaR Backtesting Model ARCH m  norm")
lines(dates[1511:N],roll_VaR_arch,type='l',col="brown")
legend("topright", inset=.05, c("ERF.PA return","VaR (ARCH - norm)"),col = c("black","red"), lty = c(1,1))

plot(dates[1511:N],rtt,type='b',xlab="Dates",ylab="Return/VaR",main = "95% 1 day VaR Backtesting Model GARCH norm")
lines(dates[1511:N],roll_VaR_g,type='l',col="blue")
legend("topright", inset=.05, c("ERF.PA return","VaR (GARCH - norm)"),col = c("black","red"), lty = c(1,1))

plot(dates[1511:N],rtt,type='b',xlab="Dates",ylab="Return/VaR",main = "95% 1 day VaR Backtesting Model iGARCH nig")
lines(dates[1511:N],roll_VaR_nig,type='l',col="blue")
legend("topright", inset=.05, c("ERF.PA return","VaR (iGARCH - nig)"),col = c("black","red"), lty = c(1,1))

############# THE END


