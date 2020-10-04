# Project #2 / Fuel Consumption
# Jaechan Park

library(ggplot2)
library(GGally)
library(MASS)
library(multcomp)
library(car)
library(polynom)
library(robustbase)
library(lmtest)
library(olsrr)

rm(list = ls())

fuel<-read.table(file.choose(), header = T)

str(fuel)
names(fuel)

par(mfrow=c(1,4))
boxplot(fuel$Income, main="boxplot for Income")
boxplot(fuel$Miles, main="boxplot for Miles")
boxplot(fuel$Tax, main="boxplot for Tax")
boxplot(fuel$Dlic, main="boxplot for Dlic")

ggpairs(fuel)

pairs(fuel)
pairs(fuel, panel = function(x,y) {points(x,y); lines(lowess(x,y), col = "red")})


################################################################################

fit.full<-lm(Fuel~Income+Miles+Tax+Dlic,data = fuel)
fit.full
summary(fit.full)
summary(fit.full, correlation = TRUE)
plot(fit.full)

fit.initial<-lm(Fuel~1,data=fuel)

################################################################################
# Backward elimination based on AIC
fit.stepwise.back<-stepAIC(fit.full,list(upper=~Income+Miles+Tax+Dlic,lower=~1),data=fuel,direction="backward")
summary(fit.stepwise.back)
# Fuel~Income+Miles+Dlic


# Forward selection based on AIC
fit.stepwise.forward<-stepAIC(fit.initial,list(upper=~Income+Miles+Tax+Dlic,lower=~1),data=fuel,direction="forward")
summary(fit.stepwise.forward)
# Fuel~Income+Miles+Dlic


# Stepwise method with "both" and starting only intercept.
fit.stepwise.initial<-stepAIC(fit.initial,list(upper=~Income+Miles+Tax+Dlic,lower=~1),data=fuel,direction="both")
summary(fit.stepwise.initial)
# Fuel~Income+Miles+Dlic


# Stepwise method with "both" and starting all predictors.
fit.stepwise.full<-stepAIC(fit.full,list(upper=~Income+Miles+Tax+Dlic,lower=~1),data=fuel,direction="both")
summary(fit.stepwise.full)
#Fuel~Income+Miles+Dlic


############################################################
# Backward elimination based on F-statistic/t-statistic
dropterm(fit.full, test = "F")
# Remove Miles
fit1 <- update(fit.full, ~ . - Tax,data=fuel)
dropterm(fit1, test = "F")
# Fuel~Income+Miles+Dlic


#############################################################
# Forward selection based on F-statistic/t-statistic
addterm(fit.initial, ~ . +Income+Miles+Tax+Dlic, test = "F", data=fuel)
# Add Dlic
fit1 <- update(fit.initial, ~ . + Dlic,data=fuel)
addterm(fit1, ~ . + Income+Miles+Tax, test = "F",data=fuel)
# Add Income
fit2 <- update(fit1, ~ . + Income,data=fuel)
addterm(fit2, ~ . + Tax+Miles, test = "F",data=fuel)
#ADD Miles
fit3<-update(fit2, ~ . + Miles,data=fuel)
addterm(fit3, ~ . +Tax, test = "F",data=fuel)
# Fuel~Income+Miles+Dlic



# model, Fuel~Income+Miles+Dlic
fit.final1<-lm(Fuel~Income+Miles+Dlic,data=fuel)
fit.final1
summary(fit.final1)
summary(fit.final1, correlation=TRUE)
plot(fit.final1)
crPlots(fit.final1)



PRESS1 <- sum((residuals(fit.final1) / (1 - lm.influence(fit.final1)$hat))^2)
PRESS1


###################################################################
### Detecting infulential outliers ##

# Standardized residuals
fit.final1.stdres<-stdres(fit.final1)
plot(fit.final1.stdres,ylab="Standardized residuals")
abline(h=c(-2.5,2.5),col="red")


# Studentized residuals
fit.final1.studress<-studres(fit.final1)
plot(fit.final1.studress,ylab="Studentized residuals")
abline(h=c(-2.5,2.5),col="red")


# Diagonal elements of hat matrix
fit.final.influence<-influence(fit.final1)
plot(fit.final1.influence$hat,ylab="Diagonal elements of hat matrix")
n<-dim(fuel)[1]
p<-dim(fuel)[2]
abline(h=2*p/n,col="red")


# DFFITS
fit.final1.dffits<-dffits(fit.final1)
plot(fit.final1.dffits,ylab="DFFITS")
abline(h=2*sqrt(p/n),col="red")


# Cook's distance
fit.final1.cook<-cooks.distance(fit.final1)
plot(fit.final1.cook,ylab="Cook's distance")
abline(h=0.1,col="red")


# DFBETAS 
fit.final1.dfbetas<-dfbetas(fit.final1)
fit.final1.dfbetas


####################################################
### With "olsrr" package for detecting outliers ####

par(mfrow=c(2,2))
# Standardized residuals
ols_plot_resid_stand(fit.final1)

# Studentized residuals
ols_plot_resid_stud(fit.final1)

# DFFITS
ols_plot_dffits(fit.final1)

# Cook's distance
ols_plot_cooksd_bar(fit.final1)

# DFBETAS 
ols_plot_dfbetas(fit.final1)

#Studentized Residuals vs Leverage Plot
ols_plot_resid_lev(fit.final1)


#########################################################
## Robust Regression #########

# RLTS (with default, 50% breakdown value)
RLTS<-ltsReg(Fuel~Income+Miles+Tax+Dlic,data=fuel)
summary(RLTS)

# Detection of outliers
plot(RLTS,which="rindex")
plot(RLTS,which="rdiag")


# Standardized residuals 
RLTS.stdres<-RLTS$residuals/RLTS$scale
plot(RLTS.stdres,ylab="Standardized residuals")
abline(h=c(-2.5,2.5),col="red")


# Diagnostic plot
plot(RLTS$RD,RLTS.stdres,xlab="Robust Distance", ylab="Standardized Residuals")
abline(v=sqrt(qchisq(0.975,p-1)),col="red")
abline(h=c(-2.5,2.5),col="red")

############################################################################

# RLTS (20% Breakdown value)
RLTS.20<-ltsReg(Fuel~Income+Miles+Tax+Dlic,data=fuel, alpha=0.80)
summary(RLTS.20)

plot(RLTS.20,which="rindex")
plot(RLTS.20,which="rdiag")

##########################################################################

