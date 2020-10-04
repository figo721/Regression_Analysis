## Regression Analysis Project #1 ###
###- Jaechan Park 


library(ggplot2)
library(MASS)
library(multcomp)
library(car)
library(polynom)
library(glmnet)
library(tidyverse)

rm(list = ls())

hospital.full<-read.table(file.choose(), header = T)

# Randomly select 20 observations for the validation set 
set.seed(0727463)
h.test <- sample(1:dim(hospital.full)[1], 20 )
hospital.val <- hospital.full[h.test, ]
hospital <- hospital.full[-h.test, ]
attach(hospital)

x.train <- model.matrix(irisk~., hospital)[,-1]
y.train <- hospital$irisk

x.val<-model.matrix(irisk~.,hospital.val)[,-1]
y.val<-hospital.val$irisk


# reorder columns, "irisk" column will be the first order of the cols.
mydata <- subset(hospital, select=c(3,2,1,4,5,6,7,8,9))


pairs(mydata)
pairs(mydata, panel = function(x,y) {points(x,y); lines(lowess(x,y), col = "red")})

p<- dim(mydata)[2]
cor(mydata[,2:p])


# Check the multicollinerity
fit.null<-lm(y.train~.,data = data.frame(x.train))
vif(fit.null)


## Variable Selection Method ##
fit.full<-lm(irisk~lstay+age+cultratio+xrayratio+nbeds+census+nnurse+facil,data = hospital)
summary(fit.full)
summary(fit.full, correlation = TRUE)
plot(fit.full)


# Backward elimination based on AIC
fit.stepwise.back<-stepAIC(fit.full,list(upper=~lstay+age+cultratio+xrayratio+nbeds+census+nnurse+facil,lower=~1),direction="backward")
# irisk~lstay+cultratio+xrayratio+facil


# Forward selection based on AIC
fit.initial<-lm(irisk~1,data = hospital)
fit.stepwise.forward<-stepAIC(fit.initial,list(upper=~lstay+age+cultratio+xrayratio+nbeds+census+nnurse+facil,lower=~1),direction="forward")
# irisk~cultratio+lstay+facil+xrayratio



# Stepwise method with "both" and starting only intercept.
fit.stepwise.initial<-stepAIC(fit.initial,list(upper=~lstay+age+cultratio+xrayratio+nbeds+census+nnurse+facil,lower=~1),direction="both")
summary(fit.stepwise.initial)
# irisk ~ lstay+cultratio+xrayratio + facil


# Stepwise method with "both" and starting all predictors.
fit.stepwise.full<-stepAIC(fit.full,list(upper=~lstay+age+cultratio+xrayratio+nbeds+census+nnurse+facil,lower=~1),direction="both")
summary(fit.stepwise.full)
# irisk ~ lstay + cultratio + xrayratio + facil


############################################################
# Backward elimination based on F-statistic/t-statistic
dropterm(fit.full, test = "F")
# Remove census
fit1 <- update(fit.full, ~ . - census)
dropterm(fit1, test = "F")
# Remove nbeds
fit2 <- update(fit1, ~ . - nbeds)
dropterm(fit2, test = "F")
# Remove age
fit3 <- update(fit2, ~ . - age)
dropterm(fit3, test = "F")
# Remove nnurse
fit4 <- update(fit3, ~ . - nnurse)
dropterm(fit4, test = "F")
# irisk ~ lstay + cultratio + xrayratio + facil



#############################################################
# Forward selection based on F-statistic/t-statistic
addterm(fit.initial, ~ . + lstay+age+cultratio+xrayratio+nbeds+census+nnurse+facil, test = "F", data=hospital)
# Add lstay
fit1 <- update(fit.initial, ~ . + lstay)
addterm(fit1, ~ . + age+cultratio+xrayratio+nbeds+census+nnurse+facil, test = "F",data=hospital)
# Add cultratio
fit2 <- update(fit1, ~ . + cultratio)
addterm(fit2, ~. +age+xrayratio+nbeds+census+nnurse+facil, test = "F",data=hospital)
# Add facil
fit3<-update(fit2,~ .+facil)
addterm(fit3, ~. +age+xrayratio+nbeds+census+nnurse, test = "F",data=hospital)
# Add xrayratio
fit4<-update(fit3,~ .+xrayratio)
addterm(fit4, ~. +age+nbeds+census+nnurse, test = "F",data=hospital)
# irisk ~ lstay+cultratio+xrayratio+facil


# Model 1
fit.final1<-lm(irisk~lstay+cultratio+xrayratio+facil,data=hospital)
fit.final1
summary(fit.final1)
summary(fit.final1,correlation = TRUE)
par(mfrow=c(2,2))
plot(fit.final1)
crPlots(fit.final1)
# irisk ~ lstay+cultratio+xrayratio+facil


# Transformation
fit.final2<-lm(irisk~lstay+sqrt(cultratio)+xrayratio+log(facil),data = hospital)
summary(fit.final2)
summary(fit.final2,correlation = TRUE)
par(mfrow=c(2,2))
plot(fit.final2)
crPlots(fit.final2)
# Accroding to the p-values of the output, the p-value for xrayratio variable has greater than 0.05.
# Therefore, irisk~lstay+sqrt(cultratio)+log(facil)


PRESS1 <- sum((residuals(fit.final1) / (1 - lm.influence(fit.final1)$hat))^2)
PRESS2 <- sum((residuals(fit.final2) / (1 - lm.influence(fit.final2)$hat))^2)
PRESS1
PRESS2

## Since a model with smaller PRESS value is more proper model, PRESS2, a transformed model,
## has less value than PRESS 1. Therefore, fit.final1 model seems to be considered as a better model.


attach(hospital.val)
# Model 1 with validation set
fit.final1.val <- lm(irisk ~ lstay + cultratio + xrayratio + facil, data = hospital.val)
summary(fit.final1.val)
summary(fit.final1)


# Model 2 with valdiation set
fit.final2.val <- lm(irisk~lstay+sqrt(cultratio)+log(facil), data = hospital.val)
summary(fit.final2.val)
summary(fit.final2)



###################################################
## Ridge Regression ##

lambda.grid<-10^seq(5,-2,length=100)

fit.ridge <- glmnet(x.train, y.train, alpha = 0,family="gaussian",lambda=lambda.grid)
plot(fit.ridge,xvar = "lambda",label=TRUE)


# To achieve reproducible results
set.seed(0727463)
fit.cv<- cv.glmnet(x.train, y.train, alpha = 0)
plot(fit.cv)
# the gray bars at each point show MSE(lambda) plus and minus one SE.

# using the corss-validation results, we now get proper values for lambda
ridge.bestlam<-fit.cv$lambda.min
ridge.lam1se<-fit.cv$lambda.1se


# To calculate the estimates, one needs to perform regression for each lambda. 
fit.ridge.best<-glmnet(x.train,y.train,alpha=0,lambda=ridge.bestlam)
coef(fit.ridge.best)

fit.ridge.1se<-glmnet(x.train,y.train,alpha=0,lambda=ridge.lam1se)
coef(fit.ridge.1se)




###################################################################
## Ridge regression with selected variables

x.train.sel <- model.matrix(irisk~ lstay+cultratio+xrayratio+facil, hospital)[,-1]
y.train.sel <- hospital$irisk

x.val.sel<-model.matrix(irisk~ lstay+cultratio+xrayratio+facil,hospital.val)[,-1]
y.val.sel<-hospital.val$irisk


fit.ridge.sel <- glmnet(x.train.sel, y.train.sel, alpha = 0,family="gaussian",lambda=lambda.grid)
plot(fit.ridge.sel,xvar = "lambda",label=TRUE)


# To achieve reproducible results
set.seed(0727463)
fit.cv.sel<- cv.glmnet(x.train.sel, y.train.sel, alpha = 0)
plot(fit.cv.sel)
# the gray bars at each point show MSE(lambda) plus and minus one SE.

# using the corss-validation results, we now get proper values for lambda
ridge.bestlam.sel<-fit.cv$lambda.min
ridge.lam1se.sel<-fit.cv$lambda.1se
ridge.bestlam.sel
ridge.lam1se.sel


# To calculate the estimates, one needs to perform regression for each lambda. 
fit.ridge.best.sel<-glmnet(x.train.sel,y.train.sel,alpha=0,lambda=ridge.bestlam.sel)
coef(fit.ridge.best.sel)

fit.ridge.1se.sel<-glmnet(x.train.sel,y.train.sel,alpha=0,lambda=ridge.lam1se.sel)
coef(fit.ridge.1se.sel)



#############
## RSS for the train set and RSS for the test(val) set

y.lm.train<-predict(fit.full, newdata = data.frame(x.train))
sum((y.lm.train-y.train)^2)
y.lm.test<-predict(fit.full,newdata=data.frame(x.val))
sum((y.lm.test-y.val)^2)

y.ridge.best.train<-predict(fit.ridge.best, newx=x.train)
sum((y.ridge.best.train-y.train)^2)
y.ridge.best.test<-predict(fit.ridge.best, newx = x.val)
sum((y.ridge.best.test- y.val)^2)

y.ridge.1se.train<-predict(fit.ridge.1se, newx=x.train)
sum((y.ridge.1se.train-y.train)^2)
y.ridge.1se.test<-predict(fit.ridge.1se, newx = x.val)
sum((y.ridge.1se.test- y.val)^2)


y.ridge.best.train.sel<-predict(fit.ridge.best.sel, newx=x.train.sel)
sum((y.ridge.best.train.sel-y.train.sel)^2)
y.ridge.best.test.sel<-predict(fit.ridge.best.sel, newx = x.val.sel)
sum((y.ridge.best.test.sel- y.val.sel)^2)


y.ridge.1se.train.sel<-predict(fit.ridge.1se.sel, newx=x.train.sel)
sum((y.ridge.1se.train.sel-y.train.sel)^2)
y.ridge.1se.test.sel<-predict(fit.ridge.1se.sel, newx = x.val.sel)
sum((y.ridge.1se.test.sel- y.val.sel)^2)

