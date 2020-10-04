# Project #3 / fossil data
# Jaechan Park

library(ggplot2)
library(MASS)
library(multcomp)
library(car)
library(polynom)
library(lmtest) # bptest for testing heteroscedasticity
library(gvlma) # global validation of linear model assumptions


rm(list = ls())

fossil<-read.table(file.choose(), header = T)

str(fossil)
attach(fossil)


ggplot(aes(age,strontium.ratio), data=fossil)+geom_point()+geom_smooth()

boxplot(strontium.ratio~age)

boxplot(fossil$age, main="Box Plot for age", xlab="age")
boxplot(fossil$strontium.ratio, main="Box Plot for strontium.ratio", xlab="strontium.ratio")

par(mfrow=c(1,1))
pairs(fossil)
pairs(fossil, panel = function(x,y) {points(x,y); lines(lowess(x,y), col = "red")})


# Linear Modlel
par(mfrow=c(1,1))
fit.linear<-lm(strontium.ratio~age,data=fossil)
plot(age,strontium.ratio,data=fossil)
abline(fit.linear, data=fossil,col="red")
crPlots(fit.linear)
plot(fit.linear)
summary(fit.linear)
shapiro.test(residuals(fit.linear))


# quadratic Model
fit.qua<-lm(strontium.ratio~age+I(age^2),data=fossil)
crPlots(fit.qua)
par(mfrow=c(2,2))
plot(fit.qua)
summary(fit.qua)
shapiro.test(residuals(fit.qua))
gvlma(fit.qua)
bptest(fit.qua)

# cubic Model
fit.cubic<-lm(strontium.ratio~age+I(age^2)+I(age^3),data=fossil)
crPlots(fit.cubic)
plot(fit.cubic)
summary(fit.cubic)
shapiro.test(residuals(fit.cubic))
bptest(fit.cubic)
gvlma(fit.cubic)
summary(fit.cubic.val)




######################################
##Nonparametric method

par(mfrow=c(1,2))
# Degree 1 Loess function polynomial
loess.model1 <- loess(strontium.ratio ~ age, span=0.90, degree=1)
loess.model2 <- loess(strontium.ratio~ age, span=0.70,degree=1)
loess.model3 <- loess(strontium.ratio~ age, span=0.50,degree=1)

# plot the curves

plot(strontium.ratio ~ age, pch=16, cex=0.6)
hat1 <- predict(loess.model1)
lines(age[order(age)], hat1[order(age)], col="red", lwd=2)
hat2 <- predict(loess.model2)
lines(age[order(age)], hat2[order(age)], col="blue", lwd=2)
hat3 <- predict(loess.model3)
lines(age[order(age)], hat3[order(age)], col="green", lwd=2)
legend("bottomleft", title = "loess() model", legend=c("span=.90","span=0.70", "span=0.5"), lwd=2, cex=1, col=c("red","blue","green"))


### Degree 2 Loess function polynomial
loess.model4 <- loess(strontium.ratio ~ age, span=0.90, degree=2)
loess.model5 <- loess(strontium.ratio~ age, span=0.70,degree=2)
loess.model6 <- loess(strontium.ratio~ age, span=0.50,degree=2)

# plot the curves
plot(strontium.ratio ~ age, pch=16, cex=0.6)
hat4 <- predict(loess.model4)
lines(age[order(age)], hat4[order(age)], col="red", lwd=2)
hat5 <- predict(loess.model5)
lines(age[order(age)], hat5[order(age)], col="blue", lwd=2)
hat6 <- predict(loess.model6)
lines(age[order(age)], hat6[order(age)], col="green", lwd=2)
legend("bottomleft", title = "loess() model", legend=c("span=.90","span=0.70", "span=0.5"), lwd=2, cex=1, col=c("red","blue","green"))


# check model assumptions
fit.loess<-loess(strontium.ratio~age, span=0.5,degree=2)
par(mfrow=c(2,2))
qqnorm(residuals(fit.loess))
qqline(residuals(fit.loess))
scatter.smooth(residuals(fit.loess),span = 1,degree=1)
scatter.smooth(fitted(fit.loess),sqrt(abs(residuals(fit.loess))),span=1,degree=1)



# Prediction
t <- c(seq(95, 125, by =5))
t.pred <- predict(fit.loess, t, se = TRUE)
t.upper <- t.pred$fit + qnorm(0.975) * t.pred$se.fit
t.lower <- t.pred$fit - qnorm(0.975) * t.pred$se.fit
data.frame("pred" = t.pred$fit, "lower" = t.lower, "upper" = t.upper)


# Test for nonlinearity
traceS <- fit.loess$trace.hat
SSE0 <- sum(residuals(fit.linear)^2)
SSE1 <- sum(residuals(fit.loess)^2)
n <- dim(fossil)[1]
Fvalue <- ((SSE0 - SSE1) / (traceS - 2)) / (SSE1 / (n - traceS))
Fvalue
Fcrit <- qf(0.95, traceS - 2, n - traceS)
Fcrit
1 - pf(Fvalue, traceS - 2, n - traceS)
# Since the p-value, linearty is rejected.




