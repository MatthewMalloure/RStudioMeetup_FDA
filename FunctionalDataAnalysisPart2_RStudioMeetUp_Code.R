#############################################################################
#############################################################################
#               Code for RStudio Webinar: Intro FDA Part II                 #
#    Functional PCA, Functional PCA Regression, & Functional Linear Model   # 
#                   Author: Matthew R. Malloure, Ph.D.                      #
#                  Associate Research Scientist, Dow Inc.                   #
#############################################################################
#############################################################################

##### Section 0: Loading Necessary Packages and Generate Example Plots #####
## Note: Code developed using R version 4.1.3
library(fda)
library(refund)
library(pls)
library(knitr)
# All Data and Plots are Provided with the fda Package
# Berkeley Growth Study Data, Canadian Weather Data, fda Handwriting Data, Pinch Force Data

par(mfrow=c(2,2))
matplot(growth$age,growth$hgtf,type='l',lwd=2,col='black',ylab='Height (in cm)',xlab='Age (Years)',
        main='Berkeley Growth Data (Girls)',cex.lab=1.5,cex.main=2,cex.axis = 1.5)
matplot(CanadianWeather$monthlyTemp,type='l',col='black',lwd=2,xlab='Month',ylab='Temperature (Celsius)',
        main='Canadian Weather Data (Monthly Average)',cex.lab=1.5,cex.main=2,cex.axis = 1.5)
plot(handwrit[, 1, 1], handwrit[, 1, 2], type="l",main='Cursive Handwriting Data',xlab='Rescaled Time (Milliseconds)',
     ylab='Rescaled Time (Milliseconds)',cex.lab=1.5,cex.main=2,cex.axis = 1.5)
for (j in 2:20) lines(handwrit[, j, 1], handwrit[, j, 2])
matplot (pinchtime, pinchraw, type="l", lty=2, cex=2,col=1, lwd=1,  xlab = "Seconds", ylab="Force (N)",
         cex.lab=1.5,cex.main=2,cex.axis = 1.5,main='Pinch Force Data')
abline(h=2, lty=2)

# Example of Smoothing Data Generated from sin Function

set.seed(3918829)
#Generating Observed Functional Observation
x = seq(-pi,pi,length.out = 200)
y = sin(x) + rnorm(length(x),0,.1)

#Construct the Fourier Basis and Penalty Term 
basis = create.fourier.basis(rangeval = range(x),nbasis=100)
Lcoef = c(0,pi^2,0)
harmaccelLfd = vec2Lfd(Lcoef,range(x))

#Plot the Functional Observation
par(mfrow=c(1,1))
plot(x,y,cex.lab=1.25,pch=1,cex=2,col='gray',cex.axis = 1.5,cex.lab = 1.5)
lines(x,sin(x),col='black',lwd=4)

#Select Moderate Smoothing Parameter, Smooth the Data, and Plot
lambda = .001
fdParobj =fdPar(basis,harmaccelLfd,lambda)
sin.fd = smooth.basis(x,y,fdParobj)$fd
lines(sin.fd,lwd=4,lty=2,col='blue')

#Select Small Smoothing Parameter, Smooth the Data, and Plot
lambda = .000000000001
fdParobj =fdPar(basis,harmaccelLfd,lambda)
sin.fd = smooth.basis(x,y,fdParobj)$fd
lines(sin.fd,lwd=4,lty=2,col='red')

#Select Large Smoothing Parameter, Smooth the Data, and Plot
lambda = .1
fdParobj =fdPar(basis,harmaccelLfd,lambda)
sin.fd = smooth.basis(x,y,fdParobj)$fd
lines(sin.fd,lwd=4,lty=2,col='green')
legend('topleft',legend=c('True Function','Under-Smoothed','Smoothed','Over-Smoothed'),
       col=c('black','red','blue','green'),
       lty=c(1,2,2,2),lwd=4,bty='n',cex=1.5)

##### Section 1: Simulation Data #####

# Randomly generate decay curves with unequally spaced & unequal number of measurements

set.seed(49183)
par(mfrow=c(1,1))
k = 41 #Total Number of Curves
control = 11 #Specify the 11th Curve as a 'Control' for Purpose of Example
coefs = runif(k,.05,.2) #Random Coefficients for Exponential Functions
full.domain = c(0,1,2,4,6,8,12,16,20,24) #Domain Specification (Hours)

# Generate, Store (in matrix Y), and Plot the 41 Curves
Y = matrix(0,nrow=k,ncol=length(full.domain))
Y[,1] = rep(100,k) #Require that All Curves Begin at 100% (no decay has occured)
#Generate Curve 1, Adding Random Normal Error: N(0,sd=3)
Y[1,2:ncol(Y)] = 100*exp(-coefs[1]*full.domain[2:ncol(Y)]) + rnorm(ncol(Y)-1,0,3)
num.to.remove = sample(2:5,1) #Determine How Many Domain Points To Exclude
#Randomly Determine Which of the Interior Points to Exclude (Keep 1st and Last)
keeps = c(T,sample(c(rep(T,8-num.to.remove),rep(F,num.to.remove))),T)
Y[1,!keeps] = NA #Set Selected as Missing
Y[1,][Y[1,] < 0] = 0 #Ensure Non-Negative Points
#Plot the 1st Curve
plot(full.domain[keeps],Y[1,keeps],type='l',col='gray',lwd=4,xlab='Time (Hours)',
     ylab='Material Characteristic',cex.lab=1.5,cex.axis=1.5)
#Repeat the Above Steps for Curves 2 to 41
for (j in 2:k) {
  Y[j,2:ncol(Y)] = 100*exp(-coefs[j]*full.domain[2:ncol(Y)]) + rnorm(ncol(Y)-1,0,3)
  num.to.remove = sample(2:5,1)
  keeps = c(T,sample(c(rep(T,8-num.to.remove),rep(F,num.to.remove))),T)
  Y[j,!keeps] = NA
  Y[j,][Y[j,] < 0] = 0
  lines(full.domain[keeps],Y[j,keeps],type='l',col='gray',lwd=4)
}
#To Differentiate the Control from the Other 40 Curves, Color it Black and Add a Legend
lines(full.domain[!is.na(Y[control,])],Y[control,!is.na(Y[control,])],type='l',col='black',lwd=4)
legend('topright',legend=c('Exp. Materials','Standard'),lty=1,lwd=2,col=c('gray','black'),bty='n',cex=2)

##### Section 2: Curve Smoothing #####

# Convert the Discrete Curves to Smooth Curves Using Cubic B-Splines
basis = create.bspline.basis(c(0,24),norder=4,nbasis=20) #Construct Basis Functions
fdobj = fdPar(basis,Lfdobj=2,lambda=10) #Define Functional Parameter with Roughness Penalty

#Individually Smooth Each Curve, Evaluate at Equally Spaced Timepoints, Store in a Matrix
smoothed.curves = matrix(0,nrow=length(seq(0,24)),ncol=nrow(Y))
for (i in 1:nrow(Y)) {
  temp.smooth = smooth.basis(argvals=full.domain[!is.na(Y[i,])],y=Y[i,!is.na(Y[i,])],fdobj)
  smoothed.curves[,i] = eval.fd(seq(0,24),temp.smooth$fd)
}

#Generate a Single Functional Data Object For Use in FPCA
Overall = smooth.basis(seq(0,24),smoothed.curves,fdobj)

## Recreate the Overlay Plot of All Profiles, Now with Smoothed Versions
plot(Overall$argvals,Overall$y[,1],xlab='Time (Hours)',ylab='Smoothed Material Characteristic',
     ylim=c(min(Overall$y),max(Overall$y)),lty=1,lwd=4,col='gray',type='l',cex.lab = 1.5,cex.main=2,cex.axis=1.5)
for (j in 2:ncol(Overall$y)) {
  lines(Overall$argvals,Overall$y[,j],lty=1,lwd=4,col='gray')
}
lines(Overall$argvals,Overall$y[,control],lty=1,lwd=4,col='black')
legend('topright',legend=c('Exp. Materials','Standard'),lty=1,lwd=4,col=c('gray','black'),bty='n',cex=2)

## For Illustration, Compare the Discrete and Smoothed Profiles for the Control Sample
plot(Overall$argvals,Overall$y[,control],xlab='Time (Hours)',ylab='Smoothed Material Characteristic',
     ylim=c(min(Overall$y),max(Overall$y)),lty=1,lwd=4,col='black',type='l',cex.lab=1.5,cex.axis=1.5)
lines(full.domain[!is.na(Y[control,])],Y[control,!is.na(Y[control,])],col='red',lty=2,lwd=4)
legend('topright',legend=c('Discrete Curve','Smooth Curve'),lty=c(2,1),lwd=4,col=c('red','black'),bty='n',cex=2)

##### Section 3: Perform Functional Principal Components Analysis #####

pc = pca.fd(Overall$fd,nharm=2) #Perform FPCA and Store Top 2 PC Functions

#Create the Interpretable Plot of PC Functions (FPCs) (Note Not Presenting Actual FPCs)
par(mfrow=c(1,2))
plot(pc)

pc$varprop #Display Proportion of Variance Explained

#In a 2x2 Grid Plot the Mean Function, Actual FPCs, and FPC Scores ([2,2] Cell Empty)
par(mfrow=c(2,2))
PC.Matrix = eval.fd(seq(0,24),pc$harmonics)
mean.function = as.vector(eval.fd(seq(0,24),pc$meanfd))

plot(seq(0,24),mean.function,lwd=3,ylab='Material Characteristic',xlab='Time (Hours)',
     main='Mean Function',col='black',ylim=c(0,100),type='l',cex.lab=1.5,cex.axis=1.5,cex.main=2)

plot(seq(0,24),PC.Matrix[,1],lwd=3,ylab='',xlab='Time (Hours)',main='PC Functions',col='blue',
     ylim=range(PC.Matrix),type='l',cex.lab=1.5,cex.axis=1.5,cex.main=2)
lines(seq(0,24),PC.Matrix[,2],lwd=3,col='red',lty=2)
legend('bottomright',legend=c('FPC 1','FPC 2'),col=c('blue','red'),lty=c(1,2),lwd=2,bty='n',cex=1.5)
abline(h=0,lwd=2,lty=3)

#Plot PC Scores in 2-Dimensions, Highlighting Control Sample
plot(pc$scores,pch=16,cex=2,col='gray',main='FPC Scores',xlab='FPC Score 1',ylab='FPC Score 2',
     cex.lab=1.5,cex.axis=1.5,cex.main=2)
points(pc$scores[control,1],pc$scores[control,2],pch=16,col='black',cex=2)
legend('top',legend=c('Exp. Materials','Standard'),pch=16,col=c('gray','black'),bty='n',cex=1.5)


##### Section 4: Clustering FPC Scores #####

#Apply K-Means Clustering with K=4 to the PC Scores
fit = kmeans(pc$scores,4)

#Update Scores Plot and Smooth Curve Overlay with Cluster Colors
colors = c('red','blue','green','orange')
par(mfrow=c(1,3))
x = seq(0,24)
plot(pc$scores,col=colors[fit$cluster],pch=16,cex=3,ylab='FPC Score 2',xlab='FPC Score 1',
     main='Clustering FPC Scores',cex.lab=1.5,cex.axis=1.5,cex.main=2)
points(pc$scores[control,1],pc$scores[control,2],col='black',pch=16,cex=3)

plot(x,smoothed.curves[,1],type='l',col=colors[fit$cluster],lwd=3,ylim=c(1,100),
     ylab='Material Characteristic',xlab='Time (Hours)',main='',cex.lab=1.5,cex.axis=1.5,cex.main=2)
for (j in 2:k) {
  lines(x,smoothed.curves[,j],type='l',col=colors[fit$cluster[j]],lwd=3,ylim=c(0,1))
}
lines(x,smoothed.curves[,control],type='l',col='black',lwd=3,ylim=c(0,1))
legend('topright',legend=c('Standard',paste('Cluster:',1:4)),col=c('black',colors),lty=1,bty='n',cex=2,lwd=3)

#Suppose We Wanted to Plot the Nearest j Observed Curves Based on Euclidean Distance in 2-D Space

Plot.Closest = function(j) {
  distances = c()
  for (i in 1:k) distances[i] = sqrt(sum((pc$scores[control,] - pc$scores[i,])^2))
  closest = seq(1:k)[order(distances)[2:(j+1)]]
  plot(full.domain[!is.na(Y[closest[1],])],Y[closest[1],!is.na(Y[closest[1],])],type='l',
       col=colors[fit$cluster[closest[1]]],
       lwd=3,ylim=c(0,100),xlab='Time (Hours)',ylab='Material Characteristic',
       main=paste('Closest',j,'Additives'),cex.lab=1.5,cex.axis=1.5,cex.main=2)
  for (l in 2:length(closest)) {
    lines(full.domain[!is.na(Y[closest[l],])],Y[closest[l],!is.na(Y[closest[l],])],lwd=2,
          col=colors[fit$cluster[closest[l]]])
  }
  lines(full.domain[!is.na(Y[control,])],Y[control,!is.na(Y[control,])],lwd=3,col='black')
}

Plot.Closest(3)

##### Section 5: Exploration of Impact of Various FPC Scores on Functional Data Shape #####

#Store the Mean Function and PC Functions (Karhunen-Loeve Expansion Building Blocks)
mean.function = as.vector(eval.fd(x,pc$meanfd))
pc.functions = eval.fd(x,pc$harmonics) #When Used in this Way, Empirical Basis Functions

#Randomly Select 100 Pairs of New Scores (Restricted to the Range of Observed Scores)

new.scores = cbind(runif(100,min(pc$scores[,1]),max(pc$scores[,1])),runif(100,min(pc$scores[,2]),max(pc$scores[,2])))
new.curves = matrix(0,nrow=100,ncol=25)

#For the First Set of Scores, Plot the following in a 2 x 2 Plot Matrix
#The 1st FPC (Blue Solid Line) with FPC1 * New Score Overlaid (Dashed Blue Line)
#The 2st FPC (Red Solid Line) with FPC2 * New Score Overlaid (Dashed Red Line)
#The Scores Plot of Original Data with Added Point for New Score Pair (Maroon)
#The Resulting New Curve Based on Karhunen-Loeve (Maroon) Overlaid with the Mean Function
par(mfrow=c(2,2))
plot(seq(0,24),pc.functions[,1],type='l',lwd=3,col='blue',ylab='Percentage Performance',
     xlab='Time (Hours)',main='PC Function 1',ylim=c(-20,30),cex.lab=1.5,cex.axis=1.5,cex.main=2)
lines(seq(0,24),pc.functions[,1] * new.scores[1,1],col='blue',lwd=3,lty=2)
plot(seq(0,24),pc.functions[,2],type='l',lwd=3,col='red',ylab='Percentage Performance',
     xlab='Time (Hours)',main='PC Function 2',ylim=c(-8,8),cex.lab=1.5,cex.axis=1.5,cex.main=2)
lines(seq(0,24),pc.functions[,2] * new.scores[1,2],col='red',lwd=3,lty=2)

plot(pc$scores,col='gray',pch=16,cex=3,ylab='FPC Score 2',xlab='FPC Score 1',main='PC Scores',
     cex.lab=1.5,cex.axis=1.5,cex.main=2)
points(pc$scores[control,1],pc$scores[control,2],col='black',pch=16,cex=3)
points(new.scores[1,1],new.scores[1,2],col='maroon',cex=3,pch=13)
new.curves[1,] = mean.function + pc.functions[,1] * new.scores[1,1] + pc.functions[,2] * new.scores[1,2]
plot(x,mean.function,ylim=c(0,100),type='l',lwd=3,ylab='Percentage Performance',xlab='Time (Hours)',
     main='New Decay Curves',cex.main=2,cex.lab=1.5,cex.axis=1.5)
lines(x,new.curves[1,],col='maroon',lwd=3)
Sys.sleep(1)

#For the Remaining New Score Pairs, Progressively Build on Above Plot.
#All Results from Past Iterations are Changed to Light Gray Color. 
#Current Iteration in Full Colors (Blue, Red, and Maroon Respectively)

for (l in 2:100) {
plot(x,pc.functions[,1],type='l',lwd=3,col='blue',ylab='Percentage Performance',
     xlab='Time (Hours)',main='PC Function 1',ylim=c(-20,30),cex.lab=1.5,cex.axis=1.5,cex.main=2)
for (h in 1:(l-1)){
  lines(x,pc.functions[,1] * new.scores[h,1],col='gray',lwd=2,lty=3)
}
lines(x,pc.functions[,1] * new.scores[l,1],col='blue',lwd=3,lty=2)
plot(x,pc.functions[,2],type='l',lwd=3,col='red',ylab='Percentage Performance',
     xlab='Time (Hours)',main='PC Function 2',ylim=c(-8,8),cex.lab=1.5,cex.axis=1.5,cex.main=2)
for (h in 1:(l-1)){
  lines(x,pc.functions[,2] * new.scores[h,2],col='gray',lwd=2,lty=3)
}
lines(x,pc.functions[,2] * new.scores[l,2],col='red',lwd=3,lty=2)
plot(pc$scores,col='gray',pch=16,cex=3,ylab='FPC Score 2',xlab='FPC Score 1',main='PC Scores',
     cex.lab=1.5,cex.axis=1.5,cex.main=2)
points(pc$scores[control,1],pc$scores[control,2],col='black',pch=16,cex=3)
points(new.scores[1:(l-1),1],new.scores[1:(l-1),2],col='gray',cex=3,pch=1)
points(new.scores[l,1],new.scores[l,2],col='maroon',cex=3,pch=13)
new.curves[l,] = mean.function + pc.functions[,1] * new.scores[l,1] + pc.functions[,2] * new.scores[l,2]
plot(x,mean.function,ylim=c(0,100),type='l',lwd=3,ylab='Percentage Performance',xlab='Time (Hours)',
     main='New Decay Curves',cex.lab=1.5,cex.axis=1.5,cex.main=2)
for (h in 1:(l-1)){
  lines(x,new.curves[h,],col='gray',lwd=2,lty=3)
}
lines(x,new.curves[l,],col='maroon',lwd=3)
lines(x,mean.function,lwd=2)
Sys.sleep(1)
}


##### Section 6: Extension of Example to Functional Regression - FPCA Regression #####

#Note the Functional Regression Model We Fit Here is y ~ \int x(t) * beta(t) dt 
# y = Scalar Response
# x(t) = Observed Decay Curves Defined Over Domain T
# beta(t) = Slope Function Defined Over Domain T
# Domain T from 0 to 24 Hours

#To Generate Scalar Response, Numerical Integration is Required and this Composite.Trapezoid Function is Used 
composite.trapezoid <- function(x,y) {
  n = length(x)
  h <- (x[n] - x[1]) / n
  approx = (h/2) * (y[1] + y[n] + 2*sum(y[2:(n-1)]))
  return(approx)
}

#Suppose the True Slope Function is f(x) = sqrt(x) / 500
par(mfrow=c(1,2))
slope.fun = function(y) sqrt(y) / 500

#Simulate Responses Using True Slope Function, True Exponential Decay Function, and Random Normal Error
set.seed(23472)
response = c()
for (i in 1:k) {
  response[i] = composite.trapezoid(x,slope.fun(x) * 100*exp(-coefs[i]*x)) + rnorm(1,0,.5)
}

#Plot the True Slope Function and the Smoothed Profiles Colored By Response Value
par(mfrow=c(1,2))
plot(seq(0,24,.01),slope.fun(seq(0,24,.01)),type='l',lwd=2,main='Slope Function',ylab='',xlab='Time (Hours)',
     cex.axis=1.5,cex.lab=1.5,cex.main=2)

#Define a Color Spectrum Based on Low (Purple) to High (Red) Response Value
pal = colorRampPalette(colors=c('purple','blue','green','yellow','orange','red'))
order.for.coloring = findInterval(response,sort(response))
plot(x,smoothed.curves[,1],type='l',col=pal(length(response))[order.for.coloring[1]],lwd=3,ylim=c(1,100),
     ylab='Material Characteristic',xlab='Time (Hours)',main='Application Performance',cex.lab=1.5,cex.axis=1.5,cex.main=2)
for (j in 2:k) {
  lines(x,smoothed.curves[,j],type='l',col=pal(length(response))[order.for.coloring[j]],lwd=3,ylim=c(0,1))
}
legend('topright',legend=round(sort(response[order.for.coloring %in% round(seq(1,41,length.out=6))]),3),col=pal(6),
       lty=1,bty='n',cex=1.25,lwd=3)


#To Perform FPCA Regression, Predict the Response Using FPC Scores in a Linear Model
par(mfrow=c(1,1))
pca.lm = lm(response ~ pc$scores)
summary(pca.lm) #Print Linear Model Summary
mean(abs(pca.lm$fitted.values - response))

##### Section 7: Extension of Example to Functional Regression - Functional Linear Model #####

#Store the Smoothed Decay Profiles as Matrix (Note We Need to Transform it)
model.matrix = as.matrix(t(Overall$y))

#Construct the Functional Linear Model - Estimate beta(t) Using 15 Cubic Regression Splines
FRModel = pfr(response ~ lf(X=model.matrix,argvals = seq(0,24),bs='cr',k=5))
summary(FRModel)

#Show the Actual by Predicted Plot For Both FPCR and FLM
plot(response,pca.lm$fitted.values,pch=16,cex=1.5,xlab='Actual Performance Property',
     ylab='Predicted Performance Property',main='Functional Regression',
     cex.lab=1.5,cex.axis=1.5,cex.main=2,col='red')
abline(a=0,b=1,lwd=2,col='darkgray')
points(response,FRModel$fitted.values,pch=0,cex=1.5)
legend('topleft',legend=c('FPCR','FLM'),col=c('red','black'),pch=c(16,0),cex=2,bty='n')

#Compare the Estimated Slope Function To the True Slope Function
par(mfrow=c(1,1))
plot(FRModel,lwd=3,ylab='',main='Slope Function',xlab='Time (Hours)',cex.lab=1.5,cex.axis=1.5,cex.main=2)
lines(x,slope.fun(x),col='red',lwd=3)
mean(abs(FRModel$fitted.values - response))

##### Section 8: Comparison of Functional Regression Models #####

#Comparison of Slope Functions Produced With Different Random Seed
#Repeat FLM Model Fitting Using Newly Simulated Data for Each Random Seed
#Store the Slope Function Estimate and MAE to Record Model Performance

set.seed(12945)
seeds = round(runif(10,1,100000),0)
maes = c()

Slope.Functions = list()
for (l in 1:length(seeds)) {
  
set.seed(seeds[l])
response = c()
for (i in 1:k) {
  response[i] = composite.trapezoid(x,slope.fun(x) * 100*exp(-coefs[i]*x)) + rnorm(1,0,.5)
}

model.matrix = as.matrix(t(Overall$y))

Sim.FRModel = pfr(response ~ lf(model.matrix,argvals = x,bs='cr',k=15))
Slope.Functions[[l]] = cbind(plot(Sim.FRModel)[[1]]$x,plot(Sim.FRModel)[[1]]$fit)
maes[l] = mean(abs(Sim.FRModel$fitted.values - response))
}

#Overlay the Slope Functions for the Initial Model and Simulated Random Models with the True Slope Function
par(mfrow=c(1,1))
plot(FRModel,lwd=3,main='Slope Functions',ylab='',xlab='Time (Hours)',
     cex.lab=1.5,cex.axis=1.5,cex.main=2,ylim=c(-.01,.02))
for (l in 1:length(seeds)) {
  lines(Slope.Functions[[l]],lwd=1,lty=1,col='gray')
}
lines(x,slope.fun(x),col='red',lwd=3)
legend('topleft',legend=c('Initial Model','Initial Model CI','True Slope','Simulated Results'),
       col=c('black','black','red','gray'),lty=c(1,2,1,1),lwd=c(3,3,3,1),bty='n')


#Compare FPCR Models with Only FPC1 and FPC1 + FPC2

par(mfrow=c(1,1))
pca.lm1 = lm(response ~ pc$scores[,1])
summary(pca.lm) #Print Linear Model Summary

#Show the Actual By Predicted Plot
plot(response,pca.lm1$fitted.values,pch=16,cex=2,xlab='Actual Performance Property',
     ylab='Predicted Performance Property',main='FPC Regression',
     cex.lab=1.5,cex.axis=1.5,cex.main=2)
points(response,pca.lm$fitted.values,pch=16,col='orange',cex=2)
abline(a=0,b=1,lwd=2,col='darkgray')
mean(abs(pca.lm1$fitted.values - response))
legend('topleft',legend=c('FPC1 & FPC2','FPC1 Only'),col=c('black','orange'),pch=16,cex=2,bty='n')


##### Section 9: Alternative Modeling Approaches (Original Data, No Missing Values) #####

#Regenerate the Data, but This Time Without Missing Values to Create a Level Playing Field.
set.seed(49183)
par(mfrow=c(1,2))
k = 41 #Total Number of Curves
#control = 11 #Specify the 11th Curve as a 'Control' for Purpose of Example
coefs = runif(k,.05,.2) #Random Coefficients for Exponential Functions
full.domain = c(0,1,2,4,6,8,12,16,20,24) #Domain Specification (Hours)
Y = matrix(0,nrow=k,ncol=length(full.domain))
Y[,1] = rep(100,k) #Require that All Curves Begin at 100% (no decay has occured)
#Generate Curve 1, Adding Random Normal Error: N(0,sd=3)
Y[1,2:ncol(Y)] = 100*exp(-coefs[1]*full.domain[2:ncol(Y)]) + rnorm(ncol(Y)-1,0,3)
Y[1,][Y[1,] < 0] = 0 #Ensure Non-Negative Points
#Plot the 1st Curve
plot(full.domain,Y[1,],type='l',col='gray',lwd=3,xlab='Time (Hours)',
     ylab='Material Characteristic',cex.lab=1.5,cex.axis=1.5)
#Repeat the Above Steps for Curves 2 to 41
for (j in 2:k) {
  Y[j,2:ncol(Y)] = 100*exp(-coefs[j]*full.domain[2:ncol(Y)]) + rnorm(ncol(Y)-1,0,3)
  Y[j,][Y[j,] < 0] = 0
  lines(full.domain,Y[j,],type='l',col='gray',lwd=3)
}
#To Differentiate the Control from the Other 40 Curves, Color it Black and Add a Legend
lines(full.domain,Y[control,],type='l',col='black',lwd=3)
legend('topright',legend=c('Exp. Materials','Standard'),lty=1,lwd=3,col=c('gray','black'),bty='n',cex=1.25)

set.seed(23472)
response = c()
for (i in 1:k) {
  response[i] = composite.trapezoid(x,slope.fun(x) * 100*exp(-coefs[i]*x)) + rnorm(1,0,.5)
}


#Compute a Single Numerical Derivative Based On Rise/Run on Raw Data
#Fit a Simple Linear Model and Store the MAE
MAE = RMSE = c()
dv = (Y[,ncol(Y)] - Y[,1]) / (full.domain[length(full.domain)] - full.domain[1])
one.dv.model = lm(response ~ dv)
MAE[1] = mean(abs(response - one.dv.model$fitted.values))
RMSE[1] = sqrt(mean((response - one.dv.model$fitted.values)^2))

#Bin the Domain Into 4, 6 Hour Windows and Compute the Rise/Run Numerical Derivatives
#Fit a Multiple Regression Model with First Order Terms
#Note the Presence of Collinearity Amongst the Derivatives
#Also, Derivatives are Computed on the Smoothed Data for Simplicity

binned.dv = data.frame(matrix(0,nrow=nrow(Y),ncol=3))
binned.dv[,1] = (Y[,5] - Y[,1]) / (full.domain[5] - full.domain[1])
binned.dv[,2] = (Y[,8] - Y[,5]) / (full.domain[8] - full.domain[5])
binned.dv[,3] = (Y[,length(full.domain)] - Y[,8]) / (full.domain[length(full.domain)] - full.domain[8])

colnames(binned.dv) = c('0-6 Hours','6-16 Hours','16-24 Hours')
binned.dv.model = lm(response ~ as.matrix(binned.dv))
MAE[2] = mean(abs(response - binned.dv.model$fitted.values))
RMSE[2] = sqrt(mean((response - binned.dv.model$fitted.values)^2))


#Compute the Area Under the Curve (Larger Value Implies Slower Decay) Using the Raw Data
#Fit a Quadratic Model to Predict Response

AUCs = c()
for (l in 1:nrow(Y)) AUCs[l] = composite.trapezoid(full.domain,Y[l,])
auc.model = lm(response ~ AUCs)
MAE[3] = mean(abs(response - auc.model$fitted.values))
RMSE[3] = sqrt(mean((response - auc.model$fitted.values)^2))

# Apply Both PC Regression and PLS, Selecting 2 Components Each Time
# Note, Both Methods Only Possible on the Smoothed Data

pcr.model = pcr(response ~ Y,validation = 'CV')
MAE[4] = mean(abs(response - pcr.model$fitted.values[,,2]))
RMSE[4] = sqrt(mean((response - pcr.model$fitted.values[,,2])^2))

#Utilize the Exponential Decay Function as a Parametric Basis Function
#Estimate the Unknown Parameter From the Raw Data
#Fit a Quadratic Model to Predict Response from Fit Parameter

fit.exp = function(x,y) {
  
  MSE = function(alpha) {
    return(sum((y - 100*exp(-alpha*x))^2))
  }
  
  return(optimize(MSE,interval=c(.01,2),maximum=F)$minimum)
}

fit.parms = c()
for (l in 1:nrow(Y)) fit.parms[l] = fit.exp(full.domain,Y[l,])
exp.parm.model = lm(response ~ fit.parms + I(fit.parms^2))
MAE[5] = mean(abs(response - exp.parm.model$fitted.values))
RMSE[5] = sqrt(mean((response - exp.parm.model$fitted.values)^2))

#Store the MAE and Relative Errors for the FPCR and FLM Models

# Convert the Discrete Curves to Smooth Curves Using Cubic B-Splines
basis = create.bspline.basis(c(0,24),norder=4,nbasis=20) #Construct Basis Functions
fdobj = fdPar(basis,Lfdobj=2,lambda=10) #Define Functional Parameter with Roughness Penalty
Overall = smooth.basis(full.domain,t(Y),fdobj)
smoothed.curves = eval.fd(seq(0,24),Overall$fd)

pal = colorRampPalette(colors=c('purple','blue','green','yellow','orange','red'))
order.for.coloring = findInterval(response,sort(response))
plot(seq(0,24),smoothed.curves[,1],type='l',col=pal(length(response))[order.for.coloring[1]],lwd=3,ylim=c(1,100),
     ylab='Material Characteristic',xlab='Time (Hours)',main='Application Performance',cex.lab=1.5,cex.axis=1.5,cex.main=2)
for (j in 2:k) {
  lines(seq(0,24),smoothed.curves[,j],type='l',col=pal(length(response))[order.for.coloring[j]],lwd=3,ylim=c(0,1))
}
legend('topright',legend=round(sort(response[order.for.coloring %in% round(seq(1,41,length.out=6))]),3),col=pal(6),
       lty=1,bty='n',cex=1.25,lwd=3)


pc = pca.fd(Overall$fd,nharm=2) #Perform FPCA and Store Top 2 PC Functions
par(mfrow=c(1,1))
pca.lm = lm(response ~ pc$scores)
summary(pca.lm) #Print Linear Model Summary
MAE[6] = mean(abs(pca.lm$fitted.values - response))
RMSE[6] = sqrt(mean((response - pca.lm$fitted.values)^2))


model.matrix = as.matrix(t(Overall$y))
FRModel = pfr(response ~ lf(model.matrix,argvals = full.domain,bs='cr',k=5))
summary(FRModel)
MAE[7] = mean(abs(FRModel$fitted.values - response))
RMSE[7] = sqrt(mean((response - FRModel$fitted.values)^2))


MAE = round(MAE,3)
RMSE = round(RMSE,3)

Methods = c('One Derivative','Binned Derivatives','AUCs','PCR','Exp. Decay','FPCR','FLM')
kable(cbind(Methods,MAE,RMSE),align='c')


##### Section 10: Alternative Modeling Approaches (Extreme Example 10 Curves with 25 Measurements) #####

set.seed(49183)
par(mfrow=c(1,2))
k = 10 
coefs = runif(k,.05,.2) 
full.domain = seq(0,24) 
Y = matrix(0,nrow=k,ncol=length(full.domain))
Y[,1] = rep(100,k) 
Y[1,2:ncol(Y)] = 100*exp(-coefs[1]*full.domain[2:ncol(Y)]) + rnorm(ncol(Y)-1,0,3)
Y[1,][Y[1,] < 0] = 0
plot(full.domain,Y[1,],type='l',col='gray',lwd=3,xlab='Time (Hours)',
     ylab='Material Characteristic',cex.lab=1.5,cex.axis=1.5)
for (j in 2:k) {
  Y[j,2:ncol(Y)] = 100*exp(-coefs[j]*full.domain[2:ncol(Y)]) + rnorm(ncol(Y)-1,0,3)
  Y[j,][Y[j,] < 0] = 0
  lines(full.domain,Y[j,],type='l',col='gray',lwd=3)
}

set.seed(23472)
response = c()
for (i in 1:k) {
  response[i] = composite.trapezoid(x,slope.fun(x) * 100*exp(-coefs[i]*x)) + rnorm(1,0,.5)
}

MAE = RMSE = c()
dv = (Y[,ncol(Y)] - Y[,1]) / (full.domain[length(full.domain)] - full.domain[1])
one.dv.model = lm(response ~ dv)
MAE[1] = mean(abs(response - one.dv.model$fitted.values))
RMSE[1] = sqrt(mean((response - one.dv.model$fitted.values)^2))

binned.dv = data.frame(matrix(0,nrow=nrow(Y),ncol=3))
binned.dv[,1] = (Y[,9] - Y[,1]) / (full.domain[9] - full.domain[1])
binned.dv[,2] = (Y[,17] - Y[,9]) / (full.domain[17] - full.domain[9])
binned.dv[,3] = (Y[,length(full.domain)] - Y[,17]) / (full.domain[length(full.domain)] - full.domain[17])

colnames(binned.dv) = c('0-6 Hours','6-16 Hours','16-24 Hours')
binned.dv.model = lm(response ~ as.matrix(binned.dv))
MAE[2] = mean(abs(response - binned.dv.model$fitted.values))
RMSE[2] = sqrt(mean((response - binned.dv.model$fitted.values)^2))

AUCs = c()
for (l in 1:nrow(Y)) AUCs[l] = composite.trapezoid(full.domain,Y[l,])
auc.model = lm(response ~ AUCs)
MAE[3] = mean(abs(response - auc.model$fitted.values))
RMSE[3] = sqrt(mean((response - auc.model$fitted.values)^2))

pcr.model = pcr(response ~ Y,validation = 'CV')
MAE[4] = mean(abs(response - pcr.model$fitted.values[,,2]))
RMSE[4] = sqrt(mean((response - pcr.model$fitted.values[,,2])^2))

fit.exp = function(x,y) {
  
  MSE = function(alpha) {
    return(sum((y - 100*exp(-alpha*x))^2))
  }
  
  return(optimize(MSE,interval=c(.01,2),maximum=F)$minimum)
}

fit.parms = c()
for (l in 1:nrow(Y)) fit.parms[l] = fit.exp(full.domain,Y[l,])
exp.parm.model = lm(response ~ fit.parms + I(fit.parms^2))
MAE[5] = mean(abs(response - exp.parm.model$fitted.values))
RMSE[5] = sqrt(mean((response - exp.parm.model$fitted.values)^2))

basis = create.bspline.basis(c(0,24),norder=4,nbasis=20) 
fdobj = fdPar(basis,Lfdobj=2,lambda=10) 
Overall = smooth.basis(full.domain,t(Y),fdobj)

pc = pca.fd(Overall$fd,nharm=2) 
pca.lm = lm(response ~ pc$scores) 
summary(pca.lm) 
MAE[6] = mean(abs(pca.lm$fitted.values - response))
RMSE[6] = sqrt(mean((response - pca.lm$fitted.values)^2))


model.matrix = as.matrix(t(Overall$y))
FRModel = pfr(response ~ lf(model.matrix,argvals = full.domain,bs='cr',k=5))
summary(FRModel)
MAE[7] = mean(abs(FRModel$fitted.values - response))
RMSE[7] = sqrt(mean((response - FRModel$fitted.values)^2))

smoothed.curves = eval.fd(seq(0,24),Overall$fd)
pal = colorRampPalette(colors=c('purple','blue','green','yellow','orange','red'))
order.for.coloring = findInterval(response,sort(response))

plot(seq(0,24),smoothed.curves[,1],type='l',col=pal(length(response))[order.for.coloring[1]],lwd=3,ylim=c(1,100),
     ylab='Material Characteristic',xlab='Time (Hours)',main='Application Performance',cex.lab=1.5,cex.axis=1.5,cex.main=2)
for (j in 2:k) {
  lines(seq(0,24),smoothed.curves[,j],type='l',col=pal(length(response))[order.for.coloring[j]],lwd=3,ylim=c(0,1))
}
legend('topright',legend=round(sort(response[order.for.coloring %in% round(seq(1,k,length.out=6))]),3),col=pal(6),
       lty=1,bty='n',cex=1.25,lwd=3)

MAE = round(MAE,3)
RMSE = round(RMSE,3)

Methods = c('One Derivative','Binned Derivatives','AUCs','PCR','Exp. Decay','FPCR','FLM')
kable(cbind(Methods,MAE,RMSE),align='c')


