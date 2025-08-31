##############################################################################
# change the following two paths to the correct location on your computer:
MY_INPUT_FILE_PATH ="path_to_data_LuHfOdata.csv" # this is the input file 
# it needs to be in csv format with header (containing the variable names) e.g.:
# Age_Ma, d18O, eHf_t
#   -500, 5.10,  0.51
#   -500, 6.70, -7.84
#   ...
#
# set the R working directory (where output files will be saved):
MY_WORKING_DIRECTORY ="path_to_data"

setwd(MY_WORKING_DIRECTORY)

# You need two packages installed: 'RobustLinearReg' and 'zoo' (for timeseries) 
# load the package for Theil Sen regression:
require(RobustLinearReg) 
# load the zoo package for handling/plotting the timeseries
require(zoo)

# set your timestep and window widths
win_width=5   	# window is +/- win_width Myr (can change to suit) NOTE FULL WINDOW IS 2 x win_width (centered on times in 'tout', see below)
increment=5		# timestep

##############################################################################
## if you are using the same data/input parameters as before (Age_Ma, d18O and
## eHf_t) then you shouldn't need to change anything after this point!


# read in the data from the csv file and make time go in the right direction:
data=read.csv(MY_INPUT_FILE_PATH , header=T, stringsAsFactors = F,strip.white = TRUE)
data$Age_Ma=-data$Age_Ma
# this should give you a data frame with "Age_Ma" "d18O"   "eHf_t" 

# find the start and end of the input data (rounding to lower/upper integer value respectively)
start=min(floor(data$Age_Ma)) + win_width # start from the very beginning if you prefer
# calculate  slope, correlation etc. at 'increment' timesteps from 'start' to 'end'
# e.g. first interval is start +/- win_width
end=max(ceiling(data$Age_Ma))-win_width 

# generate a sequence of times for calculating the moving regression 
# regular 'increment' steps, from 'start' to 'end'
# e.g. -875 -870 -865 -860 -855 ...  
tout= seq(start,end,increment)

##############################################################################
# Define our function to calculate the rolling regression parameters
# given x,y,t : vectors of input data 
# tout: vector of timesteps at which output is calculated
# width: half window width (the window at time t spans t +/- width so it is symmetric)
# filename: name of output pdf file for individual timestep regression plots 
# when you call the function, set filename=NULL if you don't want to generate the plots
# This can be useful if you are using a lot of timesteps (small increment)
# names: 'clean' variable names for plot axes
# lm: type of linear regression to use (lm - default, standard linear regression model, ts - theil sen regression)
# you shouldn't need to change this function

####################
## ROLLregression ##
####################

## rolling regression function for a symmetric window,  
## t>=tout-width; t<=tout+width etc for width=win_width (defined below)
## computed at regular timesteps along "tout"
## 
## slope, R squared plus added in a few other measures for fun 
## adjusted R2, correlations - Pearson (default) and Spearman rank correlation.
## ALSO - added in a step to print all the individual scatters and regression 
## to a PDF for checking, if you give a non-NULL filename
## can run with either standard Linear Regression or the more robust Theil Sen Regression 
## use lm="lm" for standard and "ts" for Theil Sen


rollregression_eq = function(x,y,t,tout,width, filename, names,lm) {
  xlab=names[1]
  ylab=names[2]
  xlim=c(min(x), max(x))
  ylim=c(min(y), max(y))
  
  if(lm!="ts") lm="lm"
  print(paste("regression: ", lm))
  if(!is.null(filename)) { 
    pdf(filename, width=10, height=15)
    par(mfrow = c(6,4))
  }
  out = data.frame(MidAge_Ma=tout, slope=numeric(length(tout)),R2=numeric(length(tout)),AdjR2=numeric(length(tout)), Correlation=numeric(length(tout)), var=numeric(length(tout)),Spearman=numeric(length(tout)))
  for( i in seq_along(tout) ) {
    tempy=y[(t >= (tout[i]-width)) & (t <= (tout[i]+width))]
    tempx=x[(t >= (tout[i]-width)) & (t <= (tout[i]+width))]
    if(lm=="ts") tempmodel = theil_sen_regression(tempy ~ tempx) 
    else tempmodel = lm(tempy ~ tempx)
    out$slope[i] =tempmodel$coefficients[[2]]
    out$R2[i]=summary(tempmodel)$r.squared
    out$AdjR2[i]=summary(tempmodel)$adj.r.squared
    out$Correlation[i]=cor(tempx,tempy)
    out$var[i]=var(tempx,tempy)
    out$Spearman[i]=cor(tempx,tempy, method="spearman")
    
    # TEST PLOT 
    if(!is.null(filename)){
      summarytxt=paste0("Slope = ",formatC(out$slope[i], digits = 3, format = "f"), "; R^2 = ",formatC(out$R2[i], digits = 3, format = "f"))
      plot(tempx,tempy, main=paste("From",(tout[i]-width),"to",(tout[i]+width),"; Mid:",tout[i]),sub=summarytxt, xlab=xlab,ylab=ylab , xlim=xlim, ylim=ylim)
      if(lm=="ts") abline(theil_sen_regression(tempy ~ tempx))
      else abline(lm(tempy ~ tempx))
      
    }
  }
  
  if(!is.null(filename)) dev.off()
  
  return(out)
}

####################
#################### end of function definition
##############################################################################


# call the function - here we compare the results using standard linear regression (first call) and theil sen regression (second call)
# the function is called with a list of input parameters
# this makes it easy to repeat the calcs with different inputs, timesteps, window width etc.
# the function returns a data frame with window midpoints (MidAge_Ma) and the various regression parameters for each timestep/window
# Function inputs:
# the input data are: x=data$d18O,y=data$eHf_t,t=data$Age_Ma
# our chosen output timesteps: tout
# half window width: win_width
# filename is the name of the output filename for plotting (which I have labeled according to the type of regression -LM or TS - and full window width)
# this should be saved in the current working directory
# names are just the names I want to use for the plot axes (x and y) names=c("d18O","eHf")
# lm is the regression model (="lm" or "ts")

regressionOP=rollregression_eq(x=data$d18O,y=data$eHf_t,t=data$Age_Ma,tout=tout,width=win_width, filename=paste0("TEST_LM_REGRESSIONS_5MyrData_",2*win_width,"Myr_Window.pdf"), names=c("d18O","eHf"),lm="lm" )

regressionOP_TS=rollregression_eq(x=data$d18O,y=data$eHf_t,t=data$Age_Ma,tout=tout,width=win_width, filename=paste0("TEST_TS_REGRESSIONS_5MyrData_",2*win_width,"Myr_Window.pdf"), names=c("d18O","eHf"),lm="ts" )

# you can just save regressionOP and regressionOP_TS data frames to a file if you like
# using:
#  write.csv(regressionOP, "regressionOP.csv", quote=FALSE, row.names=FALSE)
#  write.csv(regressionOP_TS, "regressionOP_TS.csv", quote=FALSE, row.names=FALSE)
# I convert the output to a timeseries object (using the zoo package) 
# as it is easier to make a nice clean plot this way


# check the output:
plot(regressionOP$MidAge_Ma,regressionOP$Correlation) # 
plot(regressionOP$MidAge_Ma,regressionOP$var) 
plot(regressionOP$MidAge_Ma,regressionOP$Spearman) 

plot(regressionOP$MidAge_Ma,regressionOP$slope, main="Linear Regression") #, type="l")
plot(regressionOP$MidAge_Ma,regressionOP$R2, main="Linear Regression") #, type="l", col="red")
plot(regressionOP$MidAge_Ma,regressionOP$AdjR2, main="Linear Regression") #, type="l", col="darkred")

plot(regressionOP_TS$MidAge_Ma,regressionOP_TS$slope, main="Theil Sen Regression") #, type="l")
plot(regressionOP_TS$MidAge_Ma,regressionOP_TS$R2, main="Theil Sen Regression") #, type="l", col="red")
plot(regressionOP_TS$MidAge_Ma,regressionOP_TS$AdjR2, main="Theil Sen Regression") #, type="l", col="darkred")


# make zoo timeseries of all the output parameters (this just makes it easier to plot, you don't have to do this):
Slope=zoo(regressionOP$slope,regressionOP$MidAge_Ma)
Rsquared=zoo(regressionOP$R2,regressionOP$MidAge_Ma)
Correlation=zoo(regressionOP$Correlation,regressionOP$MidAge_Ma)
SpearmanRankCorr=zoo(regressionOP$Spearman,regressionOP$MidAge_Ma)

# merge into a single multivariate timeseries for tidy plotting
Results=merge.zoo(Slope, Rsquared, Correlation, SpearmanRankCorr)

# save the results to a csv file:
write.zoo(Results, paste0("Linear_Regression_etc_",-1*start,"Ma_",2*win_width,"Myr_movingWindow_060923.txt"), quote=F, sep=",")

#plot all the timeseries together on one plot (saved to working directory) - plots as a smooth line
pdf(paste0("Linear_Regres_Corr_",-1*start,"Ma_",2*win_width,"Myr_movingWindow_060923.pdf"), width=6, height=10)
plot(Results, main=paste("Moving window:", 2*win_width, "Myr"), xlab="midpoint (Myr)")
dev.off()

# using  type="p" plot as points instead of a line:
pdf(paste0("Linear_Regres_Corr_",-1*start,"Ma_",2*win_width,"Myr_movingWindow_points_060923.pdf"), width=6, height=10)
plot(Results, main=paste("Moving window:", 2*win_width, "Myr"), xlab="midpoint (Myr)", type="p", pch=19)
dev.off()


## repeat plots for theil sen regression
Slope=zoo(regressionOP_TS$slope,regressionOP_TS$MidAge_Ma)
Rsquared=zoo(regressionOP_TS$R2,regressionOP_TS$MidAge_Ma)
Correlation=zoo(regressionOP_TS$Correlation,regressionOP_TS$MidAge_Ma)
SpearmanRankCorr=zoo(regressionOP_TS$Spearman,regressionOP_TS$MidAge_Ma)

Results_TS=merge.zoo(Slope, Rsquared, Correlation, SpearmanRankCorr)

# write all the theilsen results to a file:
write.zoo(Results_TS, paste0("ThielSen_Regression_etc_",-1*start,"Ma_",2*win_width,"Myr_movingWindow_060923.txt"), quote=F, sep=",")

pdf(paste0("ThielSen_Regres_Corr_",-1*start,"Ma_",2*win_width,"Myr_movingWindow_060923.pdf"), width=6, height=10)
plot(Results_TS, main=paste("Moving window:", 2*win_width, "Myr"), xlab="midpoint (Myr)")
dev.off()

pdf(paste0("ThielSen_Regres_Corr_",-1*start,"Ma_",2*win_width,"Myr_movingWindow_points_060923.pdf"), width=6, height=10)
plot(Results_TS, main=paste("Moving window:", 2*win_width, "Myr"), xlab="midpoint (Myr)", type="p", pch=19)
dev.off()



## optionally - plot the input timeseries d18O and eHf_t and save as a pdf to check the inputs are all correct:

pdf("InputDataOrig_060923.pdf", width=8, height=8)
par(mfrow = c(2,1))
plot(data$Age_Ma, data$d18O,type="l") #,type="b", pch=20, xlim=c(-700,0))
plot(data$Age_Ma, data$eHf_t, type="l") #, type="b",pch=20, xlim=c(-700,0))
dev.off()






