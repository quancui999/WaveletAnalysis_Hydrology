# this is a recap of workshop about wavelet analysis on hydrology
# at WSML, University of Alberta 
#  Quan Cui, 2020-10-30

# comparing with the original version, this one:
# 1) delete the examples about step-by-step construct data and WA
# 2) use the HYDAT station data directly as an example
# 3) add the part about MannKendall trend detect. The WA package already includes a de-trend part. 
#    So this MannKendall trend detect part will help you if you are specially interested about describe the trend of your data.



# 4 An example of importing my own data, then wavelet analysis, plot, reconstruct---------------------
setwd("~/Wavelet_Quan") #mylaptop
# setwd("D:/My_R_Projects/waveletWorkshop")
rm(list = ls()) 

library(xts)         # to transform daily data to monthly data
library(readr)       # read my downloaded .csv data
library(Kendall)     # Mann-Kendall trend test
library(WaveletComp) # wavelet analysis
library(forecast)    # detect the lag value for detrend

#  1. Import and process data ------------------------------------------------

DailyIm <- read_csv("05UE005_Daily_Flow_ts.csv")  # I download streamflow from Environment Canada, HYDAT. station name 05UE005
date.lookup <- format(seq(as.Date("1960-07-01"), as.Date("2014-12-31"), by = "1 day"))  # based on your data's time span

# import daily data, make monthly mean data
MyDayDF <- DailyIm[, 4]   # daily data itself
myXtsDay <- xts(MyDayDF, order.by=as.Date(date.lookup)) 
myMonthXts <- apply.monthly(myXtsDay, mean)    # mean or sum
MyMonthDF <- as.data.frame(myMonthXts)

# change column names. These 2 are going to be used in Wavelet analysis
colnames(MyDayDF) <- "x"               
colnames(MyMonthDF) <- "y"

plot(myXtsDay, col = "blue")                      # plot daily data - myXtsDay                          
plot(MyMonthDF[,"y"], type = "l", col = "red")    # plot monthly data - MyMonthDF



# 2 Trend detect -----------------------------------------------------------

# Mann-Kendall test 
# For the time series given in Example 1, test the significance of trend using (a)
# Mann-Kendall Test and (b) Kendall's tau test at 10% significance level.

# Null Hypothesis: Time Series does not have a trend.
# Alternative Hypothesis: Time Series has a trend.
# Level of Significance: α = 10%

mydata <- MyMonthDF[,"y"]  # use monthly data as example for trend testing
  
res = MannKendall(mydata$x) ## tau is a measure of the strength of the trend
res

mkStat = res$S
mkStat

mkStatVar = res$varS
mkStatVar

mkPValue = res$sl
mkPValue

Z = (mkStat-1)/sqrt(mkStatVar)
Z

Zlower = qnorm(0.05) # α/2
Zlower

Zupper = qnorm(0.95) # 1- α/2
Zupper

# As |Z| > 1.6450, so the null hypothesis of no trend is rejected.

# Statistically significant?


# 3 Detrend data ------------------------------------------------------------
library(forecast)

# mydata <- MyMonthDF[,"y"] # use monthly data for detrend example
mydata <- CO2$average[361:700]  # take the late half will be enough for anlaysis purpose

plot.ts(mydata)

## twice-difference 
mydata.lag2 = diff(mydata, differences = 2)   # still use monthly data as example for detrend
plot.ts(mydata.lag2 , ylab = expression("Detrended monthly data, lag=2"))


# find out the lag=? in the detrend step
Acf(mydata.lag2)
# read the plot, the distance between peaks (x=12, and x=24) = 12 time units

## difference the differenced CO2 data
mydata.lag12 <- diff(mydata.lag2, lag = 12)   # 12 was read from ACF result plot
## plot the newly differenced data
plot.ts(mydata.lag12 , ylab = expression("Detrended monthly data, lag=12"))
Acf(mydata.lag12) # not useful

# So, mydata.lag12 is the de-trended monthly data

mydata.lag122 = diff(mydata, differences = 12) # use the original monthly value, did not lag2
plot.ts(mydata.lag122 , ylab = expression("Detrended monthly data, lag=2"))
Acf(mydata.lag122)



# 4 Wavelet analysis --------------------------------------------------------
# wavelet analysis and reconstruc was independent from trend analysis. 
# If your data has a trend, set loess.span = 0.75. If no trend, loess.span = 0.

# Replace NA with 0 (you may use other values)
match("TRUE", is.na(MyMonthDF))   # replace month data NA to 0 (or change 0 to other values as you wish)
# shows how many non- NA in your dataset, if it less than data size, you have NA
MyMonthDF[is.na(MyMonthDF)] <- 0
match("TRUE", is.na(MyMonthDF))
# shows NA, means no Null value now in your

# Set up upperlimit using log2 in R
# if your target testing period is 3 years, for daily data
2^ceiling(log2(12*20))


# wavelet analysis comparing daily and monthly
my.dataD <- data.frame(x = MyDayDF)

my.dataM <- data.frame(Y = MyMonthDF)


my.wx <- analyze.wavelet(my.dataD, "x", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 8192,
                         make.pval = TRUE, n.sim = 10)
my.wy <- analyze.wavelet(my.dataM, "y", loess.span = 0.75,
                         dt = 1, dj = 1/250,
                         lowerPeriod = 2, upperPeriod = 256,
                         make.pval = TRUE, n.sim = 10)

wt.image(my.wx, n.levels = 250,
         legend.params = list(lab = "wavelet power levels") )
wt.image(my.wy, n.levels = 250,
         legend.params = list(lab = "wavelet power levels") )

# plot average power
wt.avg(my.wx)
wt.avg(my.wy)

# reconstruct
reconstruct(my.wx, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")
reconstruct(my.wy, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

# Find exact periods. use 'my.wy' the wavelet analysis result of monthly data 
Test.y <- as.numeric(my.wy$Power.avg)
Test.x <- as.numeric(2^(my.wy$axis.2))
plot(x = Test.x, y = Test.y, type="l")

library(quantmod)

MyPeaks <- findPeaks(Test.y)    # this includes all peaks (only 4 out of the 5 is significant)
MyPeriodResult <- Test.x[MyPeaks]
MyPeriodResult                  # pick out the 2rd, 3rd, 4th, 5th one as periods that significant








