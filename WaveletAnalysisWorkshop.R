# Wavelet Analysis 
# Workshop at Watershed Science & Modelling Laboratory, University of Alberta

# Quan Cui, 2020-10-30





# setwd("D:/My_R_Projects/waveletWorkshop/")               # desktop

rm(list = ls()) 


# 1. Introductory examples - combined cycles of time series data -----------

x = seq(0, 30, 0.1)
x

Y1 = sin(x)
plot(x, Y1, type = "l", col="red", xlab = "time", ylab = "Y1")
legend("topright", legend = "Y1 = sin(x)", col="red", lty=1, cex=0.8)

Y2 = sin(2*x)
plot(x, Y2, type = "l", col="blue", xlab = "time", ylab = "Y2")
legend("topright", legend = "Y2 = sin(2x)", col="blue", lty=1, cex=0.8)

Y12 = Y1 + Y2
plot(x, Y12, type = "l", xlab = "time", ylab = "Y1 + Y2")
legend("topright", legend = "Y = sin(x) + sin(2x)", lty=1, cex=0.8)

Y3 = 0.2*x
plot(x, Y3, type = "l", xlab = "time", ylab = "Y3", col=73)
legend("bottomright", legend = "Y3 = 0.2 * x", lty=1, cex=0.8, col=73)

Y123 = Y1 + Y2 + Y3
plot(x, Y123, type = "l", xlab = "time", ylab = "Y1 + Y2 + Y3")
legend("bottomright", legend = "Y' = sin(x) + sin(2x) + 0.1x", lty=1, cex=0.8)

# sin with increasing preiods, use a function
library(WaveletComp)
Ychange = periodic.series(start.period = 20, end.period = 100, length = 1000)
plot(Ychange, type = "l", xlab = "time", ylab = "Y")



# 2. Wavelet Analysis - Package 'WaveletComp' examples -------------------------
rm(list = ls()) 
library(WaveletComp)

#  2.1 Example 1: a series with constant period---------------------------------

x = periodic.series(start.period = 50, length = 1000)
x = x + 0.2*rnorm(1000)                # add some noise
my.data <- data.frame(x = x)

plot(x, type = "l")

# wavelet analysis. Function "analyze.wavelet". 
# mydata is the dataframe, x is the name of column
# loess.span = 0,     0 means we don't need to detrend this series
# dt=1/24 means 24 data per time unit 
# dj = 1/250 means testing periods was set to [2^n; 2^n+1)/250 
# lowerPeriod and upperPeriod are the limits for testing periods, has to be 2^n
# make.pval = TRUE p-value is 0.1. significant value will circled in white line

my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0.75, 
                        dt = 1, dj = 1/500, 
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)

# plot the result of wavelet analysis
wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7))     
     
# use the period to reconstruct the function, and plot it 
reconstruct(my.w, plot.waves = FALSE, lwd = c(1,2),
            legend.coords = "bottomleft", ylim = c(-1.8, 1.8))

# plot average power     
wt.avg(my.w)   

#  2.2 Example 2: a series with variable period --------------------------------

x = periodic.series(start.period = 20, end.period = 100, length = 1000)
x = x + 0.2*rnorm(1000)     
plot(x, type = "l")  
     
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 8,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
# make.pval p-value is tested against (by default) 0.1 (white line circles the red area)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))     
# black line is the a ridge in the "landscape" of power.    
     
my.rec <- reconstruct(my.w)
x.rec <- my.rec$series$x.r     # x: original series

# plot average power  
wt.avg(my.w)

#  2.3 Example 3: a series with two periods-------------------------------------

x1 <- periodic.series(start.period = 80, length = 1000)
x2 <- periodic.series(start.period = 30, length = 1000)
x <- x1 + x2 + 0.2*rnorm(1000)

plot(x, type="l")

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels") )

# reconstruct 1, use the period=80
reconstruct(my.w, sel.period = 80, plot.waves = TRUE, lwd = c(1,2),
            legend.coords = "bottomleft")

# reconstruct 2, use visually picked periods= 30 and 80
reconstruct(my.w, sel.period = c(30, 80), plot.waves = TRUE, lwd = c(1,2),
                   legend.coords = "bottomleft")

# reconstruct 3, use all relevant periods
reconstruct(my.w, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

# why period =80, but actually on plot it shows 79.9
my.w$Period[(my.w$Period > 79) & (my.w$Period < 81)]



#  2.4 Compare 2 similar series ------------------------------------------------

x1 <- periodic.series(start.period = 100, length = 500)
x2 <- periodic.series(start.period = 60, length = 500)
x <- c(x1, x2) 

y1 <- periodic.series(start.period = 100, length = 1000)
y2 <- periodic.series(start.period = 60, length = 1000)
y <- (y1 + y2)/2

plot(x, type = "l")  
plot(y, type = "l")  

my.data <- data.frame(x = x, y = y)
my.wx <- analyze.wavelet(my.data, "x", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 16, upperPeriod = 256,
                         make.pval = TRUE, n.sim = 10)
my.wy <- analyze.wavelet(my.data, "y", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 16, upperPeriod = 256,
                         make.pval = TRUE, n.sim = 10)

wt.image(my.wx, n.levels = 250,
         legend.params = list(lab = "wavelet power levels") )
wt.image(my.wy, n.levels = 250,
         legend.params = list(lab = "wavelet power levels") )

reconstruct(my.wx, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")
reconstruct(my.wx, sel.period = c(60, 100), plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")
reconstruct(my.wy, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")
reconstruct(my.wy, sel.period = c(60, 100), plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

# plot in grey scale
wt.image(my.wy, n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         color.palette = "gray((n.levels):1/n.levels)",
         col.ridge = "blue")



# 3 Analysis of a bivariate time series-----------------------------------------

#  3.1 bivariate  constant periods ---------------------------------------------

x1 <- periodic.series(start.period = 1*24, length = 24*40)
x2 <- periodic.series(start.period = 2*24, length = 24*40)
x3 <- periodic.series(start.period = 4*24, length = 24*40)
x4 <- periodic.series(start.period = 8*24, length = 24*40)
x5 <- periodic.series(start.period = 16*24, length = 24*40)
x <- x1 + x2 + 3*x3 + x4 + x5 
y <- x1 + x2 - 3*x3 + x4 + 3*x5

plot(x, type = "l", col = "red")  
lines(y, type = "l", col = "blue")  

my.data <- data.frame(x = x, y = y)
my.wc <- analyze.coherency(my.data, my.pair = c("x","y"),
                           loess.span = 0,
                           dt = 1/24, dj = 1/100,
                           lowerPeriod = 1/2,
                           make.pval = TRUE, n.sim = 10)

# plot individual wavelet analysis results.  
wt.image(my.wc, my.series = "x")
wt.image(my.wc, my.series = "y")

# Plot cross-wavelet analysis of x and y
# By taking away 'plot.ridge = TRUE,' no ridge will be ploted, because the ridge overlaps with arrows
wc.image(my.wc, n.levels = 250, plot.ridge = TRUE,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period (days)")



# 4 An example of importing my own data, then wavelet analysis, plot, reconstruct---------------------
library(xts)

# import daily data, make monthly mean data
library(readr)
DailyIm <- read_csv("05UE005_Daily_Flow_ts.csv")  # I download streamflow from CanadaEnvironment

date.lookup <- format(seq(as.Date("1960-07-01"), as.Date("2014-12-31"), by = "1 day"))

MyDayDF <- DailyIm[,4]

myXtsDay <- xts(MyDayDF, order.by=as.Date(date.lookup)) 
myMonthXts <- apply.monthly(myXtsDay, mean)    # mean or sum
MyMonthDF <- as.data.frame(myMonthXts)

# change column names. These 2 are going to be used in Wavelet analysis
colnames(MyDayDF) <- "x"
colnames(MyMonthDF) <- "y"

# plot daily and monthly data
plot(myXtsDay, col = "blue")
plot(MyMonthDF[,"y"], type = "l", col = "red")

# Replace NA with 0 (you may use other values)
match("TRUE", is.na(MyDayDF))
MyDayDF[is.na(MyDayDF)] <- 0                  # replace NA with 0, or other value
match("TRUE", is.na(MyDayDF))

match("TRUE", is.na(MyMonthDF))
MyMonthDF[is.na(MyMonthDF)] <- 0
match("TRUE", is.na(MyMonthDF))

# Set up upperlimit using log2 in R
# if your target testing period is 3 years, for daily data
2^ceiling(log2(12*20))



# wavelet analysis cmparing daily and monthly
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



# 5 Extra. Package 'biwavelet' -------------------------------------------------
library(biwavelet)

MOTHERS    
#   "morlet" "paul"   "dog"   







