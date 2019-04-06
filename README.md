# ATLAS_CheckSignificance
Script with functions to smooth data and determine the significance of a peak in a histogram.

testSmoother.cxx
=================
Sample script to run smoother.cxx


smoother.cxx
============

__TH1D* smooth(TH1D* data, int numSmooths)__ \
The “smoother.cxx” file contains a function (“smooth”) to smooth N data points a specified amount of times (numSmooths). 
The smoothing works by doing the following operation on each point:\
- x_(1 new)=(2*x_1+x_2)/3;             
- x_(i new)=(x_(i-1)+2*x_i+x_(i+1))/4;           
- x_(N new)=(x_(N-1)+2*x_N)/3  
 
__TH1D* significance(TH1D* data, TH1D* smoothed_data, int mergeBins)__ \
It then calculates the significance at each point with another function (“significance”) by the following formula:
- significance=(x_(i new)-x_(i old))/sqrt(x_(i old))

It also merges the final histogram by the specified number of bins (mergeBins). Note that the bin values cannot be empty.

__void significanceTest(TH1D* hist, int numSmooths, int mergeBins, const char* title, bool absVal=false,  bool drawLine=false, bool roundVal=false, bool printInfo=true)__ \
The function “significanceTest” can be used to smooth and obtain the significance plot of a histogram (hist). 
The number of smooths (numSmooths) and merged bins (mergeBins) needs to be specified, and the other parameters are:
- title = title of x-axis
- absVal = make significance the absolute value
- drawLine = plot curve for smoothed data instead of points
- roundVal = round the significance
- printInfo = show number of smooths and merges on plot


__void runMultipleSigTests(TH1D* hist, int lowSmooths, int highSmooths, int lowBins, int highBins, const char* title, bool absVal=false, bool round=false)__\
To test multiple numbers of smooths and merged bins on a single histogram, the function “runMultipleSigTests” can be used.
 
A range of merged bins (lowBins-highBins) and smooths (lowSmooths-highSmooths) must be specified, 
and there is once again the option to return the absolute value or rounded value of the significance.
 
A plot of each of the smoothed data and a significance plot for each combination of bins and smooths is plotted. It will also print to screen which combinations 
yielded significance greater than 3 (observation), or greater than 5 (discovery). Note that it will skip over numbers of merged bins if the number of original 
bins is indivisible by that value.


