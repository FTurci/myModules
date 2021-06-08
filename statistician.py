
# Nonlinear curve fit with confidence interval
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t
from scipy.stats import binned_statistic

def confidenceFit(func,x,y, initial_guess, alpha=0.1,maxfev=800):

	# alpha is the degree of confidence interval = 100*(1-alpha), for example for 98% take alpha=0.02
	pars, pcov = curve_fit(func, x, y, p0=initial_guess,maxfev=maxfev)



	n = len(y)    # number of data points
	p = len(pars) # number of parameters

	dof = max(0, n - p) # number of degrees of freedom
	#print pcov
	# student-t value for the dof and confidence level
	tval = t.ppf(1.0-alpha/2., dof) 
	#std from covariance	
	sigma = np.diag(pcov)**0.5
	lower_bound = pars - sigma*tval
	upper_bound = pars + sigma*tval

	return pars,lower_bound,upper_bound,sigma,tval
	


def binned_cloud(x,y, bins=None):
	if bins==None:
		bins = np.unique(x)
		bins = np.append(bins, 2*bins[-1]-bins[-2])
	bs, be, bn= binned_statistic(x,y,bins= bins, statistic="mean")
	std, be, bn= binned_statistic(x,y,bins= bins, statistic="std")
	count, be, bn = binned_statistic(x,y,bins= bins, statistic="count")

	er = std/np.sqrt(count)
	return be[:-1],bs,er