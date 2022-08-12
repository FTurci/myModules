
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


def next_pow_two(n):
	"""Bitwise shift operations to find the closest power of 2"""
	i = 1
	while i < n:
		i = i << 1
	return i

def autocorr_func_1d(x, norm=True):
	"Calculate the autocorrelation of a time signal. Code from emcee"
	x = np.atleast_1d(x)
	if len(x.shape) != 1:
		raise ValueError("invalid dimensions for 1D autocorrelation function")
	n = next_pow_two(len(x))

	# Compute the FFT and then (from that) the auto-correlation function
	f = np.fft.fft(x - np.mean(x), n=2 * n)
	acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
	acf /= 4 * n

	# Optionally normalize
	if norm:
		acf /= acf[0]

	return acf

def auto_window(taus, c):
	# Automated windowing procedure following Sokal (1989)
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1

def autocorr_time(x,c=5):
	"""Calculate the integrated autocorrelation time of a series. The input is just a 1d array (time is assumed to be discretized in equal steps of size 1)"""
	x = np.atleast_1d(x)
	if len(x.shape) != 1:
		raise ValueError("invalid dimensions for 1D autocorrelation function")
	acf = autocorr_func_1d(x)
	taus = 2*np.cumsum(acf) - 1.0
	return taus[auto_window(taus,c)]
