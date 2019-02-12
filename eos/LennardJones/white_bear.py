
# Twitter @Francesco_Turci

import numpy as np
import scipy.integrate as integrate
from scipy.optimize import fsolve
import matplotlib.pyplot as pl
import warnings
from scipy.interpolate import UnivariateSpline


def volume_fraction(rho, sigma):
		return np.pi/6*sigma**3*rho

class WhiteBearEOS(object):
	"""Computing the Lennard-Jones chemical potential, pressure and binodal within the Carnahan-Starling approximation."""
	def __init__(self, sigma, epsilon,rcut,infinity=100):
		super(WhiteBearEOS, self).__init__()
		self.sigma = sigma
		self.epsilon = epsilon
		self.rcut= rcut
		self.infinity = infinity*rcut

		self.integral_att = self.get_integral_att()
	
	def Vatt(self,r,rcut):
		"""Attractive part of the potential"""
		rmin = 2**(1./6.)*self.sigma 
		value = 4 *self.epsilon*( (self.sigma/r)**12-(self.sigma/r)**6)
		if r<rmin :
			return  - self.epsilon
		elif r>rcut:
			return 0
		return value


	def mu_CS(self,rho,T):
		"""Carnahan-Starling approximation for the hard-sphere chemical potential."""
		eta = volume_fraction(rho,self.sigma)
		return T*(np.log(rho)+(8*eta-9*eta**2+3*eta**3)/(1-eta)**3 )

	def p_CS(self,rho,T):
		"""Carnahan-Starling approximation for the hard-sphere pressure."""
		eta = volume_fraction(rho,self.sigma)
		return T*rho*(1+eta+eta**2-eta**3)/(1-eta)**3

	def get_mu(self,rho,T):
		"""Hard-sphere chemical potential plus the integral over the attractive contribution."""
		return self.mu_CS(rho,T)+4*np.pi*rho*self.integral_att

	def get_p(self, rho,T):
		"""Hard-sphere pressure plus the integral over the attractive contribution."""
		return self.p_CS(rho,T)+2*np.pi*rho**2*self.integral_att

	def get_integral_att(self):
		"""Integral over the attractive contribution of the  pair potential."""
		def integrand(r):
			return r**2*self.Vatt(r,self.rcut)
		#resolve the small r and large r regions differently
		integration_points = np.concatenate((np.linspace(0, 10,1000), np.linspace(11, self.infinity, 10000)))

		# integrate via quadrature, using QUADPACK
		value = integrate.quad(integrand,0,self.infinity, limit=len(integration_points)*2, points=integration_points)

		result = value[0]
		precision = value[1]
		print (f"::: For cutoff {self.rcut} sigma the attractive integral has value {result} with precision {precision}.")
		return result

	def find_coex(self,T,guesslo,guesshi, maxfev):
		"""Solve simultaneous equations for coexistence:

			mu(liquid) = mu(vapor)
			p (liquid) = p (vapor)

		"""
		def equations(p):
			""" Inline simultaneous equations"""
			rhov,rhol = p
			return (self.get_mu(rhov,T)-self.get_mu(rhol,T), self.get_p(rhov,T)-self.get_p(rhol,T))
		# solve the simultaneous nonlinear equations with MINPACk
		# via Powell hybrid method
		x,y = fsolve(equations,(guesslo,guesshi),maxfev=maxfev)
		return x,y

	def coexistence(self, Thigh, Tlow=0.694, npoints=1000, maxfev=10000, inital_guess_lo=0.001, inital_guess_hi=0.87):
		"""Find the binodal coexistence densities between temperature Tlow and Thigh. 

		The computation starts at low temperature with given initial guesses and climbs up to the highest temperature.

		If the temperature is too high (i.e. beyond the critical point) warning will be issued by the equation of state methods.
		"""
		Ts = np.linspace(Tlow, Thigh , npoints )
		los, his =[] ,[]
		for i,T in enumerate(Ts):
			if i==0:
				x,y = self.find_coex(T,inital_guess_lo,inital_guess_hi,maxfev)
			else:
				x,y = self.find_coex(T,x,y,maxfev)
			los.append(x)
			his.append(y)
		result={}
		result['rho_vapor']=np.array(los)
		result['rho_liquid']=np.array(his)
		result['temperature']=Ts
		self.binodal = result
		return result

	def plot_binodal(self, show=True, color="#0096ff"):
		"""Plot the binodal with dashed and continuous coloured lines."""
		pl.plot(self.binodal['rho_vapor'],self.binodal['temperature'], '--', color=color)
		pl.plot( self.binodal['rho_liquid'],self.binodal['temperature'],'-', color=color)
		pl.xlabel(r"$\rho^*$")
		pl.ylabel(r"$T^*$")
		if show:
			pl.show()

	def get_binodal_densities(self,T):
		"""Get the coexistence densities for a given temperature."""
		self.coexistence(T,npoints=100)
		return self.binodal['rho_vapor'][-1], self.binodal['rho_liquid'][-1]
	def get_binodal_temperature(self,rho,first_max_temperature=1.3,low_temperature=0.6,epsilon=1e-8):
		"""Get the temperature at which the given density crosses the binodal. 

		First the binodal is estimated in a wide range of temperatures. Then, a spline is used to fit the relation T(rho) and return the estimate for T at the given rho.
		"""
		T = first_max_temperature
		
		while True:
			T *= 1.05 #5 percent increment to quickly climb  up
			try:
				self.coexistence(Thigh=T,Tlow=low_temperature, npoints=1000)
			except Exception as e:
				break

		diff = np.abs(np.gradient(self.binodal['rho_vapor']))
		self.binodal['rho_vapor'] = self.binodal['rho_vapor'][diff>epsilon]
		self.binodal['rho_liquid'] = self.binodal['rho_liquid'][diff>epsilon]
		self.binodal['temperature'] = self.binodal['temperature'][diff>epsilon]

		if rho < max(self.binodal['rho_vapor']):
			spline = UnivariateSpline(self.binodal['rho_vapor'],self.binodal['temperature'] ,s=0)
		else:
			order = np.argsort(self.binodal['rho_liquid'])
			spline = UnivariateSpline(self.binodal['rho_liquid'][order],self.binodal['temperature'][order],s=0 )
		return spline(rho)+0.0


