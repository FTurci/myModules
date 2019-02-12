# Conversion from legacy FORTRAN code to Python
# Equation of state from J. Kolafa, I. Nezbeda, Fluid Phase Equil. 100 (1994), 1
# Twitter @Francesco_Turci

import numpy as np
from numpy import pi
import scipy.integrate as integrate
from scipy.optimize import fsolve
import matplotlib.pyplot as pl
import warnings
from scipy.interpolate import UnivariateSpline


############### Begin of legacy code ###################
def ALJ(rho,T):
	"""Helmholtz free energy (including the ideal term)"""
	eta = pi/6.*rho * (dC(T))**3
 
	return (np.log(rho)+betaAHS(eta) +rho*BC(T)/np.exp(gammaBH(T)*rho**2) )*T +DALJ(rho,T)


def ALJres(rho,T):
      """C/* Helmholtz free energy (without ideal term) */"""
      eta = pi/6. *rho*(dC(T))**3
      return (betaAHS(eta)+rho*BC(T)/np.exp(gammaBH(T)*rho**2))*T+DALJ(rho,T)


def PLJ(rho,T):
	"""Pressure"""
	eta=pi/6. *rho*(dC(T))**3
	sum=((2.01546797*2+rho*( (-28.17881636)*3+rho*( 28.28313847*4+rho* (-10.42402873)*5))) +((-19.58371655)*2+rho*( +75.62340289*3+rho*( (-120.70586598)*4+rho*( +93.92740328*5+rho* (-27.37737354)*6))))/np.sqrt(T) + ((29.34470520*2+rho*( (-112.35356937)*3+rho*( +170.64908980*4+rho*( (-123.06669187)*5+rho* 34.42288969*6))))+ ((-13.37031968)*2+rho*( 65.38059570*3+rho*( (-115.09233113)*4+rho*( 88.91973082*5+rho* (-25.62099890)*6))))/T)/T)*rho**2

	return ((zHS(eta)+ BC(T)/np.exp(gammaBH(T)*rho**2)*rho*(1-2*gammaBH(T)*rho**2))*T+sum )*rho

def ULJ( T, rho):
	"""Internal Energy"""
	dBHdT=dCdT(T)
	dB2BHdT=BCdT(T)
	d=dC(T)
	eta=pi/6. *rho*d**3
	sum= ((2.01546797+rho*( (-28.17881636)+rho*( +28.28313847+rho* (-10.42402873)))) + (-19.58371655*1.5+rho*( 75.62340289*1.5+rho*( (-120.70586598)*1.5+rho*( 93.92740328*1.5+rho* (-27.37737354)*1.5))))/np.sqrt(T) + ((29.34470520*2+rho*( -112.35356937*2+rho*(  170.64908980*2+rho*( -123.06669187*2+rho* 34.42288969*2)))) + (-13.37031968*3+rho*(  65.38059570*3+rho*(  -115.09233113*3+rho*( 88.91973082*3+rho* (-25.62099890)*3))))/T)/T) *rho*rho
	
	return 3*(zHS(eta)-1)*dBHdT/d+rho*dB2BHdT/np.exp(gammaBH(T)*rho**2) +sum

def zHS(eta):
	return (1+eta*(1+eta*(1-eta/1.5*(1+eta)))) / (1-eta)**3

def betaAHS( eta ):
	return np.log(1-eta)/0.6+ eta*( (4.0/6*eta-33.0/6)*eta+34.0/6 ) /(1.-eta)**2

def dLJ(T):     
	isT=1/np.sqrt(T)
	return ((( 0.011117524191338 *isT-0.076383859168060)*isT)*isT+0.000693129033539)/isT+1.080142247540047+0.127841935018828*np.log(isT)

def dC(T):
	sT=np.sqrt(T)
	return -0.063920968*np.log(T)+0.011117524/T-0.076383859/sT+1.080142248+0.000693129*sT

def dCdT( T):
	sT=np.sqrt(T)
	return   0.063920968*T+0.011117524+(-0.5*0.076383859-0.5*0.000693129*T)*sT

def BC( T):
	isT=1/np.sqrt(T)
	return  (((((-0.58544978*isT+0.43102052)*isT+.87361369)*isT-4.13749995)*isT+2.90616279)*isT-7.02181962)/T+0.02459877

def BCdT( T):
	isT=1/np.sqrt(T)
	return ((((-0.58544978*3.5*isT+0.43102052*3)*isT +0.87361369*2.5)*isT-4.13749995*2)*isT+2.90616279*1.5)*isT-7.02181962

def gammaBH(X):
	return 1.92907278

def DALJ(rho,T):
      return((+2.01546797+rho*(-28.17881636+rho*(+28.28313847+rho*(-10.42402873))))+(-19.58371655+rho*(75.62340289+rho*((-120.70586598)+rho*(93.92740328+rho*(-27.37737354)))))/np.sqrt(T)+ ( (29.34470520+rho*((-112.35356937)+rho*(+170.64908980+rho*((-123.06669187)+rho*34.42288969))))+(-13.37031968+rho*(65.38059570+rho*((-115.09233113)+rho*(88.91973082+rho* (-25.62099890)))))/T)/T) *rho*rho

def muLJ(rho,T):
	return ALJ(rho,T)+PLJ(rho,T)/rho-T

############### End of legacy code ###################

class KolafaNezbevdaEOS(object):
	"""Computing the Lennard-Jones chemical potential, pressure and binodal within the Kolafa-Nezbevda approximation."""
	def __init__(self):
		super(KolafaNezbevdaEOS, self).__init__()
	

	def get_mu(self,rho,T):
		"""Chemical potential."""
		return muLJ(rho,T)

	def get_p(self, rho,T):
		"""Pressure."""
		return PLJ(rho,T)


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

