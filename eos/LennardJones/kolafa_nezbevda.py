
from numpy import pi
import numpy as np
import pylab as pl
from scipy.optimize import fsolve

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


def find_coex(T,guesslo,guesshi):
	def equations(p):
		rhov,rhol = p
		return (muLJ(rhov,T)-muLJ(rhol,T), PLJ(rhov,T)-PLJ(rhol,T))
	x,y = fsolve(equations,(guesslo,guesshi),maxfev=10000)
	return x,y


def get_coex(T):
	assert T>0.5, "Only temperatures between 0.75 and 1.3 are accepted"
	Ts = np.arange(0.75,T,0.01)
	for i,T in enumerate(Ts):
		if i==0:
			x,y = find_coex(T,0.001,0.87)
		else:
			x,y = find_coex(T,x,y)
	return x,y


Ts = np.arange(0.75,1.35,0.01)


los, his =[] ,[]
for i,T in enumerate(Ts):
	if i==0:
		x,y = find_coex(T,0.001,0.87)
	else:
		x,y = find_coex(T,x,y)
	los.append(x)
	his.append(y)
	print (x,y	)

pl.plot(los,Ts,'o')
pl.plot(his,Ts,'o')
pl.xlabel("DENSITY")
pl.ylabel("Temperature")
pl.show()