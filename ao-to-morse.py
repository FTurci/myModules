import pylab as pl
import numpy as np
def reduced_second_virial(B2, sigma_eff):
    return B2/(2./3.*pl.pi*sigma_eff**2)

def second_virial(potential,params, beta,rmax=10, npoints=1e4):
    from scipy.integrate import cumtrapz
    r = pl.linspace(0,rmax,npoints)
    y = r**2*(1.-pl.exp(-beta*potential(r,params=params)))
    I = cumtrapz(y,r,initial=0)
    return I[-1]

def effective_sigma(potential,params, beta,rmax=10, npoints=1e4):
    from scipy.integrate import cumtrapz
    r = np.linspace(0,rmax,npoints)
    p = potential(r, params=params)
    repulsive  = pl.zeros(len(r))
    repulsive [p>0] = p[p>0]

    y = 1.-pl.exp(-beta*repulsive)
    I = cumtrapz(y,r,initial=0)
    return I[-1]

def morse(r,params=[1.,1.,1.]):
    rho ,sigma, epsilon= params[0],params[1], params[2]
    return epsilon*pl.exp( rho*(sigma-r) )*(pl.exp(rho*(sigma-r))-2.)

def ao_potential(r, params=[1.,1.,1]):
    q, sigma,etap = params[0], params[1],params[2]

    Rg = sigma*q/2.
    zp = 6.*etap/(pl.pi*(2*Rg)**3)

    rhard = r<= sigma
    rgood = (r>sigma)*(r< (sigma+2*Rg))
    
    p = pl.zeros(len(r))

    p[rhard] = np.nan
    p[rgood] = -pl.pi*(2*Rg)**3*zp/6.*(1+q)**3/q**3*( 1.- 3*r[rgood]/(2*(1.+q)*sigma)+(r[rgood])**3/(2.*(1+q)**3*sigma**3) )

    return p

print ( effective_sigma(morse,[25,1,3],1)) 
r = pl.linspace(0.9,1.8,1000)
pl.plot(r,morse(r,[25,1,3]))
pl.plot(r,ao_potential(r,[0.1,1,0.1]))
pl.ylim(-2,0.5)
pl.show()


