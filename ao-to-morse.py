import pylab as pl
import numpy as np
def reduced_second_virial(B2, sigma_eff):
	return B2/(2./3.*pl.pi*sigma_eff**2)

def second_virial(potential,params, beta,rmax=100, npoints=1e5):
	from scipy.integrate import cumtrapz
	r = pl.linspace(0,rmax,npoints)
	y = r**2*(1.-pl.exp(-beta*potential(r,params=params)))
	I = cumtrapz(y,r,initial=0)
	return I[-1]

def effective_sigma(potential,params, beta,rmax=100, npoints=1e5):
	from scipy.integrate import cumtrapz
	r = np.linspace(0,rmax,npoints)
	p = potential(r, params=params)
	repulsive  = pl.zeros(len(r))
	repulsive [p>0] = p[p>0]

	y = 1.-pl.exp(-beta*repulsive)
	I = cumtrapz(y,r,initial=0)
	return I[-1]

def morse(r,params=[1.,1.,1.],):
	rho ,sigma, epsilon= params[0],params[1], params[2]
	p= epsilon*pl.exp( rho*(sigma-r) )*(pl.exp(rho*(sigma-r))-2.)
	return p

def morse_cut(r,params=[1.,1.,1.],rcut=2.5):
	p = morse(r,params)
	pcut=morse(rcut,params)
	pp = p -pcut
	pp[r>rcut]=0
	return pp

def ao_potential(r, params=[1.,1.,1]):
	q, sigma,etap = params[0], params[1],params[2]

	Rg = sigma*q/2.
	zp = 6.*etap/(pl.pi*(2*Rg)**3)

	rhard = r<= sigma
	rgood = (r>sigma)*(r< (sigma+2*Rg))
	
	p = pl.zeros(len(r))

	p[rhard] = np.inf
	p[rgood] = -pl.pi*(2*Rg)**3*zp/6.*(1+q)**3/q**3*( 1.- 3*r[rgood]/(2*(1.+q)*sigma)+(r[rgood])**3/(2.*(1+q)**3*sigma**3) )

	return p
def ao_contact( etap ,q):
	return etap*(1.+3./2./q)
def etap_from_ao_contact(u_ao,q):
	return u_ao/(1.+3./2./q)

def find_beta(betas, morse_reduced_b2s, target_reduced_b2):
	pos = abs(pl.array(morse_reduced_b2s)-target_reduced_b2).argmin()
	return betas[pos], pos
# INPUT : 

q = 0.17
rho0= 25.
sigma = 1

epsilon = 1

Uao = [2.4, 4.] #well depths in kbT

etap = [ etap_from_ao_contact(u,q) for u in Uao]

ao_params = [[q,sigma,e] for e in etap]

morse_params = [rho0,sigma,epsilon]
# compute the reduced B2 for morse for a range of temperatures

reduced_B2s_morse = []
betas = pl.linspace(0.5,4,1000)
for beta in betas:	

	B2_morse = second_virial(morse_cut,morse_params,beta)
	effective_sigma_morse = effective_sigma(morse_cut, morse_params,beta)
	reduced_B2_morse = reduced_second_virial(B2_morse, effective_sigma_morse)
	reduced_B2s_morse.append(reduced_B2_morse)

for ao_p in ao_params:
	print ("====> Uao ",ao_contact(ao_p[-1],q))
	# Compute the reduced B2 for AO
	B2_ao = second_virial(ao_potential,ao_p,1)
	effective_sigma_ao= effective_sigma(ao_potential, ao_p,1)
	reduced_B2_ao = reduced_second_virial(B2_ao, effective_sigma_ao)

	beta ,pos= find_beta(betas,reduced_B2s_morse, reduced_B2_ao)
	print ("Morse BetaEpsilon", beta*epsilon, "T", 1./beta)
	print ("AO reduced B2", reduced_B2_ao)
	print ("Morse reduced B2",reduced_B2s_morse[pos] )
	print ("Morse effective sigma ",effective_sigma(morse_cut, morse_params,beta) )

	# pl.plot(betas, reduced_B2_ao*np.ones(len(betas)), '-k')
	# pl.plot(betas, reduced_B2s_morse, '-')
	# pl.show()	

	# plot the last potential
	r = pl.linspace(0.5,1.5,1000)
	pl.plot(r, ao_potential(r,ao_p))
	pl.plot(r,beta*morse_cut(r, morse_params))
	pl.ylim(np.floor(ao_potential(r,ao_params[1]).min())-1,10)
	pl.show()