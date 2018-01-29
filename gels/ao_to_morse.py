from __future__ import print_function
import numpy as np
def reduced_second_virial(B2, sigma_eff):
	return B2/(2./3.*np.pi*sigma_eff**2)

def second_virial(potential,params, beta,rmax=20.0, npoints=10000):
	from scipy.integrate import simps
	r = np.linspace(0,rmax,npoints)
	y = r**2*(1.-np.exp(-beta*potential(r,params=params)))
	I = simps(y,r)*2*np.pi #don't forget the 2pi!
	return I

def effective_sigma(potential,params, beta,rmax=20.0, npoints=10000):
	from scipy.integrate import simps
	r = np.linspace(0,rmax,npoints)
	p = potential(r, params=params)
	repulsive  = np.zeros(len(r))
	repulsive [p>0] = p[p>0]

	y = 1.-np.exp(-beta*repulsive)
	I = simps(y,r)
	return I

def morse(r,params=[1.,1.,1.],):
	rho ,sigma, epsilon= params[0],params[1], params[2]
	p= epsilon*np.exp( rho*(sigma-r) )*(np.exp(rho*(sigma-r))-2.)
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
	# print ("Rg", Rg)
	zp = 6.*etap/(np.pi*(2*Rg)**3)

	rhard = r<= sigma
	rgood = (r>sigma)*(r< (sigma+2*Rg))
	
	p = np.zeros(len(r))

	p[rhard] = np.inf
	p[rgood] = -np.pi*(2*Rg)**3*zp/6.*(1+q)**3/q**3*( 1.- 3*r[rgood]/(2*(1.+q)*sigma)+(r[rgood])**3/(2.*(1+q)**3*sigma**3) )

	return p
def ao_contact( etap ,q):
	return etap*(1.+3./2./q)
def etap_from_ao_contact(u_ao,q):
	return u_ao/(1.+3./2./q)

def find_epsilon(epsilons, morse_reduced_b2s, target_reduced_b2):
	pos = abs(np.array(morse_reduced_b2s)-target_reduced_b2).argmin()
	return epsilons[pos], pos

def square_well(r, params):
	sigma = params[0]
	epsilon = params[1]
	lmbd = params[2]
	p = np.zeros(len(r))
	p [r<(sigma*lmbd)] = -epsilon
	p [r<=sigma]=np.inf
	return p

def test_square_well():
	print ("\n\nValidation with Square Well at Criticality:\nExpected tau = -1/(4 (B2_star -1)) = 0.0765")
	# Validation of the calculation:
	# check Noro and Frenkel J. Chem. Phys., Vol. 113, No. 8, 22 August 2000
	# table 1, for the Square  Well potential

	b2_sw = second_virial(square_well,[1,1,2],1/2.61)
	sigma_eff =  effective_sigma(square_well, [1,1,2],1/2.61)
	b2_star= reduced_second_virial(b2_sw, 1)
	tau = 1./(-4*(b2_star-1))
	print ("\n*** Found : B2_star", b2_star,"tau",tau)

def optimise_morse_epsilon_from_Uao(rho0, sigma, q,Uao,numepsilons=10000):
	""" Get optimised values for Morse's epsilon assuming beta=1"""

	if type(Uao)!=list:
		Uao = [Uao]
	etap = [ etap_from_ao_contact(u,q) for u in Uao]

	ao_params = [[q,sigma,e] for e in etap]
	# compute the reduced B2 for morse for a range of temperatures

	reduced_B2s_morse = []
	epsilons = np.linspace(min(Uao)/2.,max(Uao)*2,numepsilons)
	beta = 1.0

	print ("\n ... Starting Mapping ...\n")
	for epsilon in epsilons:	
		morse_params = [rho0,sigma,epsilon]
		B2_morse = second_virial(morse_cut,morse_params,beta)
		effective_sigma_morse = effective_sigma(morse_cut, morse_params,beta)
		reduced_B2_morse = reduced_second_virial(B2_morse, effective_sigma_morse)
		reduced_B2s_morse.append(reduced_B2_morse)

	optimised_epsilons = []

	for ao_p in ao_params:
		print ("====> Uao ",ao_contact(ao_p[-1],q))
		# Compute the reduced B2 for AO
		B2_ao = second_virial(ao_potential,ao_p,beta)
		effective_sigma_ao= effective_sigma(ao_potential, ao_p,beta)
		reduced_B2_ao = reduced_second_virial(B2_ao, effective_sigma_ao)

		eps ,pos= find_epsilon(epsilons,reduced_B2s_morse, reduced_B2_ao)
		optimised_epsilons.append(eps)

		print ("*** Morse BetaEpsilon", beta*eps)
		print ("")
		print ("AO reduced B2", reduced_B2_ao)
		print ("Morse reduced B2",reduced_B2s_morse[pos] )
		print ("Morse effective sigma ",effective_sigma(morse_cut, morse_params,beta) )
		r = np.linspace(0.5,1.5,1000)
		print ("Min AO", min(ao_potential(r,ao_p)))
		print ("")
	return optimised_epsilons

def example():
	# INPUT : 
	q = 0.17
	rho0= 25.
	sigma = 1.
	Uao = [2.4, 4.] #well depths in kbT
	optimise_morse_epsilon_from_Uao(rho0, sigma, q, Uao)
