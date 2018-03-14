from __future__ import division, print_function
import numpy as np
import ctypes
import tqdm 
import pylab as pl

mollib = ctypes.cdll['./mollified.so']
c_double_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int32)

def radial_distribution(x,y,z,box, nbins):
	"""
	Compute the radial distribution fucntion.
	"""
	# deal with the pointers
	x = x.astype(np.float64)
	y = y.astype(np.float64)
	z = z.astype(np.float64)
	box = np.array(box).astype(np.float64)
	# output vectors
	gr = np.zeros(nbins).astype(np.float64)
	r = np.zeros(nbins).astype(np.float64)

	N = len(x)
	
	compute_g = mollib['radial_distribution']

	compute_g(
		gr.ctypes.data_as(c_double_p),
		r.ctypes.data_as(c_double_p),
		x.ctypes.data_as(c_double_p),
		y.ctypes.data_as(c_double_p),
		z.ctypes.data_as(c_double_p),
		box.ctypes.data_as(c_double_p),
		nbins,N)

	return r,gr

def mollified_radial_distribution(index,x,y,z,box,stdev,rcut=None, nbins=200):
	"""
	Compute mollified radial distributioon function for particle i=index.
	- stdev: sigma of the gaussian kernels
	- rcut: maximum distance for the radial distribution
	"""
	# deal with the pointers

	x = x.astype(np.float64)
	y = y.astype(np.float64)
	z = z.astype(np.float64)
	N = len(x)

	if rcut ==None:
		rcut = 300 *stdev
	
	box = np.array(box).astype(np.float64)
	int_params = np.array([index,N,nbins]).astype(np.int32)
	double_params = np.array([stdev,rcut]).astype(np.float64)
	# print (int_params)
	# output vector
	gr = np.zeros(nbins).astype(np.float64)
	r = np.zeros(nbins).astype(np.float64)
	
	
	compute_mollified_g = mollib['mollified_radial_distribution']

	compute_mollified_g(
		gr.ctypes.data_as(c_double_p),
		r.ctypes.data_as(c_double_p),
		x.ctypes.data_as(c_double_p),
		y.ctypes.data_as(c_double_p),
		z.ctypes.data_as(c_double_p),
		box.ctypes.data_as(c_double_p),
		int_params.ctypes.data_as(c_int_p),
		double_params.ctypes.data_as(c_double_p))

	return r,gr


def s2(x,y,z,box,stdev,rcut=None, nbins=200):
	"""
	Compute directly s2 and a table of mollified gr for every particle.
	"""
	x = x.astype(np.float64)
	y = y.astype(np.float64)
	z = z.astype(np.float64)
	N = len(x)

	if rcut ==None:
		rcut = 300 *stdev
	
	box = np.array(box).astype(np.float64)
	int_params = np.array([N,nbins]).astype(np.int32)
	double_params = np.array([stdev,rcut]).astype(np.float64)

	_s2 = np.zeros(N).astype(np.float64)
	gr = np.zeros(N*nbins).astype(np.float64)
	r = np.zeros(nbins).astype(np.float64)

	compute_s2 = mollib['s2']
	# print ("here")
	compute_s2(
		_s2.ctypes.data_as(c_double_p),
		r.ctypes.data_as(c_double_p),
		gr.ctypes.data_as(c_double_p),
		x.ctypes.data_as(c_double_p),
		y.ctypes.data_as(c_double_p),
		z.ctypes.data_as(c_double_p),
		box.ctypes.data_as(c_double_p),
		int_params.ctypes.data_as(c_int_p),
		double_params.ctypes.data_as(c_double_p)
		)
	return r,gr,_s2



def local_s2(x,y,z,box,stdev,neighcut,rcut=None,nbins=200,epsilon=1e-10):
	"""
	Compute s2 and its local average within a sphere of size neighcut.
	"""

	from scipy.integrate import simps
	if rcut ==None:
		rcut = 300 *stdev

	x = x.astype(np.float64)
	y = y.astype(np.float64)
	z = z.astype(np.float64)
	N = len(x)
	rho = N*1./np.prod(box)
	box = np.array(box).astype(np.float64)

	s2 = np.ones(N)
	locals2 = np.zeros(N)
	for i in tqdm.tqdm(range(int(N))):
		ri, gi = mollified_radial_distribution(i, x, y, z, box, stdev,rcut=rcut,nbins=nbins)
		# print (gi.max(),gi.min())
		# gi =np.ones(len(gi))
		# print
		integrand = (gi*np.log(gi)-gi+1.0)*ri**2
		
		# print integrand
		integrand[gi<epsilon]=0
		if i<3:
			pl.plot(integrand)
		# s2[i] = -2*np.pi *rho*simps(integrand,ri)
		s2[i] = -2.*np.pi *rho*np.trapz(integrand,ri)
		# print np.trapz(integrand,ri
	# compute the local average
	pl.show()
	compute_local_average = mollib['local_average']

	compute_local_average(
		locals2.ctypes.data_as(c_double_p),
		s2.ctypes.data_as(c_double_p),
		x.ctypes.data_as(c_double_p),
		y.ctypes.data_as(c_double_p),
		z.ctypes.data_as(c_double_p),
		box.ctypes.data_as(c_double_p),
		N, ctypes.c_double(neighcut))
	return s2, locals2

def save_xyzsl(filename,x,y,z,s,slocal):
	"""
	Store the results in an extended XYZ file. The last two columns are s2 and the locally averaged value of s2.
	"""
	with open(filename, 'w') as fw:
		fw.write("%d\nAtoms\n"%len(x))
		for i in range(len(x)):
			fw.write("A %g %g %g %g %g\n"%(x[i],y[i],z[i], s[i], slocal[i]))

# EXAMPLE:
# Choose the parameters well:
# the number of bins is critical: the more the better the integral... (but this slows the code down a lot)
# the sigma fo the mollifying gaussians needs to be small, around 5-10% of the diameter
# There are two cutoffs: 
# - rcut  defines the range of the mollified gr 
# - neighcut is the maximum distance for two neighboring particles

# LJFCC
box = [13.1486,13.1486,13.1486] 
x,y,z=np.loadtxt("fcc.xyz", skiprows=2, usecols=[1,2,3], unpack=True)
diameter =1.
s2,locals2 = local_s2 (x, y, z, box, 0.1*diameter, 1.4*diameter,rcut=3*diameter,nbins=60) 
print ("FCC",np.mean(s2),np.mean(locals2))
# LJ liquid
x,y,z=np.loadtxt("liquid.xyz", skiprows=2, usecols=[1,2,3], unpack=True)
box = [13.4061,13.4061,13.4061] 
diameter =1.
s2,locals2 = local_s2 (x, y, z, box, 0.1*diameter, 1.4*diameter,rcut=3*diameter,nbins=60) 
print ("LIQ",np.mean(s2),np.mean(locals2))
# LJ liquid
x,y,z=np.loadtxt("centers020same_kernel_6_octaves_blur2_no_overlap_-0_radfix_14_8_3_16.xyz", skiprows=2, usecols=[1,2,3], unpack=True)
box = [340.,340.,340.] 
diameter =19.
s2,locals2 = local_s2 (x, y, z, box, 0.1*diameter, 1.4*diameter,rcut=3*diameter,nbins=60) 
print ("EXP",np.mean(s2),np.mean(locals2))


