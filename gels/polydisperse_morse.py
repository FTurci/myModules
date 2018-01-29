# October 2015
# 
# Writing a LAMMPS input file for a system with Morse interaction and gaussian polydispersity
import numpy as np
from scipy.stats import norm
import sys
from __future__ import print_function

def match_packing_fraction(input_packing_fraction,mu, polydispersity,N):
    """Find cubic box size for a system of N particles of mean size mu atgiven polydispersity"""
    sigma = polydispersity*mu
    print ("    * The imposed polydispersity is", polydispersity)
    # 7 species
    Ds=np.arange(mu-3*sigma, mu+3*sigma+sigma, sigma)-sigma/2.
    diams=Ds[:-1]+(Ds[1]-Ds[0])/2
    label_types=np.arange(len(diams)+1)+1
    cumulative=norm.cdf(Ds,mu,sigma)
    cumulative[0]=0
    cumulative[-1]=1
    Types=[]
    for p in range(N):
        random_num=np.random.uniform(0,1)
        Test=cumulative<random_num
        Types.append(label_types[Test][-1])
    Num_Types=len(Ds)-1

    def packing(scale):
        # cubic box side:
        L=scale*mu
        # the volume is
        V=L**3
        packing=0
        for i in range(N):
            packing+=np.pi/6*(diams[Types[i]-1]**3)/V
        return packing

    def delta (scale):
        return abs(  ( (packing(scale)-input_packing_fraction )/input_packing_fraction ))

    from scipy.optimize import minimize_scalar
    res = minimize_scalar(delta,bracket= [2,1000])
    L = res.x*mu
    print ("    * Box size L =", L)
    print ("    * The packing fraction is", packing(res.x))
    return L,Types,Num_Types,diams

def generate_input_file(input_packing_fraction,epsilon, rho0,mu, polydispersity,N,r_cut_coeff,filename="morse_input_phi"):
    L ,Types,Num_Types, diams= match_packing_fraction(input_packing_fraction,mu, polydispersity,N)
    packing = input_packing_fraction

    V = L**3
    # random coordinates
    x=np.random.uniform(0,L,size=N)
    y=np.random.uniform(0,L,size=N)
    z=np.random.uniform(0,L,size=N)


    print ("    * The number density is",N/float(V))

    
    print ("    * Writing file",filename,"...")

    if filename == "morse_input_phi":
        filename = filename+"%g.lmp"%packing
    with open(filename,'w') as fw:
        # write the header
        fw.write("LAMMPS Description\n\n")
        fw.write("%d atoms\n\n"%N)
        # every atom is a type
        fw.write("%d atom types\n"%Num_Types)
        fw.write("""
    0 %g xlo xhi
    0 %g ylo yhi
    0 %g zlo zhi
    \n"""%(L,L,L))
        # all the atoms have different masses, but equal mass density 1
        fw.write("Masses\n\n")
        for i in range(Num_Types):
            fw.write("%d %g\n"%(i+1,np.pi/6*diams[i]**3))

        fw.write("\n")
        # specify the interaction coefficient for the Morse potential:
        # in lammps's documentation they are:
        # d0 alpha r0 cutoff
        # in Paddy's notation:
        # epsilon rho0 sigma cutoff
        fw.write("PairIJ Coeffs # morse\n\n")
        for i in range(Num_Types):
            # fw.write("%d %g %g %g \n"%(i+1,epsilon,rho0,diam[i]))#, r_cut_coeff*diam[i] ))
            for j in range(i, Num_Types):
                diam_i=diams[i]
                diam_j=diams[j]
                mixed_diam=0.5*(diam_i+diam_j)
                fw.write("%d %d %g %g %g %g\n"%(i+1,j+1,epsilon,rho0,mixed_diam, r_cut_coeff*mixed_diam ))
        
        fw.write("\nAtoms\n\n")
        # write random coordinates
        for i in range(N):
            fw.write("%d %d %g %g %g\n"%(i+1,Types[i], x[i],z[i],y[i]))

    print ("    * Compressing...",end='')
    import os
    os.system("gzip "+filename)
    print ("    ...done.")

def example(input_packing_fraction = 0.3,epsilon=1,rho0=25,mu=1,polydispersity=0.04,N=10000, r_cut_coeff=2.5):
    generate_input_file(input_packing_fraction,epsilon, rho0,mu, polydispersity,N,r_cut_coeff)