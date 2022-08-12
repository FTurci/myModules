from ovito.io import import_file
import numpy as np
from ovito.modifiers import *
from stringato import extract_floats as ef
import glob
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import sys
import time
from ovito.data import BondsEnumerator

path = sys.argv[1]


pipeline = import_file(path, multiple_frames=True)
num_frames = pipeline.source.num_frames

rc = 1.0
pipeline.modifiers.append(CreateBondsModifier(cutoff = rc))

c = []

for frame in range(50,num_frames,1) :
    t0 = time.time()  # start time
    data = pipeline.compute(frame)
    cell=data.cell[:3,:3]
    L = [cell[0,0],cell[1,1]]
    N = data.particles.count
    
    bonds = data.particles.bonds
    topology = data.particles.bonds.topology
    positions = data.particles.positions
    
    
    # bond_vectors += np.dot(data.cell[:3,:3], data.particles.bonds.pbc_vectors.T).T
    # norms = np.linalg.norm(bond_vectors,axis=1)
    # # print(norms)
    # normed = bond_vectors/norms[:,None]
    # plt.hist(normed[:,0])
    # c = np.append(c,normed)
    # print(len(bonds), len(data.particles))

    bonds_enum = BondsEnumerator(data.particles.bonds)

    for particle_index in range(data.particles.count):
        # Loop over bonds of current atom.
        bond_vectors = [] 

        for bond_index in bonds_enum.bonds_of_particle(particle_index):
            # Obtain the indices of the two particles connected by the bond:
            a = topology[bond_index, 0]
            b = topology[bond_index, 1]

            if a==particle_index:
                j = b
            else:
                j = a

            d = positions[j] - positions[particle_index]
            for i in range(2):
                if d[i]>rc:
                    d[i]-=L[i]
                elif d[i]<-rc:
                    d[i]+=L[i]
        

            bond_vectors.append(d[:2])

        bond_vectors = np.array(bond_vectors)
        if bond_vectors.shape[0]>1:
            norms = np.linalg.norm(bond_vectors,axis=1)
        #     # print(norms)
            normed = bond_vectors/norms[:,None]

            D = np.dot(normed,normed.T)
            cosines  =  D[np.triu_indices_from(D, k=1)]
            if min(np.arccos(cosines))<0.1:
                print("wARNING: too small angle !")
            c.extend(cosines)

    t1 = time.time()  # start time

    print(frame, N, bonds.count,t1-t0)
    
# np.savetxt(path+".cosines.txt", c)

H , e = np.histogram(np.arccos(c),bins = np.linspace(0,np.pi,128),density=True)

x = e[:-1]+0.58*(e[1]-e[0])
# plt.plot(x,H)
np.savetxt(path+f".bondangle-rc{rc}.txt",list(zip(x,H)))


# plt.show()

    # print(data.particles.bonds)

