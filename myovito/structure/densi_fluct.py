from ovito.io import import_file
import numpy as np
from ovito.modifiers import *
from stringato import extract_floats as ef
import glob
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import sys



path = sys.argv[1]
nbins = int(sys.argv[2])

pipeline = import_file(path, multiple_frames=True)
num_frames = pipeline.source.num_frames

A = ef(path)[0]

B = ef(path)[1]
Pe = ef(path)[2]
rho = ef(path)[3]
vs,ms = [],[]
for frame in range(50,num_frames,2) :
    data = pipeline.compute(frame)
    N = data.particles.count
    # print(frame, N)
    Ls = (np.diag(data.cell[:3,:3]))
    origin = data.cell_[:,3]
    # print(origin)
    pos = data.particles.positions.array[:,:2]
    binx = np.linspace(origin[0],origin[0]+Ls[0],nbins)
    biny = np.linspace(origin[0],origin[0]+Ls[0],nbins)
    H,e = np.histogramdd(pos,bins=(binx,biny))
    # print(H.shape,H)

    variance = np.var(H)
    mean = np.mean(H)

    vs.append(variance)
    ms.append(mean)
ms = np.array(ms)
vs = np.array(vs)

np.savetxt('data.txt', list(zip(ms,vs)))
print(A,B, Pe,rho, np.mean(vs)/np.mean(ms), np.mean(vs/ms),nbins) 
   # plt.imshow(H)
    # plt.savefig(f"h{frame}.png")

    # plt.hist(coord, bins = np.arange(10)-0.5,alpha=0.1)
    # cluster = data.particles.cluster.array
    # print(cluster)
    # cl.append(len(cluster[cluster==1])/N)

# np.savetxt(path+".coordination.txt",c)
# print(A, np.mean(cl), np.std(cl))

