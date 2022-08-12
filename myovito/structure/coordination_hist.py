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


pipeline = import_file(path, multiple_frames=True)
num_frames = pipeline.source.num_frames


pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 1.0))
# pipeline.modifiers.append(ExpressionSelectionModifier(
#     expression = 'Coordination < 14'))

# pipeline.modifiers.append(DeleteSelectedModifier())

# pipeline.modifiers.append(ClusterAnalysisModifier(
#     cutoff=1.2,
#     sort_by_size=True,
#     ))

cl = []

c = []
for frame in range(50,num_frames,2) :
    data = pipeline.compute(frame)
    N = data.particles.count
    print(frame, N)
    coord = data.particles.coordination.array
    c= np.append(c,coord)

    # plt.hist(coord, bins = np.arange(10)-0.5,alpha=0.1)
    # cluster = data.particles.cluster.array
    # print(cluster)
    # cl.append(len(cluster[cluster==1])/N)

np.savetxt(path+".coordination.txt",c)
# print(A, np.mean(cl), np.std(cl))

