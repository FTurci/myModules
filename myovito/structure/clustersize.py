from ovito.io import import_file
import numpy as np
from ovito.modifiers import *
from stringato import extract_floats as ef
import glob


files = glob.glob("tj*gz")
N = 2000
for f in files:
    pipeline = import_file(f, multiple_frames=True)
    # print(f)
    A = ef(f)[1]

    pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 2.4))
    pipeline.modifiers.append(ExpressionSelectionModifier(
        expression = 'Coordination < 14'))

    pipeline.modifiers.append(DeleteSelectedModifier())
    
    pipeline.modifiers.append(ClusterAnalysisModifier(
        cutoff=1.2,
        sort_by_size=True,
        ))

    cl = []
    for frame in [-20,-10,-1]:
        data = pipeline.compute(pipeline.source.num_frames+frame)
        cluster = data.particles.cluster.array
        # print(cluster)
        cl.append(len(cluster[cluster==1])/N)

    print(A, np.mean(cl), np.std(cl))

