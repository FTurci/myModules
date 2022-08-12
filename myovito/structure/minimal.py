from ovito.io import import_file
from ovito.modifiers import *
import time
import sys

path = sys.argv[1]


pipeline = import_file(path, multiple_frames=True)
num_frames = pipeline.source.num_frames

for frame in range(1,num_frames,1) :
    t0 = time.time()  # start time
    data = pipeline.compute(frame)
    t1 = time.time()  

    print(frame,t1-t0)

