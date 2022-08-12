import numpy as np
from natsort import realsorted
import glob
import matplotlib.pyplot as plt
from stringato import extract_floats as ef

import sys
pattern = sys.argv[1]
files = realsorted(glob.glob(pattern))
print(files)

N=8640
x = np.arange(9)
for f in files:
    data = np.loadtxt(f)
    framebyframe = data.reshape(N,int(data.shape[0]/N))
    count=0

    hists = []
    for frame in range(framebyframe.shape[1]):
        H,edges = np.histogram(framebyframe[:,frame],bins=x-0.5,density=True)
        hists.append(H)
        count+=1
        print(H)

    print(f)
    print(H)
    plt.errorbar(x[:-1],np.mean(hists,axis=0),np.std(hists,axis=0), label= ef(f)[0])
   
plt.legend()
plt.xlabel("coordination $r_{cut}<1$")
plt.savefig("coordination.pdf")
plt.show()