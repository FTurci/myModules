import numpy as np
from natsort import realsorted
import glob
import matplotlib.pyplot as plt
import sys
import re


pattern = sys.argv[1]
string = sys.argv[2]

files = realsorted(glob.glob(pattern))
for f in files:
    x,y = np.loadtxt(f, unpack=True)
    print(len(string))
    v = re.findall(r'{}\d+(?:\.\d+)?'.format(string),f)[0]
    label = v
    plt.plot(x/np.pi*3,y, label=string+"="+label[len(string):])
plt.xlabel(r"$3\theta/\pi$")
plt.legend()
plt.show()