from __future__ import print_function
import numpy as np
import pickle

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def decompress(filename,mode='rb',delete=False):
  import gzip 
  if filename[-3:]==".gz":
    with gzip.open(filename, mode) as fin:
      data = fin.read()
      filenameout = filename+".decompressed"
      with open(filenameout ,'wb') as fout:
        fout.write(data)
      return filenameout
  else:
    print ("!!! Decompression failed: only .gz files are supported.")

def clean(filename,extension = ".decompressed"):
  import os
  os.system("rm "+filename+extension)

def cast_to_array(a):
  b = tuple([np.array(i) for i in a])
  return b