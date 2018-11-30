from __future__ import print_function
import numpy as np

def get_length(a,warning=False):
  dif = np.ediff1d(a, to_begin=a[0]-a[-1])
  dif = np.ediff1d(a, to_end=a[0]-a[-1])
  idx = np.arange(len(a))

  end=idx[dif==-1]
  begin=idx[dif==1]
  if len(end)==0:
    if warning :
      print ("    ! Warning: no gap found")
    return [len(a)]

  if end[0]<begin[0]:
    end = np.roll(end, -1)
    chords = end-begin
    chords[-1] = len(a)-begin[-1]+end[-1]
    return chords
  else:
    return end-begin

class ChordLengthAnalyser(object):
  def __init__(self,array):
    self.data = array
    self.ndim = len(self.data.shape)
    self.shape = self.data.shape

  def compute(self,warning=False):
    lengths = []
    if self.ndim==1:
      lengths.append( get_length(self.data))
    if self.ndim==2:
      for i in range(self.shape[0]):
        lengths.append(get_length(self.data[i,:]))
      for j in range(self.shape[1]):
        lengths.append(get_length(self.data[:,j]))

    if self.ndim == 3:
      #along x
      for j in range(self.shape[1]):
        for k in range(self.shape[2]):
          lengths.append(get_length(self.data[:,j,k]))
      #along y
      for i in range(self.shape[0]):
        for k in range(self.shape[2]):
          lengths.append(get_length(self.data[i,:,k]))
      #along z
      for i in range(self.shape[0]):
        for j in range(self.shape[1]):
          lengths.append(get_length(self.data[i,j,:]))

    self.lengths = np.concatenate(lengths)


