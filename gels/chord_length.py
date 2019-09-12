from __future__ import print_function
import numpy as np

def get_length(a,warning=False, pbc=True):
  dif = np.ediff1d(a, to_begin=a[0]-a[-1])
  idx = np.arange(len(a))
  end=idx[dif==-1]
  begin=idx[dif==1]
  if  pbc ==True:
    if len(end)==0:
      if warning:
        print ("    ! Warning: no gap found")
      if a[-1]==1:#if it is on
        return [len(a)]
      else:
        return [0]

    if end[0]<begin[0]:
      end = np.roll(end, -1)
      chords = end-begin
      chords[-1] = len(a)-begin[-1]+end[-1]
      return chords
    else:
      return end-begin
  else:
    return end-begin

class ChordLengthAnalyser(object):
  def __init__(self,array,pbc=True):
    self.data = array
    self.ndim = len(self.data.shape)
    self.shape = self.data.shape
    self.pbc =pbc
  def compute(self,warning=False,remove_zeros=True):
    lengths = []
    if self.ndim==1:
      lengths.append( get_length(self.data,pbc=self.pbc,warning=warning))
      self.lengths = np.concatenate(lengths)
    if self.ndim==2:
      lengthx=[]
      lengthy=[]
      for i in range(self.shape[0]):
        l=get_length(self.data[i,:],pbc=self.pbc,warning=warning)
        lengthx.append(l)
      for j in range(self.shape[1]):
        lengthy.append(get_length(self.data[:,j],pbc=self.pbc,warning=warning))
      self.lengthx=np.concatenate(lengthx)
      self.lengthy=np.concatenate(lengthy)

      if remove_zeros:
        self.lengths=[l[l>0] for l in [self.lengthx,self.lengthy]]
      else:
        self.lengths=[l for l in [self.lengthx,self.lengthy]]

    if self.ndim == 3:
      #along x
      lengthx=[]
      for j in range(self.shape[1]):
        for k in range(self.shape[2]):
          lengthx.append(get_length(self.data[:,j,k],pbc=self.pbc,warning=warning))
      #along y
      lengthy=[]
      for i in range(self.shape[0]):
        for k in range(self.shape[2]):
          lengthy.append(get_length(self.data[i,:,k],pbc=self.pbc,warning=warning))
      #along z
      lengthz=[]
      for i in range(self.shape[0]):
        for j in range(self.shape[1]):
          lengthz.append(get_length(self.data[i,j,:],pbc=self.pbc,warning=warning))

      self.lengthx=np.concatenate(lengthx)
      self.lengthy=np.concatenate(lengthy)
      self.lengthz=np.concatenate(lengthz)
      if remove_zeros:
        self.lengths=[l[l>0] for l in [self.lengthx,self.lengthy,self.lengthz]]
      else:
        self.lengths=[l for l in [self.lengthx,self.lengthy,self.lengthz]]

