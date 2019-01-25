from __future__ import print_function
import numpy as np

def get_length(a,warnings=False):
  dif = np.ediff1d(a, to_begin=a[0]-a[-1])
  dif = np.ediff1d(a, to_end=a[0]-a[-1])
  idx = np.arange(len(a))

  end=idx[dif==-1]
  begin=idx[dif==1]
  if len(end)==0:
    print ("    ! Warning: no gap found")
    return [len(a)]

  if end[0]<begin[0]:
    end = np.roll(end, -1)
    chords = end-begin
    chords[-1] = len(a)-begin[-1]+end[-1]
    return chords
  else:
    return end-begin
