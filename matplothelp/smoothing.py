from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as pl
class SmoothLine(object):
  """Smoothing x,y data """
  def __init__(self, x,y,smooth=None,knots =3,npoints=100):
    super(SmoothLine, self).__init__()
    x = np.array(x)
    y = np.array(y)
    order = x.argsort()
    self.x = x[order]
    self.y = y[order]
    self.xnew = np.linspace(min(x),max(x),npoints)
    if smooth==None:
      smooth = len(x)
    self.spline = UnivariateSpline(self.x,self.y, k=knots, s=smooth)
    self.ynew = self.spline(self.xnew)

  def plot(self,axis=None,**kwargs):
    if axis is None:
      pl.plot(self.xnew, self.ynew, **kwargs)
    else:
      axis.plot(self.xnew, self.ynew, **kwargs)

class ModelLine(object):
  """Smoothing x,y data """
  def __init__(self, x,y,model='power+const',npoints=100):
    super(ModelLine, self).__init__()
    x = np.array(x)
    y = np.array(y)
    order = x.argsort()
    self.x = x[order]
    self.y = y[order]
    self.xnew = np.linspace(min(x),max(x),npoints)

    if model=='power+const':
      self.ynew= self.powerconst()

  def powerconst(self):
    def model(x,a,b,c):
      return a+b*x**c
    
    self.fitmodel(model)
    return model(self.xnew, *self.popt)

  def fitmodel(self,model):
    self.popt,self.pcov=curve_fit(model,self.x,self.y, ) 

  def plot(self,axis=None,**kwargs):
    if axis is None:
      pl.plot(self.xnew, self.ynew, **kwargs)
    else:
      axis.plot(self.xnew, self.ynew, **kwargs)