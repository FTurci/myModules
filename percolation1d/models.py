import numpy as np
from . import helpers
from numpy.random import uniform
from scipy import ndimage


class PercolationStudy(object):
  """Study percolation. 1D case."""
  def __init__(self,probability,numsites):
    super(PercolationStudy, self).__init__()
    self.probability = probability
    self.L = numsites
    self.create_configuration()

  def create_configuration(self):
    self.configuration = helpers.generate(self.L, self.probability)

  def compute_clusters(self):
    self.clusters =  helpers.clusters(self.configuration)
    self.n_s = np.bincount(self.clusters.astype(int))*1./self.L
    self.s = np.arange(len(self.n_s))

  def compute_mean(self):
    self.mean = self.clusters.mean()

  def compute_average(self,iterations=10000):
    self.average = helpers.average(self.configuration,iterations)

  def theory_mean(self):
    return helpers.theory_mean(self.probability)

  def theory_average(self):
    return helpers.theory_average(self.probability)

  def theory_n_s(self):
    return (1-self.probability)**2*self.probability**(self.s)



def test_plots_n_s(probabilities,sites):
  import pylab as pl

  for p in probabilities:
    print (p)
    study = PercolationStudy(p, sites)
    study.compute_clusters()
    pl.loglog(study.s[1:], study.n_s[1:], label=str(p))
    pl.loglog(study.s[1:], study.theory_n_s()[1:],'k--', dash_capstyle="round", lw=0.5)
  pl.xlabel("s")
  pl.ylabel("n(s,p)")
  pl.legend()
  pl.show()

