import numpy as np
import pyvoro
import matplotlib.pyplot as pl
class tessellation:
    """Builds a Voronoi tessellation from a Frame object (from fileformats.py)"""
    def __init__(self,Frame,dimension, radii=None, periodic=[True,True]):
        print "Computing Voronoi tessellation in %d-d"%dimension
        if dimension==2:
            
            if radii==None:
                radii=np.ones(Frame.N)
            self.cells=self.voronoi_tessellation_2d(Frame, radii,periodic)

    def voronoi_tessellation_2d(self,Frame,radii, periodic):
        coords=Frame.r[:,:2]
        box=Frame.boxinfo
        print coords.shape
        cells=pyvoro.compute_2d_voronoi(coords, [box[0],box[1]],coords.shape[0],  periodic=periodic)
        
        # area of the disks
        self.areas=radii**2*np.pi
        # compute area fraction per cell
        self.volumes=np.zeros(Frame.N)

        for i in xrange(Frame.N):
            # local packing fraction per cell
            self.volumes[i]=cells[i]['volume']

        return cells 

class density_map:
    def __init__(self, Frame, dimension,divisions=[10,10,10]):
        if dimension==2:
            box=Frame.boxinfo
            bins=np.array(
                # [np.linspace(box[0].min(),box[0].max(),divisions[0]),
                # np.linspace(box[1].min(),box[1].max(),divisions[1])]
                [np.linspace(0,1,divisions[0]),
                np.linspace(0,1,divisions[1])]
            )
            r=Frame.r
            self.H=np.histogram2d(r[:,0], r[:,1], bins=bins)


