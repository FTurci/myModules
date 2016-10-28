from matplotlib import pyplot as pl
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
import numpy as np

def draw_disks(axes,coords, r,colors=None,linewidth=0.5, edgecolor='none', cmap=pl.cm.viridis):
    patches=[]
    for x,y,radius in zip(coords[:,0], coords[:,1], r):
        circle = Circle((x,y), radius, linewidth=linewidth,edgecolor=edgecolor)
        patches.append(circle)
    p = PatchCollection(patches, cmap=cmap)
    if colors is not None:
        p.set_array(colors)
    axes.add_collection(p)
    axes.axis('equal')
    axes.set_xlim(round(coords[:,0].min()), round(coords[:,0].max()))
    axes.set_ylim(round(coords[:,1].min()), round(coords[:,1].max()))
