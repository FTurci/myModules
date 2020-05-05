from matplotlib import pyplot as pl
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
import numpy as np
from scipy.interpolate import UnivariateSpline
from types import SimpleNamespace  
# Beautiful colors
colors = {}
colors['red'] = '#aa2217'
colors['blue'] = '#0096ff'
colors['moss'] = '#008f51'
colors['ocean'] = '#008f51'
colors['orange'] = '#ff9300'
colors['yellow'] = '#efbc40'
colors['lavender'] ='#d783ff'
colors['coffee'] = '#d783ff'

colors = SimpleNamespace(**colors)


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
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


def confidence(x,y,err,color='#0096ff',trend=True,s=None,label=None,**kwargs):
    ax = kwargs.pop("axes", None)
    if not ax:
        fig, ax = pl.subplots()
    
    spline =UnivariateSpline(x,y,1/err,s=s)
    erspline = UnivariateSpline(x,err,s=s)

    ax.errorbar(x,y,err,color=color,elinewidth=0.8,label=label,**kwargs)
    _x = np.linspace(x.min(), x.max(), len(x)*10)
    _y = spline(_x)
    _err = erspline(_x)
    #ax.fill_between(x,y-err,y+err, color=color, alpha=0.5)
    ax.fill_between(_x,_y-_err,_y+_err, color=color, alpha=0.25)
    if trend:
        ax.plot(_x,_y,color=color) 
