import numpy as np

def comupute_isf_d(trajectory,d,lengthscale=1.,threshold=1):
    """Compute the self part of the Intermediate Scattering Function in d-dimensions."""
    from scipy.special import jv
    assert  trajectory.ensemble =='NVT', "Error: ISF can be computed on NVT data only"
    nsamples=np.zeros(trajectory.nframes)
    isf=np.zeros(trajectory.nframes)
    from tqdm import tqdm
    for t in  tqdm(xrange(0,trajectory.nframes-1)):
        for tt in xrange(t+1, trajectory.nframes):
            if nsamples[tt-t]<threshold:
                # compute centres of mass
                r1=trajectory.Frames[t].r[:,:d]
                r2=trajectory.Frames[tt].r[:,:d]
                
                r1cm=r1.mean(axis=0)
                r2cm=r2.mean(axis=0) 
                dr=np.absolute((r2-r2cm)-(r1-r1cm))
            
                for k in xrange(d):
                    dr[dr>trajectory.lengths[d]*0.5]-=trajectory.lengths[d]
                dr=np.sqrt((dr**2).sum(axis=1))
                
                nsamples[tt-t]+=1.
                isf[tt-t]+=(jv(0, 2*np.pi*dr/lengthscale)).sum()

    return (isf/(nsamples*trajectory.N))[1:], nsamples[1:]