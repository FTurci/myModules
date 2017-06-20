import numpy as np
import pandas as pd
def unfold_coordinates(trajectory,d):
    unfolded_r=np.zeros((trajectory.N, d,trajectory.nframes))
    unfolded_r[:,:,0]=trajectory.Frames[0].r[:,:d]
    for i in xrange(1,trajectory.nframes):
            dr=trajectory.Frames[i].r[:,:d]-trajectory.Frames[i-1].r[:,:d]
            for k in xrange(d):
                    L=trajectory.Frames[i].lengths[k]
                    dr[:,k]=dr[:,k]-np.round(dr[:,k]/L)*L
                    if np.any(np.abs(dr[:,k])>=L*0.5):
                        print "too large displacement"
            unfolded_r[:,:,i]=unfolded_r[:,:,i-1]+dr
    # print unfolded_r.shape
    return unfolded_r

def compute_isf_d(trajectory,d,lengthscale=1.,threshold=10, unfold=True):
    """Compute the self part of the Intermediate Scattering Function in d-dimensions."""
    if d==2:
        from scipy.special import jv
    assert  trajectory.ensemble =='NVT', "Error: ISF can be computed on NVT data only"
    nsamples=np.zeros(trajectory.nframes)
    isf=np.zeros(trajectory.nframes)
    # unfold the trajectory
    if unfold:
        r=unfold_coordinates(trajectory,d)
    else:
        r=np.zeros((trajectory.N, d,trajectory.nframes))
        for i in xrange(1,trajectory.nframes):
            r[:,:,i]=trajectory.Frames[i].r[:,:d]


    from tqdm import tqdm
    for t in  tqdm(xrange(0,trajectory.nframes-1)):
        # print t
        for tt in xrange(t+1, trajectory.nframes):
            # print tt
            if nsamples[tt-t]<threshold:
                # compute centres of mass
                r1=r[:,:,t]
                r2=r[:,:,tt]
                
                r1cm=r1.mean(axis=0)
                r2cm=r2.mean(axis=0) 
                dr=np.absolute((r2-r2cm)-(r1-r1cm))
                
                for k in xrange(d):
                    # print "k",k,trajectory.lengths[d]
                    dr[dr>trajectory.lengths[k]*0.5]-=trajectory.lengths[k]
                dr=np.sqrt((dr**2).sum(axis=1))
                
                nsamples[tt-t]+=1.
                if d==2:
                    isf[tt-t]+=(jv(0, 2*np.pi*dr/lengthscale)).sum()
                if d==3:
                    isf[tt-t]+=(np.sinc(2*dr/lengthscale)).sum() #sinc is defined as sin(pi x)/(pi x)

    return (isf/(nsamples*trajectory.N))[1:], nsamples[1:]

def read_log(filename):
    print "Reading logfile", filename
    thermo=[]
    # read the thermo output
    # go to line befor the opuput
    read_thermo=False
    with open(filename,'r') as fin:
        for line in fin:     
            split = line.split()
            if len(split)>0:
                if read_thermo:
                    thermo.append(line)

                if split[0]=='Memory':
                    read_thermo=True
                if split[0]=="Loop":
                    read_thermo=False

    # construct a pandas data structure
    thermo=thermo[:-1]
    if len(thermo)==0:
        print "===> ERROR: Logfile",filename, "empty. Aborting."
        import sys
        sys.exit(0)
    header=thermo[0].split()
    thermo=[line.split() for line in thermo]
    table=np.array(thermo[1:]).astype(float)
    ncols=table.shape[1]
    dictionary={}
    for i in range(ncols):
        dictionary[header[i]]=table[:,i]
    dataFrame=pd.DataFrame(dictionary)
    return dataFrame

