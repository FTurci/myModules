import numpy as np
import sys
import fileformats 

reload(fileformats)



def rand_rotation_matrix3d(deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
    
    if randnums is None:
        randnums = np.random.uniform(size=(3,))
        
    theta, phi, z = randnums
    
    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0*deflection  # For magnitude of pole deflection.
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )
    
    st = np.sin(theta)
    ct = np.cos(theta)
    
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M


def digitize(coords, minedge=[0,0,0],maxedge=[1.,1.,1.], divisions=[3.,3.,3.]):
    """
    Assigns ids for cuboid regions in space. Returns ids of the cuboids, their centers and their indices
    """
    bins=[]
    digitized=[]
    dim=coords.shape[1]

    for d in xrange(dim):
        bins.append(np.linspace(minedge[d], maxedge[d], divisions[d]+1))
        digitized.append(np.digitize(coords[:,d], bins[d])-1)
    # 2d
    if dim==2:
        centers=np.meshgrid(bins[0][:-1]+(bins[0][1:]-bins[0][:-1])*0.5, bins[1][:-1]+(bins[1][1:]-bins[1][:-1])*0.5)
        index= digitized[0]+divisions[1]*digitized[1]
    else:#3d
        centers=np.meshgrid(bins[0][:-1]+(bins[0][1:]-bins[0][:-1])*0.5, bins[1][:-1]+(bins[1][1:]-bins[1][:-1])*0.5, bins[2][:-1]+(bins[2][1:]-bins[2][:-1])*0.5)
        index= digitized[0]+divisions[1]*(digitized[1]+divisions[2]*digitized[2])
    return index, centers, digitized

def spherical_subsamples(coords, minedge=[0,0,0],maxedge=[1.,1.,1.], divisions=[3.,3.,3.], R=1):
    """
    Isolates spherical subsamples of radius R located on a mesh from a set of particle coordinates. It returns an array of labels for each particle.
    """
    
    dim=coords.shape[1]
    ids,c,digitized=digitize(coords, minedge,maxedge,divisions)
    from scipy.spatial.distance import cdist
    
    if dim==2:
        centers=np.array([c[0].flatten(), c[1].flatten()]).T
    if dim==3:
        centers=np.array([c[0].flatten(), c[1].flatten(),c[2].flatten()]).T

    distances=cdist(coords, centers)
    print "Radius",R
    ids[np.logical_not(np.any(distances<R,axis=1))]=-1
    
    samples=[]
    for sample in xrange(int(np.prod(divisions))):
        c=coords[ids==sample]
        baricenter=c.mean(axis=0)
        samples.append(c-baricenter)
    return samples

def subsamples_n(coords, minedge=[0,0,0],maxedge=[1.,1.,1.],divisions=[3.,3.,3.], N=64):
    """
    Isolates spherical subsamples of N particles located on a mesh from a set of particle coordinates. It returns an array of labels for each particle.
    """
    s=[]
    s_ids=[]
    ids=np.arange(coords.shape[0])
    L=np.array(maxedge)-np.array(minedge)
    dx=np.linspace(minedge[0]+L[0]/divisions[0]/2., maxedge[0]-L[0]/divisions[0]/2.,divisions[0])
    dy=np.linspace(minedge[1]+L[1]/divisions[1]/2., maxedge[1]-L[1]/divisions[1]/2.,divisions[1])
    dz=np.linspace(minedge[2]+L[2]/divisions[2]/2., maxedge[2]-L[2]/divisions[2]/2.,divisions[2])

    X,Y,Z=np.meshgrid(dx,dy,dz)
    centers = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
    # print centers

    for c in xrange(centers.shape[0]):
        # print centers[c]
        dist2=np.sum((centers[c]-coords)**2, axis=1)
        # print dist2.shape
        sel=dist2.argsort()
        s.append(coords[sel][:N]-centers[c])
        selected_ids=ids[sel][:N]
        s_ids.append(selected_ids)
        # print "lens",len(s)
    return s, s_ids





def findQ(c1,c2, a):
    from scipy.spatial.distance import cdist
    c1=np.array(c1)
    c2=np.array(c2)
    dists=cdist(c1, c2)
    # print dists
    # print dists.shape, dists.min()
    q=np.sum(np.any(dists<=a, axis=1))/float(c1.shape[0])
    return q

def overlap(c1,c2, nrotations=100, a=0.3):
    """
    Computes and returns the maximum overlap between two spherical configurations c1,c2. It performs nrotations and returns the highest overlap value.
    """
    qs=[]
    from scipy.spatial.distance import cdist
    for x in xrange(nrotations):
        M=rand_rotation_matrix3d()
        if(x>0):
            rotated_c2=(np.dot(M,c2.T)).T
        else:
            rotated_c2=c2
        dists=cdist(c1, rotated_c2)
        # print dists.shape, dists.min()
        q=np.sum(np.any(dists<=a, axis=1))/float(c1.shape[0])
        qs.append(q)
    return np.max(qs),qs

from scipy.spatial.distance import cdist


def maximise_overlap(c1,c2, a):
    A,B=np.mat(c1),np.mat(c2)
    # find a rotation matrix based on a subset of N and find the best one
    qmax=0
    f=1
    for i in range(f,A.shape[0]/f):
        R,t=rigid_transform_3D(A[:i*f],B[:i*f])
        q=findQ(B, (R*A.T+t).T, a)
        if q>qmax:
            qmax=q
    # check with all the matrix
    R,t=rigid_transform_3D(A,B)
    q=findQ(B, (R*A.T+t).T, a)
    if q>qmax:
            qmax=q
    print "    q", qmax
    if qmax==0 :
        import sys
        fileformats.array_to_xyz(c1, "A.xyz")
        fileformats.array_to_xyz(c2, "B.xyz")
        sys.exit("!!!Found zero overlap!!!")
    return qmax

def ensemble_overlap_constant_N(samples,nrotations=100,diameter=1, overlap_scale=0.3):
    qs=[]
    for i in xrange(len(samples)-1):
        for j in xrange(i+1,len(samples)):
            print "===>",i,j
            # identify particles in sphere i and sphere j
            c1,c2=samples[i],samples[j]
            q=maximise_overlap(c1, c2,overlap_scale*diameter )
            qs.append(q)
    return np.array(qs)


def ensemble_overlap(samples,nrotations=100,diameter=1, overlap_scale=0.3):
    from scipy.spatial.distance import pdist
    qs=[]
    for i in xrange(len(samples)-1):
        for j in xrange(i+1,len(samples)):
            
            # identify particles in sphere i and sphere j
            c1,c2=samples[i],samples[j]
            q,qqs=overlap(c1, c2,nrotations=nrotations,a=overlap_scale*diameter )
            print "==>",i,j,q
            qs.append(q)
    return np.array(qs)


def random_overlaps2d(rho,a=0.3,L=100., nsamples=100):
    N=int(L*L*rho)#unit surface 
    from scipy.spatial.distance import cdist
    c0=np.random.uniform(0,L, size=(N,2))
    qs=[]
    for i in xrange(nsamples):
        c=np.random.uniform(0,L, size=(N,2))
        dists=cdist(c0, c)
        q=np.sum(np.any(dists<a, axis=1))/float(N)
        qs.append(q)

    expected=rho*np.pi*a**2
    return qs, np.mean(qs), expected

def get_lengths(xyzfile,L,rotations=100 , diameter=18., divide=np.arange(3,11,2)):
    
    coords=np.loadtxt(xyzfile, skiprows=2, usecols=[1,2,3])

    qmeans=[]
    # 0.60 -> 344
    # 0.58 -> 320
    # 0.52 -301
    lengths=float(L)/divide
    for division in divide:
        samples=spherical_subsamples(coords, maxedge=[L,L,L],divisions=[division,division,division], R=L/division/2. )
        Q=ensemble_overlap(samples,nrotations=rotations,diameter=diameter)
        qmeans.append(Q.mean())
    

    np.savetxt(xyzfile+".qmeans", zip(lengths/float(diameter),qmeans))

def get_lengths_random(phi,rotations=100 , diameter=18., divide=np.arange(3,11,2)):
    N=5000
    L=(N/(phi/(4./3*np.pi*(diameter/2.)**3)))**(1./3.)
    coords=np.random.uniform(0,L,size=(N,3))
    qmeans=[]
    # 0.60 -> 344
    # 0.58 -> 320
    # 0.52 -301
    lengths=float(L)/divide
    for division in divide:
        samples=spherical_subsamples(coords, maxedge=[L,L,L],divisions=[division,division,division], R=L/division/2. )
        Q=ensemble_overlap(samples,nrotations=rotations,diameter=diameter)
        qmeans.append(Q.mean())
    

    np.savetxt("random_phi%g.qmeans"%phi, zip(lengths/float(diameter),qmeans))
    return zip(lengths/float(diameter),qmeans)




def rigid_transform_3D(A, B):
    # using singular value decomposition to optimise the rotation
    # see http://nghiaho.com/?page_id=671
    A=np.mat(A)
    B=np.mat(B)
    assert len(A) == len(B)

    N = A.shape[0]; # total points
    # print N,"points"
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    
    # centre the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = np.transpose(AA) * BB
    U, S, Vt = np.linalg.svd(H)

    R = Vt.T * U.T

    # special reflection case
    if np.linalg.det(R) < 0:
       # print "Reflection detected"
       Vt[2,:] *= -1
       R = Vt.T * U.T

    t = -R*centroid_A.T + centroid_B.T

    # print t

    return R, t


def gyration_tensor(r, Print=False):
    from scipy.linalg import eig
    import sys 
    N=r.shape[0]
    # centre of mass
    r_cm=np.mean(r,axis=0)
    # remove the centre of mass
    r-=r_cm
    # initialize the tensor to zero
    S=np.zeros((3,3))
    S=np.matrix(S)
    # compute the entries
    for m in range(3):
        for n in range(3):
            S[m,n]=np.sum(r[:,m]*r[:,n])/N
    # diagonalise and derive the eigenvalues L and the eigenvectors V
    LL,V=eig(S)
    # sort the eigenvalues
    L=np.real(np.sort(LL))
    Rg=np.sqrt(np.sum(LL))
    b=L[2]-0.5*(L[1]+L[2])
    c=L[1]-L[0]
    k2=(b**2+(3./4.)*c**2)/Rg**4
    if Print:
        print "eigenvalues ",L
        print "Radius of gyration ",Rg
        print "asphericity ",b
        print "acilindricity ",c
        print "relative shape anisotropy ",k2
    return r_cm,V,LL,Rg,b,c,k2

def matrix_to_quaternions(M):
    m=np.mat(M)
    trace=np.trace(m)
    if trace>0:
        S=np.sqrt(trace+1.)*2
        qw = 0.25 * S
        qx = (m[2,1] - m[1.2]) / S
        qy = (m[0,2] - m[2.0]) / S
        qz = (m[1,0] - m[0.1]) / S
    elif m[0,0] > m[1,1] and m[0,0] > m[2,2]:
        S = np.sqrt(1.0 + m[0,0] - m[1,1] - m[2,2]) * 2
        qw = (m[2,1] - m[1,2]) / S
        qx = 0.25 * S
        qy = (m[0,1] + m[1,0]) / S 
        qz = (m[0,2] + m[2,0]) / S 
    elif (m[1,1] > m[2,2]) : 
        S = np.sqrt(1.0 + m[1,1] - m[0,0] - m[2,2]) * 2
        qw = (m[0,2] - m[2,0]) / S
        qx = (m[0,1] + m[1,0]) / S 
        qy = 0.25 * S
        qz = (m[1,2] + m[2,1]) / S   
    else :
        S = np.sqrt(1.0 + m[2,2] - m[0,0] - m[1,1]) * 2
        qw = (m[1,0] - m[0,1]) / S
        qx = (m[0,2] + m[2,0]) / S
        qy = (m[1,2] + m[2,1]) / S
        qz = 0.25 * S;
    return qx,qy,qz,qw



def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

def voronoi2d(Frame,maximumdistance,radii, save=True):
    shelveName="frame%d.shelve"%Frame.timeframe
    import os
    import shelve

    if os.path.isfile(shelveName):
        print "File", shelveName, "already exists"
        d=shelve.open(shelveName)
        return d['cells']
    else:
        coords=Frame.r[:,:2]
        import pyvoro
        cells=pyvoro.compute_2d_voronoi(coords,Frame.boxinfo[:2],maximumdistance, periodic=[True,True], radii=radii)
        if save:
            
            d=shelve.open(shelveName)
            d['cells']=cells
            d.close()

        return cells
