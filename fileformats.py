import pickle
import numpy as np

class frame:
    def __init__(self,table, timeframe=None, boxinfo=None):
        self.N=table.shape[0]
        self.ids=table[:,0].astype(int)
        self.types=table[:,1].astype(int)
        self.r=table[:,2:]
        
        
        if timeframe is not None:
            self.timeframe=timeframe
        if boxinfo is not None:
            self.boxinfo=np.array(boxinfo)
        else:
            print "here"
            self.boxinfo=np.array([
        np.array(self.r[:,0].min(),self.r[:,0].max()),
        np.array(self.r[:,1].min(),self.r[:,1].max()),
        np.array(self.r[:,2].min(),self.r[:,2].max())
        ])

        self.lengths=self.boxinfo[:,1]-self.boxinfo[:,0] 

class trajectory:
    def __init__(self,filename,compression=None, format='atom', numframes='all',ensemble='NVT', every=1):
        if format=='atom':
            if compression is None:
                filehandle=open(filename, 'r')
            elif compression in 'gz':
                import gzip
                filehandle=gzip.open(filename, 'r')
            self.Frames=[]
            count=0
            buf=True
            while True:
                if not buf: 
                    break
                elif count%every==0:
                    Frame,buf=read_lammps_atom_frame(filehandle)
                    if not buf:
                        break
                    self.Frames.append(Frame)
                    count+=1
                    if numframes != 'all':
                        if count>=numframes:
                            break
                else:
                    # read fast without storing anything
                    for line in xrange(3):
                        buf=filehandle.readline()
                    if not buf:
                        break
                    N=int(filehandle.readline())
                    for line in xrange(5+N):
                        buf=filehandle.readline()
                    count+=1
                    if numframes != 'all':
                        if count>=numframes:
                            break

        
        self.nframes=len(self.Frames)

        if ensemble=='NVT':
            self.ensemble='NVT'
            self.boxinfo=np.array(self.Frames[0].boxinfo)
            self.N=self.Frames[0].N
            self.lengths=self.Frames[0].lengths
        if ensemble =='NPT':
            self.ensemble='NPT'
            self.N=self.Frames[0].N
            boxinfo=np.zeros((3,2))
            for frame in xrange(self.nframes): 
                boxinfo+=self.Frames[frame].boxinfo
            self.boxinfo=boxinfo/self.nframes



def read_lammps_atom_frame(filehandle):
    """Read frame from LAMMPS trajectory in the atom format."""
    buf=filehandle.readline() #ITEM: TIMESTEP
    if(len(buf)>0): 
        print "Reading",buf
    else:
        return False,False
    timeframe=int(filehandle.readline())
    print timeframe
    buf=filehandle.readline()#ITEM: NUMBER OF ATOMS
    N=int(filehandle.readline())
    buf=filehandle.readline()#ITEM: BOX BOUNDS pp pp pp
    xlims=np.array(filehandle.readline().split()).astype(float)
    ylims=np.array(filehandle.readline().split()).astype(float)
    zlims=np.array(filehandle.readline().split()).astype(float)
    buf=filehandle.readline()#ITEM: ATOMS id type xs ys zs
    scaled=buf.split()[-1]=='zs'
    table=[]
    for p in xrange(N): 
        buf=filehandle.readline()
        # print buf.split()
        table.append(np.array(buf.split(),dtype=float))
    
    table=np.array(table)
    if scaled:
        table[:,-3]*=xlims[1]-xlims[0]
        table[:,-2]*=ylims[1]-ylims[0]
        table[:,-1]*=zlims[1]-zlims[0]
    F=frame(table,timeframe=timeframe,boxinfo=[xlims, ylims, zlims])
    return F,buf


def array_to_xyz(array, filename, mode='w'):
    with open(filename, mode) as fw:
        fw.write("%d\nAtoms\n"%array.shape[0])
        for p in xrange(array.shape[0]):
            fw.write("A %g %g %g\n"%(array[p,0],array[p,1],array[p,2]))

def read_table_from_line(pathtofile,fromline,nlines, Type=float):
    """ 
    Read a ASCII table from file starting from a given line.
    """
    import numpy as np
    with open(pathtofile, 'r') as fin:
        data=[]
        for i in xrange(fromline):
            fin.readline()
        for i in xrange(nlines):
            line=fin.readline()
            data.append(np.array(line.split()).astype(Type))

        ncols=len(data[0])
        return np.array(data).reshape((nlines,ncols))

def save_obj(obj, name ):
    """Save dictionaries to file."""
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    """Load dictionaries from file."""
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)