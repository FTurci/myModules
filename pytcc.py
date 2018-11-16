import subprocess , os, inspect, glob


def to_xyz(coords, filename="tmp.xyz", types= None, mode="a"):
  N = coords.shape[0]
  if types ==None:
    types = ["A" for i in range(N)]
  with open(filename, mode) as fout:
    fout.write("%d\nAtoms\n"%N)
    for i in range(N):
      fout.write("%s %g %g %g \n"%(types[i], coords[i,0], coords[i,1],coords[i,2]))
      

def inputparameters(xyzFileName,boxType,desiredFrames,totalFrames, N,Na,rho, rcut,fc,boxName="box.txt", startFrom=0, frequency=1,bondType=1,pbc=True,bonds=False,clust=False, raw=False,_11a=False, _13a=False, pop=True):
	if len(rcut)>1:
		replacement=(boxType, boxName,xyzFileName, desiredFrames,totalFrames,N,Na, rho, totalFrames,startFrom,frequency,rcut[0], rcut[1], rcut[2],bondType, pbc,fc,bonds,clust,raw,_11a,_13a,pop)
	else:
		replacement=(boxType, boxName,xyzFileName, desiredFrames,totalFrames,N,Na, rho, totalFrames,startFrom,frequency,rcut[0], rcut[0], rcut[0],bondType, pbc,fc,bonds,clust,raw,_11a,_13a,pop)
 
	script="""[Box]	
; Specifies how to read simulation box
box_type			= %d				; 0 if cubic NVT, 1 if system non-cubic NVT, 2 if system is NPT, 3 triclinc with tilt (INTEGER)
box_name			= %s		; name of parameters file for box size (STRING)

[Run]	
; Run specific settings - these depend on your xyz file
xyzfilename			= %s		; File name of the xyz file to be analysed. (STRING)
frames				= %d				; FRAMES - frames to read from input xmol file (INTEGER)
totalframes			= %d
num_particles		= %d			; Total number of particles. (INTEGER)
numA_particles		= %d			; Number of type A particles (same as num particles if not binary) (INTEGER)
number_density		= %g			; Number of particles per unit volume (DOUBLE)
simulationstarttime = 0				; These values have no effect on the simulatin, they only serve to label the frames in the output files. (DOUBLE)
simulationtimestep	= 1			; These values have no effect on the simulation, they only serve to label the frames in the output files.(DOUBLE)
simulationendtime	= %d				; These values have no effect on the simulation, they only serve to label the frames in the output files.(DOUBLE)
start_from			= %d				; start reading from this frame in the xmol file (INTEGER)
sample_freqency		= %d				; frequency at which to take frames from the xmol file (INTEGER)

[Simulation]	
; Simulation specific settings - these depend on the type of system you are analysing
rcutAA				= %g	; A-A bond lengths (for simple bond detection) (DOUBLE)
rcutAB				= %g	; A-B bond lengths (for simple bond detection) (DOUBLE)
rcutBB				= %g	; B-B bond lengths (for simple bond detection) (DOUBLE)
bond_type			= %d		; 0 simple bond length, 1 Voronoi bond detection (BINARY INTEGER)
PBCs				= %d     ; 0 off, 1 on, Use period boundary conditions (BINARY INTEGER)
voronoi_parameter	= %g  ; Modified Voronoi Fc parameter (DOUBLE from 0 to 1)
num_bonds			= 30	; max number of bonds to one particle (INTEGER)
cell_list			= 0		; use Cell List to calculate bond network (and potential if used as well) (BINARY INTEGER)
potential_type		= 0		; 0 BLJ, 1 SFBLJ, 2 MorYuk: polydisp morse+yukawa, 3 not used, 4 IPL, 5 BLJ_WCAs, 6 SFIPL, 7 CRVT (INTEGER)

[Output]		
; Determines what the TCC will output
bonds 				= %d		; write out bonds file (BINARY INTEGER)
clusts 				= %d		; write clusts_** files containing all clusters - USES LOTS OF HDD SPACE (BINARY INTEGER)
raw 				= %d		; write raw_** xmol cluster files (BINARY INTEGER)
11a 				= %d		; write centres of 11A (BINARY INTEGER)
13a 				= %d		; write centres of 13A (BINARY INTEGER)
pop_per_frame 		= %d		; write particle fraction of each cluster per frame (BINARY INTEGER)
bin_width 			= 0.02	; bin width for bond length distributions (double)
bond_length 		= 0		; write bond length distributions (BINARY INTEGER)
bond_length_cluster	= 0		; write bond length distributions for each cluster type (BINARY INTEGER)
bond_length_dev 	= 0		; write bond length deviations from ground.state.bondlengths.dat for each cluster type (BINARY INTEGER)
neighbour_dist 		= 0		; write number of neighbour distributions  (BINARY INTEGER)
bonded_dist 		= 0		; write distributions for the number of particles bonded to the centre of each cluster (BINARY INTEGER)
cluster_composition	= 0		; write compositions of each cluster in terms of A and B species (BINARY INTEGER)
subclusters			= 0		; write subclusters of each cluster, if dynamics also done on required subcluster (BINARY INTEGER)

[Extra]		
; Special settings for extra functions
potential_energy 	= 0    ; do potential energy calculations (BINARY INTEGER)
coslovich			= 0    ; do Coslovich-style Voronoi faces analysis (BINARY INTEGER)
dodynamics			= 0    ; do Dynamics Analysis - choose which clusters and set memory sizes in static.memsize.dat (BINARY INTEGER)
alpha_time 			= 1.0  ; alpha relaxtion time (in simulation time units)(DOUBLE)
debug 				= 1    ; printing running (per frame) debug information (BINARY INTEGER)
shear				= 0    ; shear amount (for Lees-Edwards BCs) (DOUBLE)

; Potential parameters are in potentialparams.in
"""%replacement
	with open('inputparameters.ini', 'w') as fw:
		fw.write(script)
	

def prepare(fileName, rho):
	Na=0	
	import numpy as np
	if fileName[-3:]=="xyz":
		with open(fileName, 'r') as fin:
			N=int(fin.readline())	
			box=[]	
			fin.readline()	
			coords=[]
			for p in xrange(N):
				line=fin.readline().split()
				if line[0]=='A':
					Na+=1
				coords.append(np.asarray(line[1:]).astype(float) )
			coords=np.array(coords)
			mins=coords.min(axis=0)
			maxs=coords.max(axis=0)
			box=maxs-mins


	else:
		# assume atom format:
		with open(fileName, 'r') as fin:
			fin.readline()
			fin.readline()
			fin.readline()
			N=int(fin.readline())
			fin.readline()
			lims=fin.readline().split()
			Lx=float(lims[1])-float(lims[0])
			lims=fin.readline().split()
			Ly=float(lims[1])-float(lims[0])
			lims=fin.readline().split()
			Lz=float(lims[1])-float(lims[0])
			box=[Lx,Ly,Lz]
			fin.readline()
			for p in xrange(N):
				line=fin.readline().split()
				if line[1]=='A' or line[1]==1:
					Na+=1
	return N,Na, box

class TCC(object):
	"""Construct a TCC classifier for TCC. It executes the (precompiled) TCC code and returns the statistics per cluster."""
	def __init__(self, fileName, totalFrames,desiredFrames='all',rcut=2.,fc=1,rho=None,boxType=0, box=None,startFrom=0, frequency=1,bondType=1,pbc=True,bonds=False, clust=False,raw=False,_11a=False, _13a=False, pop=False):
		super(TCC, self).__init__()


		assert ((boxType<2) or boxType==3), "!!! Handling only boxTypes 0, 1 and 3"

		if boxType==0 and rho==None:
			print "!!! boxType=0 ==> Cubic, but rho set to "+str(rho)
			quit()

		self.filename = fileName
		self.boxType=boxType
		self.N, self.Na,self.box=prepare(fileName,rho)
		
		if boxType==1:
			if len(self.box)==0:
				print "!!! boxtype=1 ==> NonCubic but box is ",self.box
				quit()
			with open("box.txt",'w') as fw:
				fw.write("#iter Lx Ly Lz\n")
				fw.write("0 %g %g %g\n"%(self.box[0],self.box[1],self.box[2]))	
		if boxType==3:
			self.box = box
			if len(self.box)==0:
				print "!!! boxtype=1 ==> NonCubic but box is ",self.box
				quit()
			with open("box.txt",'w') as fw:
				fw.write("#iter Lx Ly Lz tilt\n")
				fw.write("0 %g %g %g %g\n"%(self.box[0],self.box[1],self.box[2], self.box[3]))	

		if rho==None or rho=='guess':
			self.rho=self.N/self.box[0]/self.box[1]/self.box[2]
		else:
			self.rho=rho


		if fileName[-3:]=="xyz":
			self.xyzFileName=fileName
		if desiredFrames=="all":
			desiredFrames=totalFrames

		inputparameters(self.filename,self.boxType,desiredFrames,totalFrames, self.N,self.Na,self.rho, [rcut],fc,boxName="box.txt", startFrom=0, frequency=1,bondType=1,pbc=pbc,bonds=bonds,clust=clust, raw=raw,_11a=_11a, _13a=_13a, pop=pop)
		
		output=subprocess.check_output("~/bin/TCC", shell=True)

		with open("tcc.log",'w') as fw:
			fw.write(output)
		
		# reading the results
		self.clustNumber={}
		self.gross={}
		self.net={}
		self.pop={}
		with open(glob.glob(fileName+"*fc%g*PBCs%d*static_clust"%(fc,pbc))[0], 'r') as fin:
			fin.readline()
			fin.readline()
			count=0

			for line in fin: 
				# print count, line

				cluster=line.split()
				
				key=cluster[0].split('_')[0]
				self.clustNumber[key]=int(cluster[1])
				self.gross[key]=int(cluster[2])
				self.net[key]=int(cluster[3])
				self.pop[key]=float(cluster[4])
				count+=1
				if count> 42:
					break

#classifier = TCC("/Users/ft14968/Dropbox/Francesco/StudentsProjects/WritingConnection/fig/snapshots/phi001-be10.xyz",1, rho='guess', fc=0.82)

# print classifier.filename, classifier.boxType, classifier.box, classifier.Na, classifier.rho
		

