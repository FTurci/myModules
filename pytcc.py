import subprocess , os, inspect, glob


def to_xyz(coords, filename="tmp.xyz", types= None, mode="a"):
  N = coords.shape[0]
  if types ==None:
    types = ["A" for i in range(N)]
  with open(filename, mode) as fout:
    fout.write("%d\nAtoms\n"%N)
    for i in range(N):
      fout.write("%s %g %g %g \n"%(types[i], coords[i,0], coords[i,1],coords[i,2]))
      

def inputparameters(xyzFileName, frames, frequency,boxType=1, boxName="box.txt", rcut=2.0,min_cutAA=0.0,bondType=1, pbc=True,fc=0.82,num_bonds=50,cell_list=False,all_clust=True,write_bonds=False,write_clusts=False,write_raw=False,write_xyz=False,write_pop=True):
	if len(rcut)>1:
		replacement=(boxType, boxName,xyzFileName, frames, frequency, rcut[0], rcut[1], rcut[2],min_cutAA,bondType, pbc,fc,num_bonds,cell_list,all_clust,write_bonds,write_clusts,write_raw,write_xyz,write_pop)
	else:
		replacement=(boxType, boxName,xyzFileName, frames, frequency, rcut[0], rcut[0], rcut[0],min_cutAA,bondType, pbc,fc,num_bonds,cell_list,all_clust,write_bonds,write_clusts,write_raw,write_xyz,write_pop)
 
	script="""[[Box]
; Specifies how to read simulation box
box_type      = %d       ; 1 if NVT, 2 if NPT, 3 if triclinc with tilt
box_name      = %s   ; name of parameters file for box size

[Run]
; Run specific settings - these depend on your xyz file
xyzfilename     = %s ; File name of the xyz file to be analysed.
frames        = %d      ; Frames to process
sample_freqency   = %d       ; frequency at which to take frames from the xmol file

[Simulation]
; Simulation specific settings - these depend on the type of system you are analysing
rcutAA        = %f ; maximum A-A bond lengths  // The cutoff is always applied whether Voronoi bonds are used or not
rcutAB        = %f ; maximum A-B bond lengths
rcutBB        = %f ; maximum B-B bond lengths
min_cutAA           = %f   ; minimum A-A bond length. Good for excluding overlapping particles in ideal gases.
bond_type     = %d   ; 0 simple bond length, 1 Voronoi bond detection
PBCs        = %d     ; 0 Dont use periodic boundary conditions, 1 Use PBCs,
voronoi_parameter = %g   ; Modified Voronoi Fc parameter
num_bonds     = %d  ; max number of bonds to one particle
cell_list     =  %d  ; use Cell List to calculate bond network
analyse_all_clusters = %d    ; If set to zero, read clusters to analyse from clusters_to_analyse.ini

[Output]
; Determines what the TCC will output
bonds         = %d   ; write out bonds file
clusts        = %d   ; write clusts_** files containing all clusters - USES LOTS OF HDD SPACE
raw         = %d   ; write raw_** xmol cluster files
do_XYZ              = %d     ; write clusters to xyz files
11a         = 0   ; write centres of 11A
13a         = 0   ; write centres of 13A
pop_per_frame     = %d   ; write particle fraction of each cluster per frame"""%replacement
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
		

