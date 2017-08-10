import pylab as pl
from periodic_kdtree import PeriodicCKDTree
from scipy.spatial import KDTree, cKDTree
def get_marked(xyz,labels,box, marker=True, rcut=1.4,periodic=False):
	select=xyz[labels==marker]
	# print select
	if periodic:
		T = PeriodicCKDTree(box,select)
	else:
		T=cKDTree(select)
	# Find neighbors within a fixed distance of a point
	balls = T.query_ball_point(select, r=rcut)

	visited=pl.zeros(select.shape[0])
	added=pl.zeros(select.shape[0])
	clusters = []

	def addballs(p,cluster):
	        if visited[p]==0:
	                visited[p]=1
	                cluster.append(p)
	                for e in balls[p]:
	                        addballs(e,cluster)

	for i in xrange(select.shape[0]):
	        cluster=[]
	        addballs(i,cluster)
	        if len(cluster)>0:
	                clusters.append(cluster)
	return clusters

def get_tcc(configuration, tccrawfile,box,rcut=1.,criterium="not marked"):
	"Get the connected cluster formed by  marked or not marked particles"
	xyz=pl.loadtxt(configuration,skiprows=2,usecols=[1,2,3])
	cl=pl.loadtxt(tccrawfile, skiprows=3,dtype="S1")
	
	if criterium=='not marked':
		select=xyz[(cl=='A')+(cl=='B')]
	if criterium=="marked":
		 select=xyz[(cl=='C')+(cl=='D')]

	T = PeriodicCKDTree(box,select)
	# Find neighbors within a fixed distance of a point
	balls = T.query_ball_point(select, r=rcut)

	visited=pl.zeros(select.shape[0])
	added=pl.zeros(select.shape[0])
	clusters = []

	def addballs(p,cluster):
		if visited[p]==0:
			visited[p]=1
			cluster.append(p)
			for e in balls[p]:
				addballs(e,cluster)

	for i in xrange(select.shape[0]):
		cluster=[]
		addballs(i,cluster)
		if len(cluster)>0:
			clusters.append(cluster)
	return clusters

