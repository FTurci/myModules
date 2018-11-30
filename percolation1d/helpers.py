import numpy as np
from numpy.random import uniform
from scipy import ndimage

def generate(N,p):
	a = []
	for i in range(N):
		if uniform()<=p:
			a.append(1)
		else:
			a.append(0)
	return a

def get_lengths(a):
	"Get cluster sizes using ndimage"
	labels,numlabels = ndimage.label(a)
	return ndimage.sum(a,labels, np.arange(1, numlabels+1))

def clusters(a):
	N = len(a)
	clusters=[]
	i=0
	while i <N:
		clen =0
		while i<N and a[i]==1 :
			i+=1
			clen +=1
		if clen!=0:
			clusters.append(clen)
		i+=1

	return np.array(clusters)


def average_clusters(plow,phigh,num=10, N=10000):
	ps = np.linspace(plow,phigh,num)
	means = []
	for p in ps:
		a = generate(N,p)
		c = clusters(a)
		means.append(c.mean())
	import pylab as pl
	pl.plot(ps, means,'-o')
	pl.plot(ps,theory_mean(ps))
	pl.show()

def n(p, N=10000):
	a = generate(N,p)
	c = clusters(a)
	b = np.bincount(c)
	s=np.arange(len(b))
	return s,b*1./len(a),c

def get_cluster(a,k):
	"Determine cluster containing the occupied site k"
	#forward
	i=k
	count =0
	while a[i]==1:
		count+=1
		i+=1
	#backward
	i=k-1
	while a[i] ==1:
		count+=1
		i-=1
	return count

def average(a, iterations=1000):
	# Average size of a cluster that we 
	# hit if we point to an arbitrary 
	# occupied site that is not a part of 
	# an infinite cluster
	a=np.array(a)
	ids = np.arange(len(a))
	occupied = ids[a==1]

	clusts=[] 
	for i in range(iterations):
		k=np.random.choice(occupied)
		if a[k]==1:
			clusts.append(get_cluster(a,k))
		else:
			print ("Something is wrong")
	return np.mean(clusts)

def theory_mean(p, ratio=False):
	#we compute sum(s n_s)/sum(n_s)
	#given that n_s(p)=(1-p)^2*p^sum
	# which is 1/(1-p) when considering the geometric series properties

	if ratio: #explicit ratio
		normalization = (1-p)**2*(1/(1-p)-1)
		numerator = (1-p)**2*p/(p-1)**2
		return numerator/normalization
	else:
		return 1./(1-p)

def theory_average(p):
	# Average cluster size that we hit if we pick an occcupied site at random
	# Notice: this is != for the mean cluster size (see above)
	return (1+p)/(1-p)