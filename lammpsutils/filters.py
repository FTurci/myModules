import numpy as np

def read_log(filename, begin="Step", end="Loop"):
	read = False
	data = []
	with open(filename) as fin:
		for line in fin:
			split = line.split()
			if len(split)>0 and split[0]==begin:
				read = True
				header = split
				continue
			if len(split)>0 and split[0]==end:
				read = False
				continue	
			if read :
				data.append(np.array(split).astype(float))

	data = np.array(data)
	dictionary = {}
	for i,name in enumerate(header):
		dictionary[name]=data[:,i]
	return dictionary



    
