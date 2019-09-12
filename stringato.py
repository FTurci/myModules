import re
def extract_floats(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string)

def stampa(*args,**kwargs):
	print("::",*args,**kwargs)
	

