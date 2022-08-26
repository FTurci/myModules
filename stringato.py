import re
def extract_floats(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string)

def stampa(*args,**kwargs):
	print("::",*args,**kwargs)


def num_word(word,words):
    """Returns the first number after a certain word in a string of words."""
    # try float
    floats = re.compile(f'{word}'+'(\d+\.\d+)').findall(words)
    if len(floats[0])==0:
        ints = re.compile(f'{word}'+'(\d+)').findall(words)
        return float(ints[0])
    else:
        return float(floats[0])
