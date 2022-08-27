import re
def extract_floats(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string)

def stampa(*args,symbol=None,**kwargs):
    """Overloads print  function with a prepended symbol."""
<<<<<<< HEAD
    if symbol is None:
        symbol = kwargs.get('s')
        if symbol is None:
            symbol = ":::"
        else:
            del kwargs['s']
=======
>>>>>>> e7eb6d7c2a17f2687c25a6db19e8c710b386a966
    print(symbol,*args,**kwargs)


def num_word(word,words, verbose=False):
    """Returns the first number after a certain word in a string of words."""
    # try float
    floats = re.compile(f'{word}'+'(\d+\.\d+)').findall(words)
    if verbose:
        print("<num_word> detected:",word, floats)
    if len(floats)==0:
        ints = re.compile(f'{word}'+'(\d+)').findall(words)
        return int(ints[0])
    else:
        return float(floats[0])
