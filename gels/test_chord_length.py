import numpy as np 
from chord_length import ChordLengthAnalyser

print("\n# One dimensional example\n")
a= np.random.randint(0,2,size=10)

print("The input matrix is")
print(a)
cl = ChordLengthAnalyser(a)
cl.compute(warning=False)

print("\nThe chord lengths are")
print(cl.lengths)

print("\n# Two dimensional example\n")
a= np.zeros((5,5))
a[1:3]=1

print("The input matrix is")
print(a)
cl = ChordLengthAnalyser(a)
cl.compute(warning=False)

print("\nThe chord lengths are")
print(cl.lengths)


print("\n# Three dimensional example\n")
a= np.zeros((5,5,5))
a[1:3]=1
a[3,0,0]=1 #add a bit to see if it is picked up

print("The input matrix is")
print(a)
cl = ChordLengthAnalyser(a)
cl.compute(warning=False)

print("\nThe chord lengths are")
print(cl.lengths)