import numpy as np
class filter_elem():
    x = np.empty(0,dtype=float)
    fobs = np.empty(0,dtype=float)
    alfa = 0.0
    alfa_c = np.empty(0,dtype=float)

    def __init__(self,x,f,alfa,alfa_c):
        self.x = np.array(x,dtype=float)
        self.fobs = np.array(f,dtype=float)
        self.alfa = np.double(alfa)
        self.alfa_c = np.array(alfa_c,dtype=float)

    def print(self):
        print(self.x,' ',self.fobs,' ',self.alfa,' ',self.alfa_c)

x = np.ones(4)
f = np.zeros(3)
alfa = 0.0
alfa_c = np.ones(4)

e = filter_elem(x,f,alfa,alfa_c)
print(e)
e.print()

A = np.empty(0,dtype=filter_elem)
print(A)
A = np.append(A,e)
A[0].print()

alfa = 0.5
e = filter_elem(x,f,alfa,alfa_c)
A = np.append(A,e)
A[1].print()
#A = np.delete(A,0)
#A[0].print()

print()
for i,e in enumerate(A):
    print(i,' ',end='')
    e.print()

alfa = 3.5
e = filter_elem(x,f,alfa,alfa_c)
A = np.insert(A,1,e)
print()
for i,e in enumerate(A):
    print(i,' ',end='')
    e.print()
