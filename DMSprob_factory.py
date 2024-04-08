import numpy as np

class problem_BK1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n

		self.q     = 2
		self.lb    = -5.0*np.ones(self.n)
		self.ub    = 10.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'BK1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		n = len(x)

		obj_1 = x[0]**2 + x[1]**2

		obj_2 = (x[0]-5.0)**2 + (x[1]-5.0)**2

		return np.array([obj_1, obj_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_CL1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 4
		self.n     = n

		self.q     = 2

		F = 10.0
		E = 2.e+5
		L = 200.0
		sigma = 10.0
		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		self.lb[0] = F/sigma; self.ub[0] = 3*F/sigma
		self.lb[1] = np.sqrt(2.0)*F/sigma; self.ub[1] = 3*F/sigma
		self.lb[2] = np.sqrt(2.0)*F/sigma; self.ub[3] = 3*F/sigma
		self.lb[3] = F/sigma; 	    self.ub[3] = 3*F/sigma

		x0 = self.startp()

		self.name  = 'CL1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return (self.lb+self.ub)/2.0

	def functs(self,x):

		F = 10.0
		E = 2.e+5
		L = 200.0
		sigma = 10.0

		ff_1 = L*(2.0*x[0] + np.sqrt(2.0)*x[1] + np.sqrt(x[2]) + x[3])
		ff_2 = F*L/E*(2.0/x[0] + (2.0*np.sqrt(2.0))/x[1] - (2.0*np.sqrt(2.0))/x[2] + 2.0/x[3])

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DEB41():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.lb[0] = 0.1

		x0 = self.startp()

		self.name  = 'DEB 4.1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.array([0.1, 0.0])

	def functs(self,x):

		ff_1 = x[0]
		ff_2 = ( 2.0-np.exp(-((x[1]-0.2)/0.004)**2)-0.8*np.exp(-((x[1]-0.6)/0.4)**2) )/x[0]

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DEB53():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DEB 5.3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		gamma = 1.0
		beta  = 1.0
		alpha = 4.0

		if x[1] <= 0.4:
			gx = 4.0-3.0*np.exp(-((x[1]-0.2)/0.02)**2)
		else:
			gx = 4.0-2.0*np.exp(-((x[1]-0.7)/0.2)**2)

		ff_1 = 1.0-np.exp(-4.0*x[0])*np.sin(5.0*np.pi*x[0])**4
		if ff_1 < beta*gx:
			h = (1.0-(ff_1/(beta*gx))**alpha)
		else:
			h = 0.0

		ff_2 = gx*h

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DEB512a():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DEB 5.1.2.(a)'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		beta = 1.0
		alpha = 0.25
		if x[1] <= 0.4:
			gx = 4.0-3.0*np.exp(-((x[1]-0.2)/0.02)**2)
		else:
			gx = 4.0-2.0*np.exp(-((x[1]-0.7)/0.2)**2)

		ff_1 = 4.0*x[0]

		if (ff_1 <= beta*gx):
			h = (1.0-(ff_1/(beta*gx))**alpha)
		else:
			h = 0.0

		ff_2 = gx*h

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DEB512b():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DEB 5.1.2.(b)'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		beta = 1.0
		alpha = 4.0
		if x[1] <= 0.4:
			gx = 4.0-3.0*np.exp(-((x[1]-0.2)/0.02)**2)
		else:
			gx = 4.0-2.0*np.exp(-((x[1]-0.7)/0.2)**2)

		ff_1 = 4.0*x[0]

		if (ff_1 <= beta*gx):
			h = (1.0-(ff_1/(beta*gx))**alpha)
		else:
			h = 0.0

		ff_2 = gx*h

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DEB512c():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DEB 5.1.2.(c)'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		beta = 1.0

		if x[1] <= 0.4:
			gx = 4.0-3.0*np.exp(-((x[1]-0.2)/0.02)**2)
		else:
			gx = 4.0-2.0*np.exp(-((x[1]-0.7)/0.2)**2)

		ff_1 = 4.0*x[0]
		alpha = 0.25 + 3.75*(gx-1.0)
		if (ff_1 <= beta*gx):
			h = (1.0-(ff_1/(beta*gx))**alpha)
		else:
			h = 0.0

		ff_2 = gx*h

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DEB513():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DEB 5.1.3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		beta = 1.0
		alpha = 2.0
		q = 4.0

		gx = 1.0+10.0*x[1]

		ff_1 = x[0]

		h = 1.0-(ff_1/gx)**alpha-(ff_1/gx)*np.sin(2.0*np.pi*q*ff_1)

		ff_2 = gx*h

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DEB521a():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DEB 5.2.1.(a)'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		gamma = 0.25

		gx = 1.0+x[1]**gamma

		ff_1 = x[0]

		h = 1.0-(ff_1/gx)**2

		ff_2 = gx*h

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DEB521b():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DEB 5.2.1.(b)'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		gamma = 1.0

		gx = 1.0+x[1]**gamma

		ff_1 = x[0]

		h = 1.0-(ff_1/gx)**2

		ff_2 = gx*h

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DG01():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 1
		self.n     = n
		self.q     = 2

		self.lb    = -10.0*np.ones(self.n)
		self.ub    =  13.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DG01'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		ff_1 = np.sin(x[0])
		ff_2 = np.sin(x[0]+0.7)

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DPAM1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 2

		self.lb    = -0.3*np.ones(self.n)
		self.ub    =  0.3*np.ones(self.n)

		self.A = np.array([[0.218418,    -0.620254,     0.843784,    0.914311,    -0.788548,    0.428212,     \
			    0.103064,    -0.47373,     -0.300792,     -0.185507], \
			    [0.330423,     0.151614,     0.884043,   -0.272951,   -0.993822,    0.511197,      \
			   -0.0997948,   -0.659756,     0.575496,      0.675617], \
			    [0.180332,    -0.593814,    -0.492722,    0.0646786,   -0.666503,   -0.945716,     \
			   -0.334582,     0.611894,     0.281032,      0.508749], \
			   [-0.0265389,   -0.920133,     0.308861,   -0.0437502,   -0.374203,    0.207359,     \
			   -0.219433,     0.914104,     0.184408,      0.520599], \
			   [-0.88565,     -0.375906,    -0.708948,   -0.37902,     0.576578,    0.0194674,     \
			   -0.470262,     0.572576,     0.351245,     -0.480477], \
			    [0.238261,    -0.1596,      -0.827302,    0.669248,     0.494475,    0.691715,     \
			   -0.198585,     0.0492812,    0.959669,      0.884086], \
			   [-0.218632,    -0.865161,    -0.715997,    0.220772,     0.692356,    0.646453,     \
			   -0.401724,     0.615443,    -0.0601957,    -0.748176], \
			   [-0.207987,    -0.865931,     0.613732,   -0.525712,    -0.995728,    0.389633,     \
			   -0.064173,     0.662131,    -0.707048,     -0.340423], \
			    [0.60624,      0.0951648,   -0.160446,   -0.394585,    -0.167581,    0.0679849,    \
			    0.449799,     0.733505,    -0.00918638,    0.00446808], \
			    [0.404396,     0.449996,     0.162711,    0.294454,    -0.563345,   -0.114993,     \
			    0.549589,    -0.775141,     0.677726,      0.610715]])

		x0 = self.startp()

		self.name  = 'DPAM1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		A = self.A
		y  = np.matmul(A,x)

		gx = 1.0 + 10.0*float(self.n-1)
		for i in range(1,self.n):  #i = 2,n
			gx += (y[i]**2-10.0*np.cos(4.0*np.pi*y[i]))

		ff_1 = y[0]
		ff_2 = gx*np.exp(-y[0]/gx)

		return np.array([ff_1,ff_2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 7
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += ((x[i]-0.5)**2-np.cos(20.0*np.pi*(x[i]-0.5)))

		gx = 100.0*(k + gx)

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(self.q-1): #i = 1,qq-1
			f[0] *= x[i]

		f[0] = 0.5*(1.0+gx)*f[0]

		for i in range(1,self.q): #i = 2,qq
			f[i] = 1.0
			for j in range(self.q-i): #j = 1,qq-i
				f[i] *= x[j]
			f[i] = 0.5*(1.0+gx)*f[i]*(1.0-x[self.q-i])

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ1n2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ1 n=2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += ((x[i]-0.5)**2-np.cos(20.0*np.pi*(x[i]-0.5)))

		gx = 100.0*(k + gx)

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(self.q-1): #i = 1,qq-1
			f[0] *= x[i]

		f[0] = 0.5*(1.0+gx)*f[0]

		for i in range(1,self.q): #i = 2,qq
			f[i] = 1.0
			for j in range(self.q-i): #j = 1,qq-i
				f[i] *= x[j]
			f[i] = 0.5*(1.0+gx)*f[i]*(1.0-x[self.q-i])

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 12
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += (x[i]-0.5)**2

		gx = 100.0*(k + gx)

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(self.q-1): #i = 1,qq-1
			f[0] *= np.cos(0.5*np.pi*x[i])

		f[0] = (1.0+gx)*f[0]

		for i in range(1,self.q): #i = 2,qq
			f[i] = 1.0
			for j in range(self.q-i): #j = 1,qq-i
				f[i] *= np.cos(0.5*np.pi*x[j])
			f[i] = (1.0+gx)*f[i]*(np.sin(0.5*np.pi*x[self.q-i]))

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ2n2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ2 n=2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += (x[i]-0.5)**2

		gx = 100.0*(k + gx)

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(self.q-1): #i = 1,qq-1
			f[0] *= np.cos(0.5*np.pi*x[i])

		f[0] = (1.0+gx)*f[0]

		for i in range(1,self.q): #i = 2,qq
			f[i] = 1.0
			for j in range(self.q-i): #j = 1,qq-i
				f[i] *= np.cos(0.5*np.pi*x[j])
			f[i] = (1.0+gx)*f[i]*(np.sin(0.5*np.pi*x[self.q-i]))

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 12
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += ((x[i]-0.5)**2-np.cos(20.0*np.pi*(x[i]-0.5)))

		gx = 100.0*(k + gx)

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(self.q-1): #i = 1,qq-1
			f[0] *= np.cos(0.5*np.pi*x[i])

		f[0] = (1.0+gx)*f[0]

		for i in range(1,self.q): #i = 2,qq
			f[i] = 1.0
			for j in range(self.q-i): #j = 1,qq-i
				f[i] *= np.cos(0.5*np.pi*x[j])
			f[i] = (1.0+gx)*f[i]*(np.sin(0.5*np.pi*x[self.q-i]))

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ3n2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ3 n=2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += ((x[i]-0.5)**2-np.cos(20.0*np.pi*(x[i]-0.5)))

		gx = 100.0*(k + gx)

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(self.q-1): #i = 1,qq-1
			f[0] *= np.cos(0.5*np.pi*x[i])

		f[0] = (1.0+gx)*f[0]

		for i in range(1,self.q): #i = 2,qq
			f[i] = 1.0
			for j in range(self.q-i): #j = 1,qq-i
				f[i] *= np.cos(0.5*np.pi*x[j])
			f[i] = (1.0+gx)*f[i]*(np.sin(0.5*np.pi*x[self.q-i]))

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 12
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		alpha = 100.0
		k  = float(self.n-self.q+1)

		y = np.array([x[i]**alpha for i in range(self.n)])

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += ((y[i]-0.5)**2)

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(self.q-1): #i = 1,qq-1
			f[0] *= np.cos(0.5*np.pi*y[i])

		f[0] = (1.0+gx)*f[0]

		for i in range(1,self.q): #i = 2,qq
			f[i] = 1.0
			for j in range(self.q-i): #j = 1,qq-i
				f[i] *= np.cos(0.5*np.pi*y[j])
			f[i] = (1.0+gx)*f[i]*(np.sin(0.5*np.pi*y[self.q-i]))

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ4n2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ4 n=2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		alpha = 100.0
		k  = float(self.n-self.q+1)

		y = np.array([x[i]**alpha for i in range(self.n)])

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += ((y[i]-0.5)**2)

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(self.q-1): #i = 1,qq-1
			f[0] *= np.cos(0.5*np.pi*y[i])

		f[0] = (1.0+gx)*f[0]

		for i in range(1,self.q): #i = 2,qq
			f[i] = 1.0
			for j in range(self.q-i): #j = 1,qq-i
				f[i] *= np.cos(0.5*np.pi*y[j])
			f[i] = (1.0+gx)*f[i]*(np.sin(0.5*np.pi*y[self.q-i]))

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ5():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 12
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ5'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		alpha = 100.0
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += x[i]**0.1

		theta = np.zeros(self.n)
		for i in range(1,self.n):
			theta[i] = (np.pi/2.0)*(1.0+2.0*gx*x[i])/(2.0*(1.0+gx))

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(1,self.q-1): #i = 1,qq-1
			f[0] *= np.cos(theta[i])

		f[0] = (1.0+gx)*(np.cos(0.5*np.pi*x[0]))*f[0]

		for i in range(1,self.q-1): #i = 2,qq
			f[i] = 1.0
			for j in range(1,self.q-i): #j = 1,qq-i
				f[i] *= np.cos(theta[j])
			f[i] = (1.0+gx)*(np.cos(0.5*np.pi*x[0]))*f[i]*(np.sin(theta[self.q-i]))

		f[self.q-1] = (1.0+gx)*np.sin(0.5*np.pi*x[0])
		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ5n2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ5 n=2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		alpha = 100.0
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += x[i]**0.1

		theta = np.zeros(self.n)
		for i in range(1,self.n):
			theta[i] = (np.pi/2.0)*(1.0+2.0*gx*x[i])/(2.0*(1.0+gx))

		f = np.zeros(self.q)
		f[0] = 1.0
		for i in range(1,self.q-1): #i = 1,qq-1
			f[0] *= np.cos(theta[i])

		f[0] = (1.0+gx)*(np.cos(0.5*np.pi*x[0]))*f[0]

		for i in range(1,self.q-1): #i = 2,qq
			f[i] = 1.0
			for j in range(1,self.q-i): #j = 1,qq-i
				f[i] *= np.cos(theta[j])
			f[i] = (1.0+gx)*(np.cos(0.5*np.pi*x[0]))*f[i]*(np.sin(theta[self.q-i]))

		f[self.q-1] = (1.0+gx)*np.sin(0.5*np.pi*x[0])
		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ6():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 22
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ6'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		alpha = 100.0
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += x[i]
		gx = 1.0 + (9.0/k)*gx

		f = np.zeros(self.q)
		for i in range(self.q-1):
			f[i] = x[i]

		f[self.q-1] = 0.0
		for i in range(self.q-1):
			f[self.q-1] += (x[i]/(1.0+gx)*(1.0+np.sin(3.0*np.pi*x[i])))

		f[self.q-1] = (1.0+gx) * (float(self.q) - f[self.q-1])

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_DTLZ6n2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'DTLZ6 n=2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):
		k  = float(self.n-self.q+1)

		gx = 0.0
		for i in range(self.q-1,self.n): #i = qq,n
			gx += x[i]
		gx = 1.0 + (9.0/k)*gx

		f = np.zeros(self.q)
		for i in range(self.q-1):
			f[i] = x[i]

		f[self.q-1] = 0.0
		for i in range(self.q-1):
			f[self.q-1] += (x[i]/(1.0+gx)*(1.0+np.sin(3.0*np.pi*x[i])))

		f[self.q-1] = (1.0+gx) * (float(self.q) - f[self.q-1])

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_EX005():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.lb[0] = -1.0
		self.lb[1] =  1.0
		self.ub    = 2.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'EX. 005'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return (self.lb+self.ub)/2.0

	def functs(self,x):

		f1 = x[0]**2 - x[1]**2
		f2 = x[0]/x[1]

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_FAR1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -1.0*np.ones(self.n)
		self.ub    =      np.ones(self.n)

		x0 = self.startp()

		self.name  = 'FAR 1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = -2.0*np.exp(15.0*(-(x[0]-0.1)**2-x[1]**2))		\
			     -np.exp(20.0*(-(x[0]-0.6)**2-(x[1]-0.6)**2))	\
			     +np.exp(20.0*(-(x[0]+0.6)**2-(x[1]-0.6)**2))	\
			     +np.exp(20.0*(-(x[0]-0.6)**2-(x[1]+0.6)**2))	\
			     +np.exp(20.0*(-(x[0]+0.6)**2-(x[1]+0.6)**2))

		f2 = 2.0*np.exp(20.0*(-x[0]**2-x[1]**2))			\
			    +np.exp(20.0*(-(x[0]-0.4)**2-(x[1]-0.6)**2))	\
			    -np.exp(20.0*(-(x[0]+0.5)**2-(x[1]-0.7)**2))	\
			    -np.exp(20.0*(-(x[0]-0.5)**2-(x[1]+0.7)**2))	\
			    +np.exp(20.0*(-(x[0]+0.4)**2-(x[1]+0.8)**2))

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_FES1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'FES 1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.0
		f2 = 0.0
		for i in range(self.n): # i = 1,n
			f1 += (np.abs(x[i]-np.exp((float(i+1)/float(self.n))**2)/3.0)**0.5)
			f2 += ((x[i]-0.5*np.cos(10.0*np.pi*(float(i+1)/float(self.n)))-0.5)**2)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_FES2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'FES 2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.0
		f2 = 0.0
		f3 = 0.0
		for i in range(self.n): # i = 1,n
			f1 += ((x[i]-0.5*np.cos(10.0*np.pi*(float(i+1)/float(self.n)))-0.5)**2)
			f2 += (np.abs(x[i]-np.sin(float(i))**2*np.cos(float(i))**2)**0.5)
			f3 += (np.abs(x[i]-0.25*np.cos(float(i))*np.cos(float(2*(i+1)-2))-0.5)**0.5)


		return np.array([f1, f2, f3])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_FES3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 4

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'FES 3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.0
		f2 = 0.0
		f3 = 0.0
		f4 = 0.0
		for i in range(self.n): # i = 1,n
			f1 += (np.abs(x[i]-np.exp((float(i+1)/float(self.n))**2)/3.0)**0.5)
			f2 += (np.abs(x[i]-np.sin(float(i))**2*np.cos(float(i))**2)**0.5)
			f3 += (np.abs(x[i]-0.25*np.cos(float(i))*np.cos(float(2*(i+1)-2))-0.5)**0.5)
			f4 += ((x[i]-0.5*np.sin(1000.0*np.pi*(float(i+1)/float(self.n)))-0.5)**2)


		return np.array([f1, f2, f3, f4])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_FONSECA():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -4.0*np.ones(self.n)
		self.ub    =  4.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'Fonseca'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 1.0-np.exp(-(x[0]-1.0)**2-(x[1]+1.0)**2)
		f2 = 1.0-np.exp(-(x[0]+1.0)**2-(x[1]-1.0)**2)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_I1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'I1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		k = 4
		l = 4
		kr = 4.0
		lr = 4.0
		pi = np.pi
		pi2= np.pi/2.0
		S = np.ones(self.q)
		Av = np.ones(self.q-1)
		zmax = np.ones(self.n)
		y = np.array([z[i]/zmax[i] for i in range(self.n)])
		t1 = np.array([y[i] for i in range(self.n)])

		t2 = np.zeros(self.n)
		for i in range(k): # i = 1,k
			t2[i] = t1[i]
		for i in range(k,self.n): #i = k+1,n
			t2[i] = np.abs(t1[i]-0.35)/np.abs(np.floor(0.35-t1[i])+0.35)

		w = np.ones(self.n)
		t3 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			a = 0.0
			b = 0.0
			for j in range(2*i,2*(i+1)):  #j = ((i-1)*k/(M-1)+1), (i*k/(M-1))  2*(i-1)+1, 2*i
				a += (w[j]*t2[j])
				b += w[j]
			t3[i] = a/b

		a = 0.0
		b = 0.0
		for j in range(k,self.n): #j = k+1,n
			a += (w[j]*t2[j])
			b += w[j]
		t3[self.q-1] = a/b

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],Av[i])*(t3[i]-0.5)+0.5

		x[self.q-1] = t3[self.q-1]

		h = np.zeros(self.q)
		h[0] = 1.0
		for j in range(self.q-1): #j = 1,M-1
			h[0] *= np.sin(x[j]*pi2)

		for i in range(1,self.q-1): #i = 2,M-1
			h[i] = 1.0
			for j in range(self.q-i): #j = 1,M-i
				h[i] *= np.sin(x[j]*pi2)

			h[i] *= np.cos(x[self.q-i]*pi2)

		h[self.q-1] = np.cos(x[0]*pi2)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_I2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'I2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		k = 4
		l = 4
		kr = 4.0
		lr = 4.0
		AA = 0.98/49.98
		BB = 0.02
		CC = 50.0
		pi = np.pi
		pi2= np.pi/2.0
		S = np.ones(self.q)
		Av = np.ones(self.q-1)
		zmax = np.ones(self.n)
		y = np.array([z[i]/zmax[i] for i in range(self.n)])

		w = np.ones(self.n)
		r_sum = np.zeros(self.n-1)
		for i in range(self.n-1):
			a = 0.0
			b = 0.0
			for j in range(i+1,self.n):
				a += w[j]*y[j]
				b += w[j]
			r_sum[i] = a/b

		t1 = np.zeros(self.n)
		for i in range(self.n-1):
			t1[i] = y[i]**(BB+(CC-BB)*(AA-(1.0-2.0*r_sum[i])*np.abs(np.floor(0.5-r_sum[i])+AA)))
		t1[self.n-1] = y[self.n-1]

		t2 = np.zeros(self.n)
		for i in range(k): # i = 1,k
			t2[i] = t1[i]
		for i in range(k,self.n): #i = k+1,n
			t2[i] = np.abs(t1[i]-0.35)/np.abs(np.floor(0.35-t1[i])+0.35)

		t3 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			a = 0.0
			b = 0.0
			for j in range(2*i,2*(i+1)):  #j = ((i-1)*k/(M-1)+1), (i*k/(M-1))  2*(i-1)+1, 2*i
				a += (w[j]*t2[j])
				b += w[j]
			t3[i] = a/b

		a = 0.0
		b = 0.0
		for j in range(k,self.n): #j = k+1,n
			a += (w[j]*t2[j])
			b += w[j]
		t3[self.q-1] = a/b

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],Av[i])*(t3[i]-0.5)+0.5

		x[self.q-1] = t3[self.q-1]

		h = np.zeros(self.q)
		h[0] = 1.0
		for j in range(self.q-1): #j = 1,M-1
			h[0] *= np.sin(x[j]*pi2)

		for i in range(1,self.q-1): #i = 2,M-1
			h[i] = 1.0
			for j in range(self.q-i): #j = 1,M-i
				h[i] *= np.sin(x[j]*pi2)

			h[i] *= np.cos(x[self.q-i]*pi2)

		h[self.q-1] = np.cos(x[0]*pi2)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_I3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'I3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		k = 4
		l = 4
		kr = 4.0
		lr = 4.0
		AA = 0.98/49.98
		BB = 0.02
		CC = 50.0
		pi = np.pi
		pi2= np.pi/2.0
		S = np.ones(self.q)
		Av = np.ones(self.q-1)
		zmax = np.ones(self.n)
		y = np.array([z[i]/zmax[i] for i in range(self.n)])

		w = np.ones(self.n)
		r_sum = np.zeros(self.n)
		for i in range(1,self.n):
			a = 0.0
			b = 0.0
			for j in range(i):
				a += w[j]*y[j]
				b += w[j]
			r_sum[i] = a/b

		t1 = np.zeros(self.n)
		t1[0] = y[0]
		for i in range(1,self.n):
			t1[i] = y[i]**(BB+(CC-BB)*(AA-(1.0-2.0*r_sum[i])*np.abs(np.floor(0.5-r_sum[i])+AA)))
		t1[self.n-1] = y[self.n-1]

		t2 = np.zeros(self.n)
		for i in range(k): # i = 1,k
			t2[i] = t1[i]
		for i in range(k,self.n): #i = k+1,n
			t2[i] = np.abs(t1[i]-0.35)/np.abs(np.floor(0.35-t1[i])+0.35)

		t3 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			a = 0.0
			b = 0.0
			for j in range(2*i,2*(i+1)):  #j = ((i-1)*k/(M-1)+1), (i*k/(M-1))  2*(i-1)+1, 2*i
				a += (w[j]*t2[j])
				b += w[j]
			t3[i] = a/b

		a = 0.0
		b = 0.0
		for j in range(k,self.n): #j = k+1,n
			a += (w[j]*t2[j])
			b += w[j]
		t3[self.q-1] = a/b

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],Av[i])*(t3[i]-0.5)+0.5

		x[self.q-1] = t3[self.q-1]

		h = np.zeros(self.q)
		h[0] = 1.0
		for j in range(self.q-1): #j = 1,M-1
			h[0] *= np.sin(x[j]*pi2)

		for i in range(1,self.q-1): #i = 2,M-1
			h[i] = 1.0
			for j in range(self.q-i): #j = 1,M-i
				h[i] *= np.sin(x[j]*pi2)

			h[i] *= np.cos(x[self.q-i]*pi2)

		h[self.q-1] = np.cos(x[0]*pi2)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_I4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'I4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		k = 4
		l = 4
		kr = 4.0
		lr = 4.0

		pi = np.pi
		pi2= np.pi/2.0
		S = np.ones(self.q)
		Av = np.ones(self.q-1)
		zmax = np.ones(self.n)
		y = np.array([z[i]/zmax[i] for i in range(self.n)])

		w = np.ones(self.n)

		t1 = np.zeros(self.n)
		for i in range(self.n):
			t1[i] = y[i]

		t2 = np.zeros(self.n)
		for i in range(k): # i = 1,k
			t2[i] = t1[i]
		for i in range(k,self.n): #i = k+1,n
			t2[i] = np.abs(t1[i]-0.35)/np.abs(np.floor(0.35-t1[i])+0.35)

		t3 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			a = 0.0
			b = 0.0
			for ii in range(2*i,2*(i+1)):
				asum = 0.0
				for j in range(1): #range(k/(self.q-1)-1):
					ind = ((i)*k/(self.q-1)+1)+( (ii+j-((i)*k/(self.q-1)+1)+1) % \
						   (((i+1)*k/(self.q-1))-((i)*k/(self.q-1)+1)+1) ) - 1

					asum = asum + np.abs(t2[ii]-t2[int(ind)])
				a += t2[ii] + asum

			t3[i] = a/(((i*kr/float(self.q-1))-((i-1)*kr/float(self.q-1)+1.0)+1.0)/	\
			(kr/float(self.q-1))*np.ceil(kr/float(self.q-1)/2.0)*(1.0+2.0*kr/	\
		 	float(self.q-1)-2.0*np.ceil(kr/float(self.q-1)/2.0)))

		a = 0.0
		b = 0.0
		for ii in range(k,self.n): #j = k+1,n
			asum = 0.0
			for j in range(l-1):

				asum += np.abs(t2[ii]-t2[k+((ii+j-(k+1)+1) % (self.n-k))])
			a += t2[ii] + asum

		t3[self.q-1] = a/((float(self.n-k)/lr)*np.ceil(lr/2.0)*(1.0+2.0*lr-2.0*np.ceil(lr/2.0)))

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],Av[i])*(t3[i]-0.5)+0.5

		x[self.q-1] = t3[self.q-1]

		h = np.zeros(self.q)
		h[0] = 1.0
		for j in range(self.q-1): #j = 1,M-1
			h[0] *= np.sin(x[j]*pi2)

		for i in range(1,self.q-1): #i = 2,M-1
			h[i] = 1.0
			for j in range(self.q-i): #j = 1,M-i
				h[i] *= np.sin(x[j]*pi2)

			h[i] *= np.cos(x[self.q-i]*pi2)

		h[self.q-1] = np.cos(x[0]*pi2)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_I5():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'I5'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		k = 4
		l = 4
		kr = 4.0
		lr = 4.0
		AA = 0.98/49.98
		BB = 0.02
		CC = 50.0

		pi = np.pi
		pi2= np.pi/2.0
		S = np.ones(self.q)
		Av = np.ones(self.q-1)
		zmax = np.ones(self.n)
		y = np.array([z[i]/zmax[i] for i in range(self.n)])

		w = np.ones(self.n)

		t1 = np.zeros(self.n)
		r_sum = np.zeros(self.n)
		for i in range(1,self.n):
			a = 0.0
			b = 0.0
			for j in range(i): #j = 1,i-1
				a += (w[j]*y[j])
				b += w[j]

			r_sum[i] = a/b

		t1[0] = y[0]
		for i in range(1,self.n): #i = 2,n
			t1[i] = y[i]**(BB+(CC-BB)*(AA-(1.0-2.0*r_sum[i])*np.abs(np.floor(0.5-r_sum[i])+AA)))

		t2 = np.zeros(self.n)
		for i in range(k): # i = 1,k
			t2[i] = t1[i]
		for i in range(k,self.n): #i = k+1,n
			t2[i] = np.abs(t1[i]-0.35)/np.abs(np.floor(0.35-t1[i])+0.35)

		t3 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			a = 0.0
			b = 0.0
			for ii in range(2*i,2*(i+1)):
				asum = 0.0
				for j in range(1): #range(k/(self.q-1)-1):
					ind = ((i)*k/(self.q-1)+1)+( (ii+j-((i)*k/(self.q-1)+1)+1) % \
						   (((i+1)*k/(self.q-1))-((i)*k/(self.q-1)+1)+1) ) - 1

					asum = asum + np.abs(t2[ii]-t2[int(ind)])
				a += t2[ii] + asum

			t3[i] = a/(((i*kr/float(self.q-1))-((i-1)*kr/float(self.q-1)+1.0)+1.0)/	\
			(kr/float(self.q-1))*np.ceil(kr/float(self.q-1)/2.0)*(1.0+2.0*kr/	\
		 	float(self.q-1)-2.0*np.ceil(kr/float(self.q-1)/2.0)))

		a = 0.0
		b = 0.0
		for ii in range(k,self.n): #j = k+1,n
			asum = 0.0
			for j in range(l-1):

				asum += np.abs(t2[ii]-t2[k+((ii+j-(k+1)+1) % (self.n-k))])
			a += t2[ii] + asum

		t3[self.q-1] = a/((float(self.n-k)/lr)*np.ceil(lr/2.0)*(1.0+2.0*lr-2.0*np.ceil(lr/2.0)))

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],Av[i])*(t3[i]-0.5)+0.5

		x[self.q-1] = t3[self.q-1]

		h = np.zeros(self.q)
		h[0] = 1.0
		for j in range(self.q-1): #j = 1,M-1
			h[0] *= np.sin(x[j]*pi2)

		for i in range(1,self.q-1): #i = 2,M-1
			h[i] = 1.0
			for j in range(self.q-i): #j = 1,M-i
				h[i] *= np.sin(x[j]*pi2)

			h[i] *= np.cos(x[self.q-i]*pi2)

		h[self.q-1] = np.cos(x[0]*pi2)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_IKK1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 3

		self.lb    = -50.0*np.ones(self.n)
		self.ub    =  50.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'IKK1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]**2
		f2 = (x[0]-20.0)**2
		f3 = x[1]**2
		return np.array([f1, f2, f3])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_IM1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.ones(self.n)
		self.ub    = 2.0*np.ones(self.n)
		self.ub[1] = 4.0

		x0 = self.startp()

		self.name  = 'IM1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 2.0*np.sqrt(x[0])
		f2 = x[0]*(1.0-x[1])+5.0
		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_JIN1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'JIN1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.0
		f2 = 0.0
		for i in range(self.n): #i = 1,n
			f1 += x[0]**2
			f2 += (x[0]-2.0)**2

		f1 /= float(self.n)
		f2 /= float(self.n)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_JIN2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'JIN2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += x[i]

		gx = 1.0 + (9.0*gx)/float(self.n-1)

		f1 = x[0]
		f2 = gx*(1.0-np.sqrt(x[0]/gx))

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_JIN3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'JIN3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += x[i]

		gx = 1.0 + (9.0*gx)/float(self.n-1)

		f1 = x[0]
		f2 = gx*(1.0-(x[0]/gx)**2)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_JIN4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'JIN4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += x[i]

		gx = 1.0 + (9.0*gx)/float(self.n-1)

		f1 = x[0]
		f2 = gx*(1.0-(x[0]/gx)**0.25-(x[0]/gx)**4)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_KURSAWE():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 3
		self.n     = n
		self.q     = 2

		self.lb    = -5.0*np.ones(self.n)
		self.ub    =  5.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'KURSAWE'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.0
		for i in range(self.n-1): #i = 1,n-1
			f1 += (-10.0*np.exp(-0.2*np.sqrt(x[i]**2+x[i+1]**2)))

		f2 = 0.0
		for i in range(self.n): #i = 1,n
			f2 += (np.abs(x[i])**0.8+5.0*np.sin(x[i])**3.0)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L1ZDT4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 2

		self.lb    = -5.0*np.ones(self.n); self.lb[0] = 0.0
		self.ub    =  5.0*np.ones(self.n); self.ub[0] = 1.0
		self.A     = np.array([[1.  ,   0.        ,   0.         ,	   0.  ,   0.          ,  0.,	\
			     0.  ,   0.        ,     0.       ,      0.]	,\
			     [0.  ,   0.884043  , -0.272951    ,	-0.993822  ,  0.511197 ,   -0.0997948,	\
			-0.659756,     0.575496,      0.675617,      0.180332],	\
			     [0.  ,  -0.492722  ,  0.0646786   ,	-0.666503  , -0.945716 ,   -0.334582,	\
			0.611894 ,    0.281032 ,     0.508749 ,    -0.0265389],	\
			     [0.  ,   0.308861  , -0.0437502   ,	-0.374203  ,  0.207359 ,   -0.219433,	\
			0.914104 ,    0.184408 ,     0.520599 ,    -0.88565],	\
			     [0.  ,  -0.708948  , -0.37902     , 	 0.576578 ,   0.0194674,   -0.470262,	\
			0.572576 ,    0.351245 ,    -0.480477 ,     0.238261],	\
			     [0.  ,  -0.827302  ,  0.669248    ,	 0.494475  ,  0.691715 ,   -0.198585,	\
			0.0492812,    0.959669 ,     0.884086 ,    -0.218632],	\
			     [0.  ,  -0.715997  ,  0.220772    ,	 0.692356  ,  0.646453 ,   -0.401724,	\
			0.615443 ,   -0.0601957,    -0.748176 ,    -0.207987],	\
			     [0.  ,   0.613732  , -0.525712    ,	-0.995728  ,  0.389633 ,   -0.064173,	\
			0.662131 ,   -0.707048 ,    -0.340423 ,     0.60624],	\
			     [0.  ,  -0.160446  , -0.394585    ,	-0.167581 ,   0.0679849,    0.449799,	\
			0.733505 ,   -0.00918638,    0.00446808,    0.404396],	\
			     [0.  ,   0.162711   , 0.294454    ,	-0.563345  , -0.114993 ,    0.549589,	\
			-0.775141,     0.677726 ,     0.610715,    0.0850755]])

		x0 = self.startp()

		self.name  = 'L1ZDT4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		A = self.A
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += A[j,i]*x[j]

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += (y[i]**2-10.0*np.cos(4.0*np.pi*y[i]))

		gx = 1.0 + 10.0*float(self.n-1) + gx

		h = 1.0 - np.sqrt(f1/gx)
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L2ZDT1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507	,\
		0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617	,\
		0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749]	,\
		[-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599	,\
		-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477	,\
		0.238261	,	-0.1596		,	-0.827302	,	0.669248	,	0.494475	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086]	,\
		[-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176	,\
		-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423	,\
		0.60624		,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581	,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808]	,\
		[0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715	,\
		0.0850755	,	0.0419388	,	-0.323614	,	-0.973719	,	-0.680238	,\
		-0.270873	,	-0.209617	,	0.968436	,	0.908798	,	0.975851	,\
		-0.994918	,	-0.0621977	,	0.628171	,	0.761228	,	0.34372		,\
		-0.792042	,	-0.144765	,	-0.965748	,	0.0133606	,	-0.0260565]	,\
		[-0.742377	,	0.426391	,	0.408202	,	0.633885	,	-0.0351053	,\
		-0.723444	,	-0.577654	,	0.0276004	,	0.0712472	,	-0.622791	,\
		0.155451	,	0.442717	,	-0.792786	,	0.925785	,	0.670266	,\
		-0.865566	,	-0.638281	,	0.333094	,	0.477628	,	0.47261		,\
		0.23151		,	0.82132		,	-0.589803	,	0.796275	,	0.57713		,\
		0.101149	,	0.970191	,	0.532821	,	0.814769	,	-0.0687269]	,\
		[0.712758	,	-0.191812	,	-0.390938	,	0.952828	,	0.921519	,\
		0.923094	,	0.93011		,	-0.945394	,	-0.0934027	,	0.964123	,\
		-0.795609	,	-0.289563	,	0.614236	,	-0.670585	,	0.466877	,\
		0.144597	,	-0.206416	,	0.6937		,	-0.967958	,	-0.0951247	,\
		-0.942473	,	-0.610767	,	-0.655472	,	-0.0960986	,	-0.302779	,\
		-0.734976	,	-0.342188	,	-0.315861	,	-0.912834	,	0.24499]	,\
		[0.0969326	,	0.089775	,	-0.241157	,	0.0835558	,	-0.420236	,\
		-0.686633	,	-0.711276	,	-0.00325377	,	0.435196	,	-0.710002	,\
		0.00283691	,	-0.168757	,	-0.134045	,	-0.655235	,	0.172361	,\
		0.998291	,	0.376291	,	-0.962215	,	-0.363174	,	-0.88777	,\
		-0.519929	,	-0.560554	,	-0.984415	,	0.601529	,	-0.984103	,\
		-0.228237	,	-0.578066	,	0.307023	,	0.606123	,	0.959635]	,\
		[0.00225943	,	0.0101814	,	0.441456	,	0.0633629	,	0.406631	,\
		-0.0100638	,	-0.177972	,	-0.491075	,	0.537035	,	-0.924987	,\
		-0.699424	,	0.742285	,	0.0181443	,	0.718971	,	-0.0308272	,\
		0.086931	,	0.524476	,	0.956457	,	0.143024	,	0.616481	,\
		0.217909	,	-0.128427	,	-0.262427	,	-0.938208	,	-0.52479	,\
		0.12919		,	0.721925	,	0.766492	,	0.470845	,	-0.0976466]	,\
		[0.507807	,	0.804148	,	0.963269	,	0.357128	,	-0.832565	,\
		-0.312441	,	0.327779	,	0.184745	,	0.246139	,	-0.936814	,\
		-0.931734	,	-0.0327827	,	0.319293	,	0.044473	,	-0.641645	,\
		0.596118	,	-0.293934	,	-0.63373	,	0.409658	,	0.759892	,\
		-0.257078	,	0.939616	,	-0.227661	,	0.115754	,	0.10964		,\
		-0.240557	,	0.66842		,	0.855535	,	-0.451536	,	0.264961]	,\
		[-0.61366	,	-0.204783	,	-0.842476	,	-0.249524	,	-0.0985226	,\
		0.0671501	,	-0.527707	,	-0.509489	,	-0.883254	,	0.14851		,\
		-0.906465	,	0.496238	,	-0.853211	,	-0.779234	,	-0.979515	,\
		0.827175	,	0.228969	,	-0.402829	,	-0.970118	,	0.762559	,\
		0.506495	,	0.460303	,	0.897304	,	0.686003	,	0.739986	,\
		0.15731		,	0.281697	,	-0.922955	,	-0.780824	,	0.449716]	,\
		[0.125225	,	0.487649	,	0.147046	,	0.679639	,	0.593707	,\
		-0.311828	,	-0.797099	,	-0.35815	,	0.95808		,	0.907244	,\
		0.772426	,	0.720574	,	-0.873217	,	0.371431	,	-0.826029	,\
		0.942716	,	0.70609		,	-0.658158	,	-0.782185	,	-0.806743	,\
		-0.627986	,	-0.405551	,	-0.258495	,	-0.796524	,	0.222498	,\
		0.087545	,	-0.0917108	,	-0.62542	,	-0.110256	,	0.0417576]	,\
		[0.24476		,	0.941339	,	-0.613783	,	0.402772	,	0.300775,\
		-0.820314	,	-0.894233	,	-0.405896	,	0.0735439	,	0.486645	,\
		-0.394355	,	0.125097	,	-0.316386	,	-0.701215	,	-0.845742	,\
		0.2065		,	-0.413743	,	0.406725	,	-0.423813	,	-0.941255	,\
		-0.558804	,	0.312326	,	0.345314	,	0.319143	,	-0.644653	,\
		-0.0408415	,	0.176461	,	0.740113	,	0.470737	,	-0.914927]	,\
		[-0.591523	,	-0.606614	,	-0.181873	,	0.692975	,	0.50208		,\
		-0.536704	,	0.359652	,	0.839082	,	0.56817		,	-0.0776788	,\
		-0.00332785	,	0.459538	,	-0.518313	,	-0.270738	,	-0.629958	,\
		-0.755084	,	-0.721573	,	0.431107	,	-0.221877	,	0.32543		,\
		0.163743	,	0.0759916	,	0.695064	,	-0.656856	,	0.074163	,\
		0.264319	,	-0.73174	,	0.731548	,	-0.489341	,	0.678946]	,\
		[0.0271269	,	0.804879	,	-0.402973	,	0.800373	,	0.760082	,\
		-0.878364	,	0.176801	,	-0.548932	,	-0.225601	,	-0.164912	,\
		-0.208143	,	0.7768		,	-0.542743	,	-0.156021	,	0.671736	,\
		0.878648	,	-0.419588	,	-0.0752896	,	0.0299447	,	-0.494459	,\
		-0.72415	,	0.35978		,	-0.32646	,	-0.96605	,	0.0127605	,\
		0.563174	,	-0.814853	,	-0.949609	,	-0.526794	,	-0.801902]	,\
		[-0.753397	,	0.617418	,	0.689874	,	0.983384	,	0.668786	,\
		0.0304653	,	-0.625221	,	-0.13318	,	0.827343	,	-0.101358	,\
		-0.999522	,	-0.0525574	,	-0.458319	,	0.587409	,	-0.334639	,\
		0.0759643	,	0.0255827	,	0.128944	,	0.17317		,	-0.284309	,\
		0.287161	,	-0.550725	,	-0.433083	,	-0.242821	,	0.878879	,\
		0.691699	,	-0.660499	,	0.389985	,	0.599856	,	-0.711442]	,\
		[-0.798697	,	-0.244945	,	-0.942649	,	0.402856	,	-0.494672	,\
		0.439941	,	-0.88216	,	0.170196	,	0.650734	,	-0.0982391	,\
		-0.468732	,	0.342133	,	-0.838071	,	-0.832362	,	0.658177	,\
		-0.565361	,	0.149473	,	0.69331		,	-0.491848	,	0.74916		,\
		0.526025	,	-0.155339	,	0.0998096	,	0.468761	,	0.324649	,\
		0.128488	,	0.544144	,	-0.495222	,	0.965229	,	-0.79314]	,\
		[-0.545421	,	-0.500243	,	0.154371	,	0.170017	,	-0.259108	,\
		-0.868862	,	-0.50731	,	-0.848317	,	0.835712	,	0.616391	,\
		-0.442608	,	-0.158		,	0.313451	,	0.703748	,	-0.755984	,\
		-0.249443	,	0.491564	,	0.985068	,	0.678644	,	0.808324	,\
		0.81975		,	-0.435823	,	-0.839855	,	0.00282368	,	-0.569165	,\
		0.0884339	,	-0.222144	,	0.499412	,	-0.565198	,	0.64824]	,\
		[0.956914	,	-0.0620912	,	0.634479	,	0.928617	,	0.464664	,\
		0.377022	,	0.63047		,	-0.198619	,	-0.576153	,	0.565373	,\
		-0.524245	,	-0.187299	,	-0.614524	,	0.429316	,	-0.491171	,\
		0.399495	,	-0.333898	,	-0.646636	,	-0.0189709	,	-0.339605	,\
		-0.798791	,	0.0494081	,	0.367012	,	0.852545	,	0.43557		,\
		0.150039	,	-0.0454542	,	0.604861	,	-0.598288	,	-0.500696]	,\
		[0.249008	,	0.370711	,	-0.633174	,	-0.0121906	,	0.42006		,\
		0.169373	,	-0.975542	,	-0.0297852	,	0.80481		,	0.638317	,\
		-0.670967	,	0.935792	,	-0.35605	,	0.175773	,	0.878601	,\
		-0.275168	,	-0.932517	,	-0.372497	,	-0.0732907	,	-0.185493	,\
		-0.357004	,	0.314786	,	-0.229239	,	0.530256	,	-0.51327	,\
		0.44187		,	0.940309	,	-0.240334	,	-0.0276121	,	0.74383]	,\
		[-0.630329	,	-0.763376	,	0.62538		,	0.818945	,	0.891598	,\
		0.680494	,	0.471868	,	-0.769787	,	-0.878099	,	-0.973724	,\
		0.354362	,	-0.1792		,	-0.225034	,	-0.44548	,	0.598865	,\
		0.544005	,	-0.478962	,	0.327193	,	-0.525784	,	0.903179	,\
		-0.899248	,	0.156514	,	0.154329	,	0.499808	,	-0.836327	,\
		-0.802627	,	0.378082	,	-0.112673	,	-0.47926	,	-0.3355]	,\
		[-0.699445	,	0.237731	,	-0.324597	,	-0.800406	,	-0.42585	,\
		-0.710739	,	-0.144068	,	-0.828545	,	-0.800912	,	0.184654	,\
		-0.63675	,	-0.16696	,	0.240427	,	-0.513443	,	0.812664	,\
		0.744943	,	0.970612	,	0.00172899	,	-0.726378	,	-0.0985012	,\
		0.224232	,	0.16495		,	0.560077	,	-0.813112	,	0.112894	,\
		-0.0955366	,	0.0187107	,	0.913887	,	0.123076	,	0.550338]	,\
		[0.400334	,	-0.367816	,	0.198455	,	-0.983183	,	0.278976	,\
		0.714817	,	0.307911	,	0.812861	,	-0.403497	,	-0.784382	,\
		-0.161823	,	-0.120835	,	0.323172	,	0.583739	,	0.732924	,\
		-0.220603	,	-0.594121	,	0.935093	,	-0.216736	,	0.659318	,\
		-0.750417	,	-0.284773	,	-0.271496	,	0.491731	,	-0.712174	,\
		-0.763681	,	0.0781023	,	0.951666	,	0.734031	,	0.826912]	,\
		[0.57488		,	-0.361951	,	-0.0739728	,	0.91438		,	-0.391653,\
		0.0193887	,	0.412634	,	-0.169813	,	0.471794	,	0.660792	,\
		-0.350906	,	-0.612644	,	0.347876	,	0.112573	,	-0.501126	,\
		0.456761	,	-0.109004	,	0.289352	,	-0.566504	,	0.585042	,\
		0.584934	,	0.923676	,	0.895312	,	-0.161036	,	-0.995895	,\
		0.0853141	,	-0.583368	,	-0.157612	,	0.234119	,	0.875043]	,\
		[0.430805	,	0.706102	,	0.423887	,	0.296828	,	-0.265607	,\
		0.338806	,	-0.15829	,	0.642516	,	0.355126	,	0.174447	,\
		-0.975015	,	0.869905	,	-0.145285	,	-0.484002	,	-0.475966	,\
		-0.67704	,	0.996452	,	-0.0685748	,	-0.851985	,	0.416498	,\
		0.791047	,	-0.211323	,	-0.302819	,	0.640735	,	-0.317908	,\
		-0.116586	,	-0.896382	,	-0.817317	,	-0.948837	,	-0.597427]	,\
		[0.975863	,	-0.971371	,	-0.124115	,	0.4339		,	-0.254671	,\
		0.298068	,	-0.349803	,	-0.73185	,	0.488106	,	-0.0495073	,\
		0.253969	,	0.168116	,	0.148772	,	0.889593	,	-0.512213	,\
		-0.165437	,	0.666773	,	-0.976304	,	-0.170024	,	0.905794	,\
		0.473908	,	-0.855725	,	-0.0413591	,	-0.508661	,	0.443453	,\
		0.842925	,	-0.144503	,	0.936699	,	-0.443935	,	-0.182996]	,\
		[0.803564	,	0.960386	,	-0.0323329	,	0.638181	,	-0.895684	,\
		-0.360502	,	0.0646258	,	-0.202449	,	-0.717228	,	0.970489	,\
		0.404608	,	-0.0861868	,	-0.879417	,	-0.866462	,	-0.938336	,\
		-0.799685	,	0.213464	,	-0.932344	,	-0.668236	,	0.751366	,\
		-0.22712	,	-0.407783	,	0.657463	,	0.0970092	,	-0.579109	,\
		-0.868866	,	-0.504041	,	0.926483	,	0.169978	,	-0.00842563],\
		[-0.530324	,	0.282745	,	0.0255867	,	0.287686	,	0.410417	,\
		-0.766576	,	-0.536566	,	-0.628288	,	0.69665		,	0.820713	,\
		-0.506162	,	-0.404114	,	0.640099	,	-0.956123	,	-0.576586	,\
		0.435502	,	-0.470676	,	-0.367062	,	-0.831765	,	-0.294942	,\
		0.518991	,	0.922338	,	0.337886	,	-0.67474	,	-0.725667	,\
		0.916684	,	0.39175		,	0.759081	,	0.496979	,	-0.200691]	,\
		[0.0417966	,	-0.687391	,	0.438773	,	0.287357	,	0.316636	,\
		-0.262311	,	-0.0755541	,	-0.442313	,	0.621378	,	0.670105	,\
		0.060982	,	0.944162	,	0.643442	,	-0.750684	,	-0.639973	,\
		0.217424	,	0.592823	,	0.929094	,	-0.239135	,	-0.41628	,\
		0.570893	,	-0.0798988	,	-0.917135	,	-0.749545	,	-0.982047	,\
		0.0626998	,	-0.977963	,	0.660401	,	0.470569	,	-0.0528868]	,\
		[-0.00138645	,	0.931065	,	-0.748519	,	0.304188	,	-0.266153,\
		0.672524	,	-0.105179	,	-0.874749	,	-0.154355	,	-0.774656	,\
		-0.69654	,	0.433098	,	0.615897	,	-0.387919	,	-0.429779	,\
		0.650202	,	0.122306	,	-0.237727	,	0.626817	,	-0.227929	,\
		0.405916	,	0.483328	,	0.282047	,	-0.262206	,	0.784123	,\
		0.83125		,	-0.662272	,	0.702768	,	0.875814	,	-0.701221]	,\
		[0.553793	,	0.471795	,	0.769147	,	0.059668	,	-0.841617	,\
		-0.191179	,	-0.972471	,	-0.825361	,	0.779826	,	-0.917201	,\
		0.43272		,	0.10301		,	0.358771	,	0.793448	,	-0.0379954	,\
		-0.870112	,	0.600442	,	-0.990603	,	0.549151	,	0.512146	,\
		-0.795843	,	0.490091	,	0.372046	,	-0.549437	,	0.0964285	,\
		0.753047	,	-0.86284	,	-0.589688	,	0.178612	,	-0.720358]])

		x0 = self.startp()

		self.name  = 'L2ZDT1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += y[i]**2

		gx = 1.0 + 9.0*float(self.n-1) * gx

		h = 1.0 - np.sqrt(f1/gx)
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L2ZDT2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507	,\
		0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617	,\
		0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749]	,\
		[-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599	,\
		-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477	,\
		0.238261	,	-0.1596		,	-0.827302	,	0.669248	,	0.494475	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086]	,\
		[-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176	,\
		-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423	,\
		0.60624		,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581	,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808]	,\
		[0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715	,\
		0.0850755	,	0.0419388	,	-0.323614	,	-0.973719	,	-0.680238	,\
		-0.270873	,	-0.209617	,	0.968436	,	0.908798	,	0.975851	,\
		-0.994918	,	-0.0621977	,	0.628171	,	0.761228	,	0.34372		,\
		-0.792042	,	-0.144765	,	-0.965748	,	0.0133606	,	-0.0260565]	,\
		[-0.742377	,	0.426391	,	0.408202	,	0.633885	,	-0.0351053	,\
		-0.723444	,	-0.577654	,	0.0276004	,	0.0712472	,	-0.622791	,\
		0.155451	,	0.442717	,	-0.792786	,	0.925785	,	0.670266	,\
		-0.865566	,	-0.638281	,	0.333094	,	0.477628	,	0.47261		,\
		0.23151		,	0.82132		,	-0.589803	,	0.796275	,	0.57713		,\
		0.101149	,	0.970191	,	0.532821	,	0.814769	,	-0.0687269]	,\
		[0.712758	,	-0.191812	,	-0.390938	,	0.952828	,	0.921519	,\
		0.923094	,	0.93011		,	-0.945394	,	-0.0934027	,	0.964123	,\
		-0.795609	,	-0.289563	,	0.614236	,	-0.670585	,	0.466877	,\
		0.144597	,	-0.206416	,	0.6937		,	-0.967958	,	-0.0951247	,\
		-0.942473	,	-0.610767	,	-0.655472	,	-0.0960986	,	-0.302779	,\
		-0.734976	,	-0.342188	,	-0.315861	,	-0.912834	,	0.24499]	,\
		[0.0969326	,	0.089775	,	-0.241157	,	0.0835558	,	-0.420236	,\
		-0.686633	,	-0.711276	,	-0.00325377	,	0.435196	,	-0.710002	,\
		0.00283691	,	-0.168757	,	-0.134045	,	-0.655235	,	0.172361	,\
		0.998291	,	0.376291	,	-0.962215	,	-0.363174	,	-0.88777	,\
		-0.519929	,	-0.560554	,	-0.984415	,	0.601529	,	-0.984103	,\
		-0.228237	,	-0.578066	,	0.307023	,	0.606123	,	0.959635]	,\
		[0.00225943	,	0.0101814	,	0.441456	,	0.0633629	,	0.406631	,\
		-0.0100638	,	-0.177972	,	-0.491075	,	0.537035	,	-0.924987	,\
		-0.699424	,	0.742285	,	0.0181443	,	0.718971	,	-0.0308272	,\
		0.086931	,	0.524476	,	0.956457	,	0.143024	,	0.616481	,\
		0.217909	,	-0.128427	,	-0.262427	,	-0.938208	,	-0.52479	,\
		0.12919		,	0.721925	,	0.766492	,	0.470845	,	-0.0976466]	,\
		[0.507807	,	0.804148	,	0.963269	,	0.357128	,	-0.832565	,\
		-0.312441	,	0.327779	,	0.184745	,	0.246139	,	-0.936814	,\
		-0.931734	,	-0.0327827	,	0.319293	,	0.044473	,	-0.641645	,\
		0.596118	,	-0.293934	,	-0.63373	,	0.409658	,	0.759892	,\
		-0.257078	,	0.939616	,	-0.227661	,	0.115754	,	0.10964		,\
		-0.240557	,	0.66842		,	0.855535	,	-0.451536	,	0.264961]	,\
		[-0.61366	,	-0.204783	,	-0.842476	,	-0.249524	,	-0.0985226	,\
		0.0671501	,	-0.527707	,	-0.509489	,	-0.883254	,	0.14851		,\
		-0.906465	,	0.496238	,	-0.853211	,	-0.779234	,	-0.979515	,\
		0.827175	,	0.228969	,	-0.402829	,	-0.970118	,	0.762559	,\
		0.506495	,	0.460303	,	0.897304	,	0.686003	,	0.739986	,\
		0.15731		,	0.281697	,	-0.922955	,	-0.780824	,	0.449716]	,\
		[0.125225	,	0.487649	,	0.147046	,	0.679639	,	0.593707	,\
		-0.311828	,	-0.797099	,	-0.35815	,	0.95808		,	0.907244	,\
		0.772426	,	0.720574	,	-0.873217	,	0.371431	,	-0.826029	,\
		0.942716	,	0.70609		,	-0.658158	,	-0.782185	,	-0.806743	,\
		-0.627986	,	-0.405551	,	-0.258495	,	-0.796524	,	0.222498	,\
		0.087545	,	-0.0917108	,	-0.62542	,	-0.110256	,	0.0417576]	,\
		[0.24476		,	0.941339	,	-0.613783	,	0.402772	,	0.300775,\
		-0.820314	,	-0.894233	,	-0.405896	,	0.0735439	,	0.486645	,\
		-0.394355	,	0.125097	,	-0.316386	,	-0.701215	,	-0.845742	,\
		0.2065		,	-0.413743	,	0.406725	,	-0.423813	,	-0.941255	,\
		-0.558804	,	0.312326	,	0.345314	,	0.319143	,	-0.644653	,\
		-0.0408415	,	0.176461	,	0.740113	,	0.470737	,	-0.914927]	,\
		[-0.591523	,	-0.606614	,	-0.181873	,	0.692975	,	0.50208		,\
		-0.536704	,	0.359652	,	0.839082	,	0.56817		,	-0.0776788	,\
		-0.00332785	,	0.459538	,	-0.518313	,	-0.270738	,	-0.629958	,\
		-0.755084	,	-0.721573	,	0.431107	,	-0.221877	,	0.32543		,\
		0.163743	,	0.0759916	,	0.695064	,	-0.656856	,	0.074163	,\
		0.264319	,	-0.73174	,	0.731548	,	-0.489341	,	0.678946]	,\
		[0.0271269	,	0.804879	,	-0.402973	,	0.800373	,	0.760082	,\
		-0.878364	,	0.176801	,	-0.548932	,	-0.225601	,	-0.164912	,\
		-0.208143	,	0.7768		,	-0.542743	,	-0.156021	,	0.671736	,\
		0.878648	,	-0.419588	,	-0.0752896	,	0.0299447	,	-0.494459	,\
		-0.72415	,	0.35978		,	-0.32646	,	-0.96605	,	0.0127605	,\
		0.563174	,	-0.814853	,	-0.949609	,	-0.526794	,	-0.801902]	,\
		[-0.753397	,	0.617418	,	0.689874	,	0.983384	,	0.668786	,\
		0.0304653	,	-0.625221	,	-0.13318	,	0.827343	,	-0.101358	,\
		-0.999522	,	-0.0525574	,	-0.458319	,	0.587409	,	-0.334639	,\
		0.0759643	,	0.0255827	,	0.128944	,	0.17317		,	-0.284309	,\
		0.287161	,	-0.550725	,	-0.433083	,	-0.242821	,	0.878879	,\
		0.691699	,	-0.660499	,	0.389985	,	0.599856	,	-0.711442]	,\
		[-0.798697	,	-0.244945	,	-0.942649	,	0.402856	,	-0.494672	,\
		0.439941	,	-0.88216	,	0.170196	,	0.650734	,	-0.0982391	,\
		-0.468732	,	0.342133	,	-0.838071	,	-0.832362	,	0.658177	,\
		-0.565361	,	0.149473	,	0.69331		,	-0.491848	,	0.74916		,\
		0.526025	,	-0.155339	,	0.0998096	,	0.468761	,	0.324649	,\
		0.128488	,	0.544144	,	-0.495222	,	0.965229	,	-0.79314]	,\
		[-0.545421	,	-0.500243	,	0.154371	,	0.170017	,	-0.259108	,\
		-0.868862	,	-0.50731	,	-0.848317	,	0.835712	,	0.616391	,\
		-0.442608	,	-0.158		,	0.313451	,	0.703748	,	-0.755984	,\
		-0.249443	,	0.491564	,	0.985068	,	0.678644	,	0.808324	,\
		0.81975		,	-0.435823	,	-0.839855	,	0.00282368	,	-0.569165	,\
		0.0884339	,	-0.222144	,	0.499412	,	-0.565198	,	0.64824]	,\
		[0.956914	,	-0.0620912	,	0.634479	,	0.928617	,	0.464664	,\
		0.377022	,	0.63047		,	-0.198619	,	-0.576153	,	0.565373	,\
		-0.524245	,	-0.187299	,	-0.614524	,	0.429316	,	-0.491171	,\
		0.399495	,	-0.333898	,	-0.646636	,	-0.0189709	,	-0.339605	,\
		-0.798791	,	0.0494081	,	0.367012	,	0.852545	,	0.43557		,\
		0.150039	,	-0.0454542	,	0.604861	,	-0.598288	,	-0.500696]	,\
		[0.249008	,	0.370711	,	-0.633174	,	-0.0121906	,	0.42006		,\
		0.169373	,	-0.975542	,	-0.0297852	,	0.80481		,	0.638317	,\
		-0.670967	,	0.935792	,	-0.35605	,	0.175773	,	0.878601	,\
		-0.275168	,	-0.932517	,	-0.372497	,	-0.0732907	,	-0.185493	,\
		-0.357004	,	0.314786	,	-0.229239	,	0.530256	,	-0.51327	,\
		0.44187		,	0.940309	,	-0.240334	,	-0.0276121	,	0.74383]	,\
		[-0.630329	,	-0.763376	,	0.62538		,	0.818945	,	0.891598	,\
		0.680494	,	0.471868	,	-0.769787	,	-0.878099	,	-0.973724	,\
		0.354362	,	-0.1792		,	-0.225034	,	-0.44548	,	0.598865	,\
		0.544005	,	-0.478962	,	0.327193	,	-0.525784	,	0.903179	,\
		-0.899248	,	0.156514	,	0.154329	,	0.499808	,	-0.836327	,\
		-0.802627	,	0.378082	,	-0.112673	,	-0.47926	,	-0.3355]	,\
		[-0.699445	,	0.237731	,	-0.324597	,	-0.800406	,	-0.42585	,\
		-0.710739	,	-0.144068	,	-0.828545	,	-0.800912	,	0.184654	,\
		-0.63675	,	-0.16696	,	0.240427	,	-0.513443	,	0.812664	,\
		0.744943	,	0.970612	,	0.00172899	,	-0.726378	,	-0.0985012	,\
		0.224232	,	0.16495		,	0.560077	,	-0.813112	,	0.112894	,\
		-0.0955366	,	0.0187107	,	0.913887	,	0.123076	,	0.550338]	,\
		[0.400334	,	-0.367816	,	0.198455	,	-0.983183	,	0.278976	,\
		0.714817	,	0.307911	,	0.812861	,	-0.403497	,	-0.784382	,\
		-0.161823	,	-0.120835	,	0.323172	,	0.583739	,	0.732924	,\
		-0.220603	,	-0.594121	,	0.935093	,	-0.216736	,	0.659318	,\
		-0.750417	,	-0.284773	,	-0.271496	,	0.491731	,	-0.712174	,\
		-0.763681	,	0.0781023	,	0.951666	,	0.734031	,	0.826912]	,\
		[0.57488		,	-0.361951	,	-0.0739728	,	0.91438		,	-0.391653,\
		0.0193887	,	0.412634	,	-0.169813	,	0.471794	,	0.660792	,\
		-0.350906	,	-0.612644	,	0.347876	,	0.112573	,	-0.501126	,\
		0.456761	,	-0.109004	,	0.289352	,	-0.566504	,	0.585042	,\
		0.584934	,	0.923676	,	0.895312	,	-0.161036	,	-0.995895	,\
		0.0853141	,	-0.583368	,	-0.157612	,	0.234119	,	0.875043]	,\
		[0.430805	,	0.706102	,	0.423887	,	0.296828	,	-0.265607	,\
		0.338806	,	-0.15829	,	0.642516	,	0.355126	,	0.174447	,\
		-0.975015	,	0.869905	,	-0.145285	,	-0.484002	,	-0.475966	,\
		-0.67704	,	0.996452	,	-0.0685748	,	-0.851985	,	0.416498	,\
		0.791047	,	-0.211323	,	-0.302819	,	0.640735	,	-0.317908	,\
		-0.116586	,	-0.896382	,	-0.817317	,	-0.948837	,	-0.597427]	,\
		[0.975863	,	-0.971371	,	-0.124115	,	0.4339		,	-0.254671	,\
		0.298068	,	-0.349803	,	-0.73185	,	0.488106	,	-0.0495073	,\
		0.253969	,	0.168116	,	0.148772	,	0.889593	,	-0.512213	,\
		-0.165437	,	0.666773	,	-0.976304	,	-0.170024	,	0.905794	,\
		0.473908	,	-0.855725	,	-0.0413591	,	-0.508661	,	0.443453	,\
		0.842925	,	-0.144503	,	0.936699	,	-0.443935	,	-0.182996]	,\
		[0.803564	,	0.960386	,	-0.0323329	,	0.638181	,	-0.895684	,\
		-0.360502	,	0.0646258	,	-0.202449	,	-0.717228	,	0.970489	,\
		0.404608	,	-0.0861868	,	-0.879417	,	-0.866462	,	-0.938336	,\
		-0.799685	,	0.213464	,	-0.932344	,	-0.668236	,	0.751366	,\
		-0.22712	,	-0.407783	,	0.657463	,	0.0970092	,	-0.579109	,\
		-0.868866	,	-0.504041	,	0.926483	,	0.169978	,	-0.00842563],\
		[-0.530324	,	0.282745	,	0.0255867	,	0.287686	,	0.410417	,\
		-0.766576	,	-0.536566	,	-0.628288	,	0.69665		,	0.820713	,\
		-0.506162	,	-0.404114	,	0.640099	,	-0.956123	,	-0.576586	,\
		0.435502	,	-0.470676	,	-0.367062	,	-0.831765	,	-0.294942	,\
		0.518991	,	0.922338	,	0.337886	,	-0.67474	,	-0.725667	,\
		0.916684	,	0.39175		,	0.759081	,	0.496979	,	-0.200691]	,\
		[0.0417966	,	-0.687391	,	0.438773	,	0.287357	,	0.316636	,\
		-0.262311	,	-0.0755541	,	-0.442313	,	0.621378	,	0.670105	,\
		0.060982	,	0.944162	,	0.643442	,	-0.750684	,	-0.639973	,\
		0.217424	,	0.592823	,	0.929094	,	-0.239135	,	-0.41628	,\
		0.570893	,	-0.0798988	,	-0.917135	,	-0.749545	,	-0.982047	,\
		0.0626998	,	-0.977963	,	0.660401	,	0.470569	,	-0.0528868]	,\
		[-0.00138645	,	0.931065	,	-0.748519	,	0.304188	,	-0.266153,\
		0.672524	,	-0.105179	,	-0.874749	,	-0.154355	,	-0.774656	,\
		-0.69654	,	0.433098	,	0.615897	,	-0.387919	,	-0.429779	,\
		0.650202	,	0.122306	,	-0.237727	,	0.626817	,	-0.227929	,\
		0.405916	,	0.483328	,	0.282047	,	-0.262206	,	0.784123	,\
		0.83125		,	-0.662272	,	0.702768	,	0.875814	,	-0.701221]	,\
		[0.553793	,	0.471795	,	0.769147	,	0.059668	,	-0.841617	,\
		-0.191179	,	-0.972471	,	-0.825361	,	0.779826	,	-0.917201	,\
		0.43272		,	0.10301		,	0.358771	,	0.793448	,	-0.0379954	,\
		-0.870112	,	0.600442	,	-0.990603	,	0.549151	,	0.512146	,\
		-0.795843	,	0.490091	,	0.372046	,	-0.549437	,	0.0964285	,\
		0.753047	,	-0.86284	,	-0.589688	,	0.178612	,	-0.720358]])

		x0 = self.startp()

		self.name  = 'L2ZDT2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += y[i]**2

		gx = 1.0 + 9.0*float(self.n-1) * gx

		h = 1.0-(f1/gx)**2
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L2ZDT3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507	,\
		0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617	,\
		0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749]	,\
		[-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599	,\
		-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477	,\
		0.238261	,	-0.1596		,	-0.827302	,	0.669248	,	0.494475	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086]	,\
		[-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176	,\
		-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423	,\
		0.60624		,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581	,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808]	,\
		[0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715	,\
		0.0850755	,	0.0419388	,	-0.323614	,	-0.973719	,	-0.680238	,\
		-0.270873	,	-0.209617	,	0.968436	,	0.908798	,	0.975851	,\
		-0.994918	,	-0.0621977	,	0.628171	,	0.761228	,	0.34372		,\
		-0.792042	,	-0.144765	,	-0.965748	,	0.0133606	,	-0.0260565]	,\
		[-0.742377	,	0.426391	,	0.408202	,	0.633885	,	-0.0351053	,\
		-0.723444	,	-0.577654	,	0.0276004	,	0.0712472	,	-0.622791	,\
		0.155451	,	0.442717	,	-0.792786	,	0.925785	,	0.670266	,\
		-0.865566	,	-0.638281	,	0.333094	,	0.477628	,	0.47261		,\
		0.23151		,	0.82132		,	-0.589803	,	0.796275	,	0.57713		,\
		0.101149	,	0.970191	,	0.532821	,	0.814769	,	-0.0687269]	,\
		[0.712758	,	-0.191812	,	-0.390938	,	0.952828	,	0.921519	,\
		0.923094	,	0.93011		,	-0.945394	,	-0.0934027	,	0.964123	,\
		-0.795609	,	-0.289563	,	0.614236	,	-0.670585	,	0.466877	,\
		0.144597	,	-0.206416	,	0.6937		,	-0.967958	,	-0.0951247	,\
		-0.942473	,	-0.610767	,	-0.655472	,	-0.0960986	,	-0.302779	,\
		-0.734976	,	-0.342188	,	-0.315861	,	-0.912834	,	0.24499]	,\
		[0.0969326	,	0.089775	,	-0.241157	,	0.0835558	,	-0.420236	,\
		-0.686633	,	-0.711276	,	-0.00325377	,	0.435196	,	-0.710002	,\
		0.00283691	,	-0.168757	,	-0.134045	,	-0.655235	,	0.172361	,\
		0.998291	,	0.376291	,	-0.962215	,	-0.363174	,	-0.88777	,\
		-0.519929	,	-0.560554	,	-0.984415	,	0.601529	,	-0.984103	,\
		-0.228237	,	-0.578066	,	0.307023	,	0.606123	,	0.959635]	,\
		[0.00225943	,	0.0101814	,	0.441456	,	0.0633629	,	0.406631	,\
		-0.0100638	,	-0.177972	,	-0.491075	,	0.537035	,	-0.924987	,\
		-0.699424	,	0.742285	,	0.0181443	,	0.718971	,	-0.0308272	,\
		0.086931	,	0.524476	,	0.956457	,	0.143024	,	0.616481	,\
		0.217909	,	-0.128427	,	-0.262427	,	-0.938208	,	-0.52479	,\
		0.12919		,	0.721925	,	0.766492	,	0.470845	,	-0.0976466]	,\
		[0.507807	,	0.804148	,	0.963269	,	0.357128	,	-0.832565	,\
		-0.312441	,	0.327779	,	0.184745	,	0.246139	,	-0.936814	,\
		-0.931734	,	-0.0327827	,	0.319293	,	0.044473	,	-0.641645	,\
		0.596118	,	-0.293934	,	-0.63373	,	0.409658	,	0.759892	,\
		-0.257078	,	0.939616	,	-0.227661	,	0.115754	,	0.10964		,\
		-0.240557	,	0.66842		,	0.855535	,	-0.451536	,	0.264961]	,\
		[-0.61366	,	-0.204783	,	-0.842476	,	-0.249524	,	-0.0985226	,\
		0.0671501	,	-0.527707	,	-0.509489	,	-0.883254	,	0.14851		,\
		-0.906465	,	0.496238	,	-0.853211	,	-0.779234	,	-0.979515	,\
		0.827175	,	0.228969	,	-0.402829	,	-0.970118	,	0.762559	,\
		0.506495	,	0.460303	,	0.897304	,	0.686003	,	0.739986	,\
		0.15731		,	0.281697	,	-0.922955	,	-0.780824	,	0.449716]	,\
		[0.125225	,	0.487649	,	0.147046	,	0.679639	,	0.593707	,\
		-0.311828	,	-0.797099	,	-0.35815	,	0.95808		,	0.907244	,\
		0.772426	,	0.720574	,	-0.873217	,	0.371431	,	-0.826029	,\
		0.942716	,	0.70609		,	-0.658158	,	-0.782185	,	-0.806743	,\
		-0.627986	,	-0.405551	,	-0.258495	,	-0.796524	,	0.222498	,\
		0.087545	,	-0.0917108	,	-0.62542	,	-0.110256	,	0.0417576]	,\
		[0.24476		,	0.941339	,	-0.613783	,	0.402772	,	0.300775,\
		-0.820314	,	-0.894233	,	-0.405896	,	0.0735439	,	0.486645	,\
		-0.394355	,	0.125097	,	-0.316386	,	-0.701215	,	-0.845742	,\
		0.2065		,	-0.413743	,	0.406725	,	-0.423813	,	-0.941255	,\
		-0.558804	,	0.312326	,	0.345314	,	0.319143	,	-0.644653	,\
		-0.0408415	,	0.176461	,	0.740113	,	0.470737	,	-0.914927]	,\
		[-0.591523	,	-0.606614	,	-0.181873	,	0.692975	,	0.50208		,\
		-0.536704	,	0.359652	,	0.839082	,	0.56817		,	-0.0776788	,\
		-0.00332785	,	0.459538	,	-0.518313	,	-0.270738	,	-0.629958	,\
		-0.755084	,	-0.721573	,	0.431107	,	-0.221877	,	0.32543		,\
		0.163743	,	0.0759916	,	0.695064	,	-0.656856	,	0.074163	,\
		0.264319	,	-0.73174	,	0.731548	,	-0.489341	,	0.678946]	,\
		[0.0271269	,	0.804879	,	-0.402973	,	0.800373	,	0.760082	,\
		-0.878364	,	0.176801	,	-0.548932	,	-0.225601	,	-0.164912	,\
		-0.208143	,	0.7768		,	-0.542743	,	-0.156021	,	0.671736	,\
		0.878648	,	-0.419588	,	-0.0752896	,	0.0299447	,	-0.494459	,\
		-0.72415	,	0.35978		,	-0.32646	,	-0.96605	,	0.0127605	,\
		0.563174	,	-0.814853	,	-0.949609	,	-0.526794	,	-0.801902]	,\
		[-0.753397	,	0.617418	,	0.689874	,	0.983384	,	0.668786	,\
		0.0304653	,	-0.625221	,	-0.13318	,	0.827343	,	-0.101358	,\
		-0.999522	,	-0.0525574	,	-0.458319	,	0.587409	,	-0.334639	,\
		0.0759643	,	0.0255827	,	0.128944	,	0.17317		,	-0.284309	,\
		0.287161	,	-0.550725	,	-0.433083	,	-0.242821	,	0.878879	,\
		0.691699	,	-0.660499	,	0.389985	,	0.599856	,	-0.711442]	,\
		[-0.798697	,	-0.244945	,	-0.942649	,	0.402856	,	-0.494672	,\
		0.439941	,	-0.88216	,	0.170196	,	0.650734	,	-0.0982391	,\
		-0.468732	,	0.342133	,	-0.838071	,	-0.832362	,	0.658177	,\
		-0.565361	,	0.149473	,	0.69331		,	-0.491848	,	0.74916		,\
		0.526025	,	-0.155339	,	0.0998096	,	0.468761	,	0.324649	,\
		0.128488	,	0.544144	,	-0.495222	,	0.965229	,	-0.79314]	,\
		[-0.545421	,	-0.500243	,	0.154371	,	0.170017	,	-0.259108	,\
		-0.868862	,	-0.50731	,	-0.848317	,	0.835712	,	0.616391	,\
		-0.442608	,	-0.158		,	0.313451	,	0.703748	,	-0.755984	,\
		-0.249443	,	0.491564	,	0.985068	,	0.678644	,	0.808324	,\
		0.81975		,	-0.435823	,	-0.839855	,	0.00282368	,	-0.569165	,\
		0.0884339	,	-0.222144	,	0.499412	,	-0.565198	,	0.64824]	,\
		[0.956914	,	-0.0620912	,	0.634479	,	0.928617	,	0.464664	,\
		0.377022	,	0.63047		,	-0.198619	,	-0.576153	,	0.565373	,\
		-0.524245	,	-0.187299	,	-0.614524	,	0.429316	,	-0.491171	,\
		0.399495	,	-0.333898	,	-0.646636	,	-0.0189709	,	-0.339605	,\
		-0.798791	,	0.0494081	,	0.367012	,	0.852545	,	0.43557		,\
		0.150039	,	-0.0454542	,	0.604861	,	-0.598288	,	-0.500696]	,\
		[0.249008	,	0.370711	,	-0.633174	,	-0.0121906	,	0.42006		,\
		0.169373	,	-0.975542	,	-0.0297852	,	0.80481		,	0.638317	,\
		-0.670967	,	0.935792	,	-0.35605	,	0.175773	,	0.878601	,\
		-0.275168	,	-0.932517	,	-0.372497	,	-0.0732907	,	-0.185493	,\
		-0.357004	,	0.314786	,	-0.229239	,	0.530256	,	-0.51327	,\
		0.44187		,	0.940309	,	-0.240334	,	-0.0276121	,	0.74383]	,\
		[-0.630329	,	-0.763376	,	0.62538		,	0.818945	,	0.891598	,\
		0.680494	,	0.471868	,	-0.769787	,	-0.878099	,	-0.973724	,\
		0.354362	,	-0.1792		,	-0.225034	,	-0.44548	,	0.598865	,\
		0.544005	,	-0.478962	,	0.327193	,	-0.525784	,	0.903179	,\
		-0.899248	,	0.156514	,	0.154329	,	0.499808	,	-0.836327	,\
		-0.802627	,	0.378082	,	-0.112673	,	-0.47926	,	-0.3355]	,\
		[-0.699445	,	0.237731	,	-0.324597	,	-0.800406	,	-0.42585	,\
		-0.710739	,	-0.144068	,	-0.828545	,	-0.800912	,	0.184654	,\
		-0.63675	,	-0.16696	,	0.240427	,	-0.513443	,	0.812664	,\
		0.744943	,	0.970612	,	0.00172899	,	-0.726378	,	-0.0985012	,\
		0.224232	,	0.16495		,	0.560077	,	-0.813112	,	0.112894	,\
		-0.0955366	,	0.0187107	,	0.913887	,	0.123076	,	0.550338]	,\
		[0.400334	,	-0.367816	,	0.198455	,	-0.983183	,	0.278976	,\
		0.714817	,	0.307911	,	0.812861	,	-0.403497	,	-0.784382	,\
		-0.161823	,	-0.120835	,	0.323172	,	0.583739	,	0.732924	,\
		-0.220603	,	-0.594121	,	0.935093	,	-0.216736	,	0.659318	,\
		-0.750417	,	-0.284773	,	-0.271496	,	0.491731	,	-0.712174	,\
		-0.763681	,	0.0781023	,	0.951666	,	0.734031	,	0.826912]	,\
		[0.57488		,	-0.361951	,	-0.0739728	,	0.91438		,	-0.391653,\
		0.0193887	,	0.412634	,	-0.169813	,	0.471794	,	0.660792	,\
		-0.350906	,	-0.612644	,	0.347876	,	0.112573	,	-0.501126	,\
		0.456761	,	-0.109004	,	0.289352	,	-0.566504	,	0.585042	,\
		0.584934	,	0.923676	,	0.895312	,	-0.161036	,	-0.995895	,\
		0.0853141	,	-0.583368	,	-0.157612	,	0.234119	,	0.875043]	,\
		[0.430805	,	0.706102	,	0.423887	,	0.296828	,	-0.265607	,\
		0.338806	,	-0.15829	,	0.642516	,	0.355126	,	0.174447	,\
		-0.975015	,	0.869905	,	-0.145285	,	-0.484002	,	-0.475966	,\
		-0.67704	,	0.996452	,	-0.0685748	,	-0.851985	,	0.416498	,\
		0.791047	,	-0.211323	,	-0.302819	,	0.640735	,	-0.317908	,\
		-0.116586	,	-0.896382	,	-0.817317	,	-0.948837	,	-0.597427]	,\
		[0.975863	,	-0.971371	,	-0.124115	,	0.4339		,	-0.254671	,\
		0.298068	,	-0.349803	,	-0.73185	,	0.488106	,	-0.0495073	,\
		0.253969	,	0.168116	,	0.148772	,	0.889593	,	-0.512213	,\
		-0.165437	,	0.666773	,	-0.976304	,	-0.170024	,	0.905794	,\
		0.473908	,	-0.855725	,	-0.0413591	,	-0.508661	,	0.443453	,\
		0.842925	,	-0.144503	,	0.936699	,	-0.443935	,	-0.182996]	,\
		[0.803564	,	0.960386	,	-0.0323329	,	0.638181	,	-0.895684	,\
		-0.360502	,	0.0646258	,	-0.202449	,	-0.717228	,	0.970489	,\
		0.404608	,	-0.0861868	,	-0.879417	,	-0.866462	,	-0.938336	,\
		-0.799685	,	0.213464	,	-0.932344	,	-0.668236	,	0.751366	,\
		-0.22712	,	-0.407783	,	0.657463	,	0.0970092	,	-0.579109	,\
		-0.868866	,	-0.504041	,	0.926483	,	0.169978	,	-0.00842563],\
		[-0.530324	,	0.282745	,	0.0255867	,	0.287686	,	0.410417	,\
		-0.766576	,	-0.536566	,	-0.628288	,	0.69665		,	0.820713	,\
		-0.506162	,	-0.404114	,	0.640099	,	-0.956123	,	-0.576586	,\
		0.435502	,	-0.470676	,	-0.367062	,	-0.831765	,	-0.294942	,\
		0.518991	,	0.922338	,	0.337886	,	-0.67474	,	-0.725667	,\
		0.916684	,	0.39175		,	0.759081	,	0.496979	,	-0.200691]	,\
		[0.0417966	,	-0.687391	,	0.438773	,	0.287357	,	0.316636	,\
		-0.262311	,	-0.0755541	,	-0.442313	,	0.621378	,	0.670105	,\
		0.060982	,	0.944162	,	0.643442	,	-0.750684	,	-0.639973	,\
		0.217424	,	0.592823	,	0.929094	,	-0.239135	,	-0.41628	,\
		0.570893	,	-0.0798988	,	-0.917135	,	-0.749545	,	-0.982047	,\
		0.0626998	,	-0.977963	,	0.660401	,	0.470569	,	-0.0528868]	,\
		[-0.00138645	,	0.931065	,	-0.748519	,	0.304188	,	-0.266153,\
		0.672524	,	-0.105179	,	-0.874749	,	-0.154355	,	-0.774656	,\
		-0.69654	,	0.433098	,	0.615897	,	-0.387919	,	-0.429779	,\
		0.650202	,	0.122306	,	-0.237727	,	0.626817	,	-0.227929	,\
		0.405916	,	0.483328	,	0.282047	,	-0.262206	,	0.784123	,\
		0.83125		,	-0.662272	,	0.702768	,	0.875814	,	-0.701221]	,\
		[0.553793	,	0.471795	,	0.769147	,	0.059668	,	-0.841617	,\
		-0.191179	,	-0.972471	,	-0.825361	,	0.779826	,	-0.917201	,\
		0.43272		,	0.10301		,	0.358771	,	0.793448	,	-0.0379954	,\
		-0.870112	,	0.600442	,	-0.990603	,	0.549151	,	0.512146	,\
		-0.795843	,	0.490091	,	0.372046	,	-0.549437	,	0.0964285	,\
		0.753047	,	-0.86284	,	-0.589688	,	0.178612	,	-0.720358]])

		x0 = self.startp()

		self.name  = 'L2ZDT3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += y[i]**2

		gx = 1.0 + 9.0*float(self.n-1) * gx

		h = 1.0-np.sqrt(f1/gx)-(f1/gx)*np.sin(10.0*np.pi*f1)
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L2ZDT4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507	,\
		0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617	,\
		0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749]	,\
		[-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599	,\
		-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477	,\
		0.238261	,	-0.1596		,	-0.827302	,	0.669248	,	0.494475	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086]	,\
		[-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176	,\
		-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423	,\
		0.60624		,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581	,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808]	,\
		[0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715	,\
		0.0850755	,	0.0419388	,	-0.323614	,	-0.973719	,	-0.680238	,\
		-0.270873	,	-0.209617	,	0.968436	,	0.908798	,	0.975851	,\
		-0.994918	,	-0.0621977	,	0.628171	,	0.761228	,	0.34372		,\
		-0.792042	,	-0.144765	,	-0.965748	,	0.0133606	,	-0.0260565]	,\
		[-0.742377	,	0.426391	,	0.408202	,	0.633885	,	-0.0351053	,\
		-0.723444	,	-0.577654	,	0.0276004	,	0.0712472	,	-0.622791	,\
		0.155451	,	0.442717	,	-0.792786	,	0.925785	,	0.670266	,\
		-0.865566	,	-0.638281	,	0.333094	,	0.477628	,	0.47261		,\
		0.23151		,	0.82132		,	-0.589803	,	0.796275	,	0.57713		,\
		0.101149	,	0.970191	,	0.532821	,	0.814769	,	-0.0687269]	,\
		[0.712758	,	-0.191812	,	-0.390938	,	0.952828	,	0.921519	,\
		0.923094	,	0.93011		,	-0.945394	,	-0.0934027	,	0.964123	,\
		-0.795609	,	-0.289563	,	0.614236	,	-0.670585	,	0.466877	,\
		0.144597	,	-0.206416	,	0.6937		,	-0.967958	,	-0.0951247	,\
		-0.942473	,	-0.610767	,	-0.655472	,	-0.0960986	,	-0.302779	,\
		-0.734976	,	-0.342188	,	-0.315861	,	-0.912834	,	0.24499]	,\
		[0.0969326	,	0.089775	,	-0.241157	,	0.0835558	,	-0.420236	,\
		-0.686633	,	-0.711276	,	-0.00325377	,	0.435196	,	-0.710002	,\
		0.00283691	,	-0.168757	,	-0.134045	,	-0.655235	,	0.172361	,\
		0.998291	,	0.376291	,	-0.962215	,	-0.363174	,	-0.88777	,\
		-0.519929	,	-0.560554	,	-0.984415	,	0.601529	,	-0.984103	,\
		-0.228237	,	-0.578066	,	0.307023	,	0.606123	,	0.959635]	,\
		[0.00225943	,	0.0101814	,	0.441456	,	0.0633629	,	0.406631	,\
		-0.0100638	,	-0.177972	,	-0.491075	,	0.537035	,	-0.924987	,\
		-0.699424	,	0.742285	,	0.0181443	,	0.718971	,	-0.0308272	,\
		0.086931	,	0.524476	,	0.956457	,	0.143024	,	0.616481	,\
		0.217909	,	-0.128427	,	-0.262427	,	-0.938208	,	-0.52479	,\
		0.12919		,	0.721925	,	0.766492	,	0.470845	,	-0.0976466]	,\
		[0.507807	,	0.804148	,	0.963269	,	0.357128	,	-0.832565	,\
		-0.312441	,	0.327779	,	0.184745	,	0.246139	,	-0.936814	,\
		-0.931734	,	-0.0327827	,	0.319293	,	0.044473	,	-0.641645	,\
		0.596118	,	-0.293934	,	-0.63373	,	0.409658	,	0.759892	,\
		-0.257078	,	0.939616	,	-0.227661	,	0.115754	,	0.10964		,\
		-0.240557	,	0.66842		,	0.855535	,	-0.451536	,	0.264961]	,\
		[-0.61366	,	-0.204783	,	-0.842476	,	-0.249524	,	-0.0985226	,\
		0.0671501	,	-0.527707	,	-0.509489	,	-0.883254	,	0.14851		,\
		-0.906465	,	0.496238	,	-0.853211	,	-0.779234	,	-0.979515	,\
		0.827175	,	0.228969	,	-0.402829	,	-0.970118	,	0.762559	,\
		0.506495	,	0.460303	,	0.897304	,	0.686003	,	0.739986	,\
		0.15731		,	0.281697	,	-0.922955	,	-0.780824	,	0.449716]	,\
		[0.125225	,	0.487649	,	0.147046	,	0.679639	,	0.593707	,\
		-0.311828	,	-0.797099	,	-0.35815	,	0.95808		,	0.907244	,\
		0.772426	,	0.720574	,	-0.873217	,	0.371431	,	-0.826029	,\
		0.942716	,	0.70609		,	-0.658158	,	-0.782185	,	-0.806743	,\
		-0.627986	,	-0.405551	,	-0.258495	,	-0.796524	,	0.222498	,\
		0.087545	,	-0.0917108	,	-0.62542	,	-0.110256	,	0.0417576]	,\
		[0.24476		,	0.941339	,	-0.613783	,	0.402772	,	0.300775,\
		-0.820314	,	-0.894233	,	-0.405896	,	0.0735439	,	0.486645	,\
		-0.394355	,	0.125097	,	-0.316386	,	-0.701215	,	-0.845742	,\
		0.2065		,	-0.413743	,	0.406725	,	-0.423813	,	-0.941255	,\
		-0.558804	,	0.312326	,	0.345314	,	0.319143	,	-0.644653	,\
		-0.0408415	,	0.176461	,	0.740113	,	0.470737	,	-0.914927]	,\
		[-0.591523	,	-0.606614	,	-0.181873	,	0.692975	,	0.50208		,\
		-0.536704	,	0.359652	,	0.839082	,	0.56817		,	-0.0776788	,\
		-0.00332785	,	0.459538	,	-0.518313	,	-0.270738	,	-0.629958	,\
		-0.755084	,	-0.721573	,	0.431107	,	-0.221877	,	0.32543		,\
		0.163743	,	0.0759916	,	0.695064	,	-0.656856	,	0.074163	,\
		0.264319	,	-0.73174	,	0.731548	,	-0.489341	,	0.678946]	,\
		[0.0271269	,	0.804879	,	-0.402973	,	0.800373	,	0.760082	,\
		-0.878364	,	0.176801	,	-0.548932	,	-0.225601	,	-0.164912	,\
		-0.208143	,	0.7768		,	-0.542743	,	-0.156021	,	0.671736	,\
		0.878648	,	-0.419588	,	-0.0752896	,	0.0299447	,	-0.494459	,\
		-0.72415	,	0.35978		,	-0.32646	,	-0.96605	,	0.0127605	,\
		0.563174	,	-0.814853	,	-0.949609	,	-0.526794	,	-0.801902]	,\
		[-0.753397	,	0.617418	,	0.689874	,	0.983384	,	0.668786	,\
		0.0304653	,	-0.625221	,	-0.13318	,	0.827343	,	-0.101358	,\
		-0.999522	,	-0.0525574	,	-0.458319	,	0.587409	,	-0.334639	,\
		0.0759643	,	0.0255827	,	0.128944	,	0.17317		,	-0.284309	,\
		0.287161	,	-0.550725	,	-0.433083	,	-0.242821	,	0.878879	,\
		0.691699	,	-0.660499	,	0.389985	,	0.599856	,	-0.711442]	,\
		[-0.798697	,	-0.244945	,	-0.942649	,	0.402856	,	-0.494672	,\
		0.439941	,	-0.88216	,	0.170196	,	0.650734	,	-0.0982391	,\
		-0.468732	,	0.342133	,	-0.838071	,	-0.832362	,	0.658177	,\
		-0.565361	,	0.149473	,	0.69331		,	-0.491848	,	0.74916		,\
		0.526025	,	-0.155339	,	0.0998096	,	0.468761	,	0.324649	,\
		0.128488	,	0.544144	,	-0.495222	,	0.965229	,	-0.79314]	,\
		[-0.545421	,	-0.500243	,	0.154371	,	0.170017	,	-0.259108	,\
		-0.868862	,	-0.50731	,	-0.848317	,	0.835712	,	0.616391	,\
		-0.442608	,	-0.158		,	0.313451	,	0.703748	,	-0.755984	,\
		-0.249443	,	0.491564	,	0.985068	,	0.678644	,	0.808324	,\
		0.81975		,	-0.435823	,	-0.839855	,	0.00282368	,	-0.569165	,\
		0.0884339	,	-0.222144	,	0.499412	,	-0.565198	,	0.64824]	,\
		[0.956914	,	-0.0620912	,	0.634479	,	0.928617	,	0.464664	,\
		0.377022	,	0.63047		,	-0.198619	,	-0.576153	,	0.565373	,\
		-0.524245	,	-0.187299	,	-0.614524	,	0.429316	,	-0.491171	,\
		0.399495	,	-0.333898	,	-0.646636	,	-0.0189709	,	-0.339605	,\
		-0.798791	,	0.0494081	,	0.367012	,	0.852545	,	0.43557		,\
		0.150039	,	-0.0454542	,	0.604861	,	-0.598288	,	-0.500696]	,\
		[0.249008	,	0.370711	,	-0.633174	,	-0.0121906	,	0.42006		,\
		0.169373	,	-0.975542	,	-0.0297852	,	0.80481		,	0.638317	,\
		-0.670967	,	0.935792	,	-0.35605	,	0.175773	,	0.878601	,\
		-0.275168	,	-0.932517	,	-0.372497	,	-0.0732907	,	-0.185493	,\
		-0.357004	,	0.314786	,	-0.229239	,	0.530256	,	-0.51327	,\
		0.44187		,	0.940309	,	-0.240334	,	-0.0276121	,	0.74383]	,\
		[-0.630329	,	-0.763376	,	0.62538		,	0.818945	,	0.891598	,\
		0.680494	,	0.471868	,	-0.769787	,	-0.878099	,	-0.973724	,\
		0.354362	,	-0.1792		,	-0.225034	,	-0.44548	,	0.598865	,\
		0.544005	,	-0.478962	,	0.327193	,	-0.525784	,	0.903179	,\
		-0.899248	,	0.156514	,	0.154329	,	0.499808	,	-0.836327	,\
		-0.802627	,	0.378082	,	-0.112673	,	-0.47926	,	-0.3355]	,\
		[-0.699445	,	0.237731	,	-0.324597	,	-0.800406	,	-0.42585	,\
		-0.710739	,	-0.144068	,	-0.828545	,	-0.800912	,	0.184654	,\
		-0.63675	,	-0.16696	,	0.240427	,	-0.513443	,	0.812664	,\
		0.744943	,	0.970612	,	0.00172899	,	-0.726378	,	-0.0985012	,\
		0.224232	,	0.16495		,	0.560077	,	-0.813112	,	0.112894	,\
		-0.0955366	,	0.0187107	,	0.913887	,	0.123076	,	0.550338]	,\
		[0.400334	,	-0.367816	,	0.198455	,	-0.983183	,	0.278976	,\
		0.714817	,	0.307911	,	0.812861	,	-0.403497	,	-0.784382	,\
		-0.161823	,	-0.120835	,	0.323172	,	0.583739	,	0.732924	,\
		-0.220603	,	-0.594121	,	0.935093	,	-0.216736	,	0.659318	,\
		-0.750417	,	-0.284773	,	-0.271496	,	0.491731	,	-0.712174	,\
		-0.763681	,	0.0781023	,	0.951666	,	0.734031	,	0.826912]	,\
		[0.57488		,	-0.361951	,	-0.0739728	,	0.91438		,	-0.391653,\
		0.0193887	,	0.412634	,	-0.169813	,	0.471794	,	0.660792	,\
		-0.350906	,	-0.612644	,	0.347876	,	0.112573	,	-0.501126	,\
		0.456761	,	-0.109004	,	0.289352	,	-0.566504	,	0.585042	,\
		0.584934	,	0.923676	,	0.895312	,	-0.161036	,	-0.995895	,\
		0.0853141	,	-0.583368	,	-0.157612	,	0.234119	,	0.875043]	,\
		[0.430805	,	0.706102	,	0.423887	,	0.296828	,	-0.265607	,\
		0.338806	,	-0.15829	,	0.642516	,	0.355126	,	0.174447	,\
		-0.975015	,	0.869905	,	-0.145285	,	-0.484002	,	-0.475966	,\
		-0.67704	,	0.996452	,	-0.0685748	,	-0.851985	,	0.416498	,\
		0.791047	,	-0.211323	,	-0.302819	,	0.640735	,	-0.317908	,\
		-0.116586	,	-0.896382	,	-0.817317	,	-0.948837	,	-0.597427]	,\
		[0.975863	,	-0.971371	,	-0.124115	,	0.4339		,	-0.254671	,\
		0.298068	,	-0.349803	,	-0.73185	,	0.488106	,	-0.0495073	,\
		0.253969	,	0.168116	,	0.148772	,	0.889593	,	-0.512213	,\
		-0.165437	,	0.666773	,	-0.976304	,	-0.170024	,	0.905794	,\
		0.473908	,	-0.855725	,	-0.0413591	,	-0.508661	,	0.443453	,\
		0.842925	,	-0.144503	,	0.936699	,	-0.443935	,	-0.182996]	,\
		[0.803564	,	0.960386	,	-0.0323329	,	0.638181	,	-0.895684	,\
		-0.360502	,	0.0646258	,	-0.202449	,	-0.717228	,	0.970489	,\
		0.404608	,	-0.0861868	,	-0.879417	,	-0.866462	,	-0.938336	,\
		-0.799685	,	0.213464	,	-0.932344	,	-0.668236	,	0.751366	,\
		-0.22712	,	-0.407783	,	0.657463	,	0.0970092	,	-0.579109	,\
		-0.868866	,	-0.504041	,	0.926483	,	0.169978	,	-0.00842563],\
		[-0.530324	,	0.282745	,	0.0255867	,	0.287686	,	0.410417	,\
		-0.766576	,	-0.536566	,	-0.628288	,	0.69665		,	0.820713	,\
		-0.506162	,	-0.404114	,	0.640099	,	-0.956123	,	-0.576586	,\
		0.435502	,	-0.470676	,	-0.367062	,	-0.831765	,	-0.294942	,\
		0.518991	,	0.922338	,	0.337886	,	-0.67474	,	-0.725667	,\
		0.916684	,	0.39175		,	0.759081	,	0.496979	,	-0.200691]	,\
		[0.0417966	,	-0.687391	,	0.438773	,	0.287357	,	0.316636	,\
		-0.262311	,	-0.0755541	,	-0.442313	,	0.621378	,	0.670105	,\
		0.060982	,	0.944162	,	0.643442	,	-0.750684	,	-0.639973	,\
		0.217424	,	0.592823	,	0.929094	,	-0.239135	,	-0.41628	,\
		0.570893	,	-0.0798988	,	-0.917135	,	-0.749545	,	-0.982047	,\
		0.0626998	,	-0.977963	,	0.660401	,	0.470569	,	-0.0528868]	,\
		[-0.00138645	,	0.931065	,	-0.748519	,	0.304188	,	-0.266153,\
		0.672524	,	-0.105179	,	-0.874749	,	-0.154355	,	-0.774656	,\
		-0.69654	,	0.433098	,	0.615897	,	-0.387919	,	-0.429779	,\
		0.650202	,	0.122306	,	-0.237727	,	0.626817	,	-0.227929	,\
		0.405916	,	0.483328	,	0.282047	,	-0.262206	,	0.784123	,\
		0.83125		,	-0.662272	,	0.702768	,	0.875814	,	-0.701221]	,\
		[0.553793	,	0.471795	,	0.769147	,	0.059668	,	-0.841617	,\
		-0.191179	,	-0.972471	,	-0.825361	,	0.779826	,	-0.917201	,\
		0.43272		,	0.10301		,	0.358771	,	0.793448	,	-0.0379954	,\
		-0.870112	,	0.600442	,	-0.990603	,	0.549151	,	0.512146	,\
		-0.795843	,	0.490091	,	0.372046	,	-0.549437	,	0.0964285	,\
		0.753047	,	-0.86284	,	-0.589688	,	0.178612	,	-0.720358]])

		x0 = self.startp()

		self.name  = 'L2ZDT4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += (y[i]**2-10.0*np.cos(4.0*np.pi*y[i]))

		gx += 1.0 + 10.0*float(self.n-1)

		h = 1.0-np.sqrt(f1/gx)
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L2ZDT6():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507]	,\
		[0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617]	,\
		[0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749]	,\
		[-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599]	,\
		[-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477]	,\
		[0.238261	,	-0.1596		,	-0.827302	,	0.669248	,	0.494475	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086]	,\
		[-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176]	,\
		[-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423]	,\
		[0.60624		,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808]	,\
		[0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715]	\
		])

		x0 = self.startp()

		self.name  = 'L2ZDT6'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += (y[i]**2)

		gx = 1.0 + 9.0*(gx/float(self.n-1))**0.25

		h = 1.0-(f1/gx)**2
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L3ZDT1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.742377	,	0.426391	,	0.408202	,	0.633885	,	-0.0351053	,\
		0.712758	,	-0.191812	,	-0.390938	,	0.952828	,	0.921519]	,\
		[0.0969326	,	0.089775	,	-0.241157	,	0.0835558	,	-0.420236	,\
		0.00225943	,	0.0101814	,	0.441456	,	0.0633629	,	0.406631	,\
		0.507807	,	0.804148	,	0.963269	,	0.357128	,	-0.832565	,\
		-0.61366	,	-0.204783	,	-0.842476	,	-0.249524	,	-0.0985226	,\
		0.125225	,	0.487649	,	0.147046	,	0.679639	,	0.593707	,\
		0.24476	,	0.941339	,	-0.613783	,	0.402772	,	0.300775]	,\
		[-0.591523	,	-0.606614	,	-0.181873	,	0.692975	,	0.50208	,\
		0.0271269	,	0.804879	,	-0.402973	,	0.800373	,	0.760082	,\
		-0.753397	,	0.617418	,	0.689874	,	0.983384	,	0.668786	,\
		-0.798697	,	-0.244945	,	-0.942649	,	0.402856	,	-0.494672	,\
		-0.545421	,	-0.500243	,	0.154371	,	0.170017	,	-0.259108	,\
		0.956914	,	-0.0620912	,	0.634479	,	0.928617	,	0.464664]	,\
		[0.249008	,	0.370711	,	-0.633174	,	-0.0121906	,	0.42006	,\
		-0.630329	,	-0.763376	,	0.62538	,	0.818945	,	0.891598	,\
		-0.699445	,	0.237731	,	-0.324597	,	-0.800406	,	-0.42585	,\
		0.400334	,	-0.367816	,	0.198455	,	-0.983183	,	0.278976	,\
		0.57488	,	-0.361951	,	-0.0739728	,	0.91438	,	-0.391653	,\
		0.430805	,	0.706102	,	0.423887	,	0.296828	,	-0.265607]	,\
		[0.975863	,	-0.971371	,	-0.124115	,	0.4339	,	-0.254671	,\
		0.803564	,	0.960386	,	-0.0323329	,	0.638181	,	-0.895684	,\
		-0.530324	,	0.282745	,	0.0255867	,	0.287686	,	0.410417	,\
		0.0417966	,	-0.687391	,	0.438773	,	0.287357	,	0.316636	,\
		-0.00138645	,	0.931065	,	-0.748519	,	0.304188	,	-0.266153	,\
		0.553793	,	0.471795	,	0.769147	,	0.059668	,	-0.841617]	,\
		[0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715	,\
		-0.723444	,	-0.577654	,	0.0276004	,	0.0712472	,	-0.622791	,\
		0.923094	,	0.93011	,	-0.945394	,	-0.0934027	,	0.964123]	,\
		[-0.686633	,	-0.711276	,	-0.00325377	,	0.435196	,	-0.710002	,\
		-0.0100638	,	-0.177972	,	-0.491075	,	0.537035	,	-0.924987	,\
		-0.312441	,	0.327779	,	0.184745	,	0.246139	,	-0.936814	,\
		0.0671501	,	-0.527707	,	-0.509489	,	-0.883254	,	0.14851	,\
		-0.311828	,	-0.797099	,	-0.35815	,	0.95808	,	0.907244	,\
		-0.820314	,	-0.894233	,	-0.405896	,	0.0735439	,	0.486645]	,\
		[-0.536704	,	0.359652	,	0.839082	,	0.56817	,	-0.0776788	,\
		-0.878364	,	0.176801	,	-0.548932	,	-0.225601	,	-0.164912	,\
		0.0304653	,	-0.625221	,	-0.13318	,	0.827343	,	-0.101358	,\
		0.439941	,	-0.88216	,	0.170196	,	0.650734	,	-0.0982391	,\
		-0.868862	,	-0.50731	,	-0.848317	,	0.835712	,	0.616391	,\
		0.377022	,	0.63047	,	-0.198619	,	-0.576153	,	0.565373]	,\
		[0.169373	,	-0.975542	,	-0.0297852	,	0.80481	,	0.638317	,\
		0.680494	,	0.471868	,	-0.769787	,	-0.878099	,	-0.973724	,\
		-0.710739	,	-0.144068	,	-0.828545	,	-0.800912	,	0.184654	,\
		0.714817	,	0.307911	,	0.812861	,	-0.403497	,	-0.784382	,\
		0.0193887	,	0.412634	,	-0.169813	,	0.471794	,	0.660792	,\
		0.338806	,	-0.15829	,	0.642516	,	0.355126	,	0.174447]	,\
		[0.298068	,	-0.349803	,	-0.73185	,	0.488106	,	-0.0495073	,\
		-0.360502	,	0.0646258	,	-0.202449	,	-0.717228	,	0.970489	,\
		-0.766576	,	-0.536566	,	-0.628288	,	0.69665	,	0.820713	,\
		-0.262311	,	-0.0755541	,	-0.442313	,	0.621378	,	0.670105	,\
		0.672524	,	-0.105179	,	-0.874749	,	-0.154355	,	-0.774656	,\
		-0.191179	,	-0.972471	,	-0.825361	,	0.779826	,	-0.917201]	,\
		[0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.0850755	,	0.0419388	,	-0.323614	,	-0.973719	,	-0.680238	,\
		0.155451	,	0.442717	,	-0.792786	,	0.925785	,	0.670266	,\
		-0.795609	,	-0.289563	,	0.614236	,	-0.670585	,	0.466877]	,\
		[0.00283691	,	-0.168757	,	-0.134045	,	-0.655235	,	0.172361	,\
		-0.699424	,	0.742285	,	0.0181443	,	0.718971	,	-0.0308272	,\
		-0.931734	,	-0.0327827	,	0.319293	,	0.044473	,	-0.641645	,\
		-0.906465	,	0.496238	,	-0.853211	,	-0.779234	,	-0.979515	,\
		0.772426	,	0.720574	,	-0.873217	,	0.371431	,	-0.826029	,\
		-0.394355	,	0.125097	,	-0.316386	,	-0.701215	,	-0.845742]	,\
		[-0.00332785	,	0.459538	,	-0.518313	,	-0.270738	,	-0.629958	,\
		-0.208143	,	0.7768	,	-0.542743	,	-0.156021	,	0.671736	,\
		-0.999522	,	-0.0525574	,	-0.458319	,	0.587409	,	-0.334639	,\
		-0.468732	,	0.342133	,	-0.838071	,	-0.832362	,	0.658177	,\
		-0.442608	,	-0.158	,	0.313451	,	0.703748	,	-0.755984	,\
		-0.524245	,	-0.187299	,	-0.614524	,	0.429316	,	-0.491171]	,\
		[-0.670967	,	0.935792	,	-0.35605	,	0.175773	,	0.878601	,\
		0.354362	,	-0.1792	,	-0.225034	,	-0.44548	,	0.598865	,\
		-0.63675	,	-0.16696	,	0.240427	,	-0.513443	,	0.812664	,\
		-0.161823	,	-0.120835	,	0.323172	,	0.583739	,	0.732924	,\
		-0.350906	,	-0.612644	,	0.347876	,	0.112573	,	-0.501126	,\
		-0.975015	,	0.869905	,	-0.145285	,	-0.484002	,	-0.475966]	,\
		[0.253969	,	0.168116	,	0.148772	,	0.889593	,	-0.512213	,\
		0.404608	,	-0.0861868	,	-0.879417	,	-0.866462	,	-0.938336	,\
		-0.506162	,	-0.404114	,	0.640099	,	-0.956123	,	-0.576586	,\
		0.060982	,	0.944162	,	0.643442	,	-0.750684	,	-0.639973	,\
		-0.69654	,	0.433098	,	0.615897	,	-0.387919	,	-0.429779	,\
		0.43272	,	0.10301	,	0.358771	,	0.793448	,	-0.0379954]	,\
		[0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423	,\
		-0.270873	,	-0.209617	,	0.968436	,	0.908798	,	0.975851	,\
		-0.865566	,	-0.638281	,	0.333094	,	0.477628	,	0.47261	,\
		0.144597	,	-0.206416	,	0.6937	,	-0.967958	,	-0.0951247]	,\
		[0.998291	,	0.376291	,	-0.962215	,	-0.363174	,	-0.88777	,\
		0.086931	,	0.524476	,	0.956457	,	0.143024	,	0.616481	,\
		0.596118	,	-0.293934	,	-0.63373	,	0.409658	,	0.759892	,\
		0.827175	,	0.228969	,	-0.402829	,	-0.970118	,	0.762559	,\
		0.942716	,	0.70609	,	-0.658158	,	-0.782185	,	-0.806743	,\
		0.2065	,	-0.413743	,	0.406725	,	-0.423813	,	-0.941255]	,\
		[-0.755084	,	-0.721573	,	0.431107	,	-0.221877	,	0.32543	,\
		0.878648	,	-0.419588	,	-0.0752896	,	0.0299447	,	-0.494459	,\
		0.0759643	,	0.0255827	,	0.128944	,	0.17317	,	-0.284309	,\
		-0.565361	,	0.149473	,	0.69331	,	-0.491848	,	0.74916	,\
		-0.249443	,	0.491564	,	0.985068	,	0.678644	,	0.808324	,\
		0.399495	,	-0.333898	,	-0.646636	,	-0.0189709	,	-0.339605]	,\
		[-0.275168	,	-0.932517	,	-0.372497	,	-0.0732907	,	-0.185493	,\
		0.544005	,	-0.478962	,	0.327193	,	-0.525784	,	0.903179	,\
		0.744943	,	0.970612	,	0.00172899	,	-0.726378	,	-0.0985012	,\
		-0.220603	,	-0.594121	,	0.935093	,	-0.216736	,	0.659318	,\
		0.456761	,	-0.109004	,	0.289352	,	-0.566504	,	0.585042	,\
		-0.67704	,	0.996452	,	-0.0685748	,	-0.851985	,	0.416498]	,\
		[-0.165437	,	0.666773	,	-0.976304	,	-0.170024	,	0.905794	,\
		-0.799685	,	0.213464	,	-0.932344	,	-0.668236	,	0.751366	,\
		0.435502	,	-0.470676	,	-0.367062	,	-0.831765	,	-0.294942	,\
		0.217424	,	0.592823	,	0.929094	,	-0.239135	,	-0.41628	,\
		0.650202	,	0.122306	,	-0.237727	,	0.626817	,	-0.227929	,\
		-0.870112	,	0.600442	,	-0.990603	,	0.549151	,	0.512146]	,\
		[0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		0.238261	,	-0.1596	,	-0.827302	,	0.669248	,	0.494475	,\
		0.60624	,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581	,\
		-0.994918	,	-0.0621977	,	0.628171	,	0.761228	,	0.34372	,\
		0.23151	,	0.82132	,	-0.589803	,	0.796275	,	0.57713	,\
		-0.942473	,	-0.610767	,	-0.655472	,	-0.0960986	,	-0.302779]	,\
		[-0.519929	,	-0.560554	,	-0.984415	,	0.601529	,	-0.984103	,\
		0.217909	,	-0.128427	,	-0.262427	,	-0.938208	,	-0.52479	,\
		-0.257078	,	0.939616	,	-0.227661	,	0.115754	,	0.10964	,\
		0.506495	,	0.460303	,	0.897304	,	0.686003	,	0.739986	,\
		-0.627986	,	-0.405551	,	-0.258495	,	-0.796524	,	0.222498	,\
		-0.558804	,	0.312326	,	0.345314	,	0.319143	,	-0.644653]	,\
		[0.163743	,	0.0759916	,	0.695064	,	-0.656856	,	0.074163	,\
		-0.72415	,	0.35978	,	-0.32646	,	-0.96605	,	0.0127605	,\
		0.287161	,	-0.550725	,	-0.433083	,	-0.242821	,	0.878879	,\
		0.526025	,	-0.155339	,	0.0998096	,	0.468761	,	0.324649	,\
		0.81975	,	-0.435823	,	-0.839855	,	0.00282368	,	-0.569165	,\
		-0.798791	,	0.0494081	,	0.367012	,	0.852545	,	0.43557]	,\
		[-0.357004	,	0.314786	,	-0.229239	,	0.530256	,	-0.51327	,\
		-0.899248	,	0.156514	,	0.154329	,	0.499808	,	-0.836327	,\
		0.224232	,	0.16495	,	0.560077	,	-0.813112	,	0.112894	,\
		-0.750417	,	-0.284773	,	-0.271496	,	0.491731	,	-0.712174	,\
		0.584934	,	0.923676	,	0.895312	,	-0.161036	,	-0.995895	,\
		0.791047	,	-0.211323	,	-0.302819	,	0.640735	,	-0.317908]	,\
		[0.473908	,	-0.855725	,	-0.0413591	,	-0.508661	,	0.443453	,\
		-0.22712	,	-0.407783	,	0.657463	,	0.0970092	,	-0.579109	,\
		0.518991	,	0.922338	,	0.337886	,	-0.67474	,	-0.725667	,\
		0.570893	,	-0.0798988	,	-0.917135	,	-0.749545	,	-0.982047	,\
		0.405916	,	0.483328	,	0.282047	,	-0.262206	,	0.784123	,\
		-0.795843	,	0.490091	,	0.372046	,	-0.549437	,	0.0964285]	,\
		[-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086	,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808	,\
		-0.792042	,	-0.144765	,	-0.965748	,	0.0133606	,	-0.0260565	,\
		0.101149	,	0.970191	,	0.532821	,	0.814769	,	-0.0687269	,\
		-0.734976	,	-0.342188	,	-0.315861	,	-0.912834	,	0.24499]	,\
		[-0.228237	,	-0.578066	,	0.307023	,	0.606123	,	0.959635	,\
		0.12919	,	0.721925	,	0.766492	,	0.470845	,	-0.0976466	,\
		-0.240557	,	0.66842	,	0.855535	,	-0.451536	,	0.264961	,\
		0.15731	,	0.281697	,	-0.922955	,	-0.780824	,	0.449716	,\
		0.087545	,	-0.0917108	,	-0.62542	,	-0.110256	,	0.0417576	,\
		-0.0408415	,	0.176461	,	0.740113	,	0.470737	,	-0.914927]	,\
		[0.264319	,	-0.73174	,	0.731548	,	-0.489341	,	0.678946	,\
		0.563174	,	-0.814853	,	-0.949609	,	-0.526794	,	-0.801902	,\
		0.691699	,	-0.660499	,	0.389985	,	0.599856	,	-0.711442	,\
		0.128488	,	0.544144	,	-0.495222	,	0.965229	,	-0.79314	,\
		0.0884339	,	-0.222144	,	0.499412	,	-0.565198	,	0.64824	,\
		0.150039	,	-0.0454542	,	0.604861	,	-0.598288	,	-0.500696]	,\
		[0.44187	,	0.940309	,	-0.240334	,	-0.0276121	,	0.74383	,\
		-0.802627	,	0.378082	,	-0.112673	,	-0.47926	,	-0.3355	,\
		-0.0955366	,	0.0187107	,	0.913887	,	0.123076	,	0.550338	,\
		-0.763681	,	0.0781023	,	0.951666	,	0.734031	,	0.826912	,\
		0.0853141	,	-0.583368	,	-0.157612	,	0.234119	,	0.875043	,\
		-0.116586	,	-0.896382	,	-0.817317	,	-0.948837	,	-0.597427]	,\
		[0.842925	,	-0.144503	,	0.936699	,	-0.443935	,	-0.182996	,\
		-0.868866	,	-0.504041	,	0.926483	,	0.169978	,	-0.00842563	,\
		0.916684	,	0.39175	,	0.759081	,	0.496979	,	-0.200691	,\
		0.0626998	,	-0.977963	,	0.660401	,	0.470569	,	-0.0528868	,\
		0.83125	,	-0.662272	,	0.702768	,	0.875814	,	-0.701221	, \
		0.753047	,	-0.86284	,	-0.589688	,	0.178612	,	-0.720358]
		])

		x0 = self.startp()

		self.name  = 'L3ZDT1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]**2

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += (y[i]**2)

		gx = 1.0 + 9.0/float(self.n-1) * gx

		h = 1.0-np.sqrt(f1/gx)
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L3ZDT2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.742377	,	0.426391	,	0.408202	,	0.633885	,	-0.0351053	,\
		0.712758	,	-0.191812	,	-0.390938	,	0.952828	,	0.921519]	,\
		[0.0969326	,	0.089775	,	-0.241157	,	0.0835558	,	-0.420236	,\
		0.00225943	,	0.0101814	,	0.441456	,	0.0633629	,	0.406631	,\
		0.507807	,	0.804148	,	0.963269	,	0.357128	,	-0.832565	,\
		-0.61366	,	-0.204783	,	-0.842476	,	-0.249524	,	-0.0985226	,\
		0.125225	,	0.487649	,	0.147046	,	0.679639	,	0.593707	,\
		0.24476	,	0.941339	,	-0.613783	,	0.402772	,	0.300775]	,\
		[-0.591523	,	-0.606614	,	-0.181873	,	0.692975	,	0.50208	,\
		0.0271269	,	0.804879	,	-0.402973	,	0.800373	,	0.760082	,\
		-0.753397	,	0.617418	,	0.689874	,	0.983384	,	0.668786	,\
		-0.798697	,	-0.244945	,	-0.942649	,	0.402856	,	-0.494672	,\
		-0.545421	,	-0.500243	,	0.154371	,	0.170017	,	-0.259108	,\
		0.956914	,	-0.0620912	,	0.634479	,	0.928617	,	0.464664]	,\
		[0.249008	,	0.370711	,	-0.633174	,	-0.0121906	,	0.42006	,\
		-0.630329	,	-0.763376	,	0.62538	,	0.818945	,	0.891598	,\
		-0.699445	,	0.237731	,	-0.324597	,	-0.800406	,	-0.42585	,\
		0.400334	,	-0.367816	,	0.198455	,	-0.983183	,	0.278976	,\
		0.57488	,	-0.361951	,	-0.0739728	,	0.91438	,	-0.391653	,\
		0.430805	,	0.706102	,	0.423887	,	0.296828	,	-0.265607]	,\
		[0.975863	,	-0.971371	,	-0.124115	,	0.4339	,	-0.254671	,\
		0.803564	,	0.960386	,	-0.0323329	,	0.638181	,	-0.895684	,\
		-0.530324	,	0.282745	,	0.0255867	,	0.287686	,	0.410417	,\
		0.0417966	,	-0.687391	,	0.438773	,	0.287357	,	0.316636	,\
		-0.00138645	,	0.931065	,	-0.748519	,	0.304188	,	-0.266153	,\
		0.553793	,	0.471795	,	0.769147	,	0.059668	,	-0.841617]	,\
		[0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715	,\
		-0.723444	,	-0.577654	,	0.0276004	,	0.0712472	,	-0.622791	,\
		0.923094	,	0.93011	,	-0.945394	,	-0.0934027	,	0.964123]	,\
		[-0.686633	,	-0.711276	,	-0.00325377	,	0.435196	,	-0.710002	,\
		-0.0100638	,	-0.177972	,	-0.491075	,	0.537035	,	-0.924987	,\
		-0.312441	,	0.327779	,	0.184745	,	0.246139	,	-0.936814	,\
		0.0671501	,	-0.527707	,	-0.509489	,	-0.883254	,	0.14851	,\
		-0.311828	,	-0.797099	,	-0.35815	,	0.95808	,	0.907244	,\
		-0.820314	,	-0.894233	,	-0.405896	,	0.0735439	,	0.486645]	,\
		[-0.536704	,	0.359652	,	0.839082	,	0.56817	,	-0.0776788	,\
		-0.878364	,	0.176801	,	-0.548932	,	-0.225601	,	-0.164912	,\
		0.0304653	,	-0.625221	,	-0.13318	,	0.827343	,	-0.101358	,\
		0.439941	,	-0.88216	,	0.170196	,	0.650734	,	-0.0982391	,\
		-0.868862	,	-0.50731	,	-0.848317	,	0.835712	,	0.616391	,\
		0.377022	,	0.63047	,	-0.198619	,	-0.576153	,	0.565373]	,\
		[0.169373	,	-0.975542	,	-0.0297852	,	0.80481	,	0.638317	,\
		0.680494	,	0.471868	,	-0.769787	,	-0.878099	,	-0.973724	,\
		-0.710739	,	-0.144068	,	-0.828545	,	-0.800912	,	0.184654	,\
		0.714817	,	0.307911	,	0.812861	,	-0.403497	,	-0.784382	,\
		0.0193887	,	0.412634	,	-0.169813	,	0.471794	,	0.660792	,\
		0.338806	,	-0.15829	,	0.642516	,	0.355126	,	0.174447]	,\
		[0.298068	,	-0.349803	,	-0.73185	,	0.488106	,	-0.0495073	,\
		-0.360502	,	0.0646258	,	-0.202449	,	-0.717228	,	0.970489	,\
		-0.766576	,	-0.536566	,	-0.628288	,	0.69665	,	0.820713	,\
		-0.262311	,	-0.0755541	,	-0.442313	,	0.621378	,	0.670105	,\
		0.672524	,	-0.105179	,	-0.874749	,	-0.154355	,	-0.774656	,\
		-0.191179	,	-0.972471	,	-0.825361	,	0.779826	,	-0.917201]	,\
		[0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.0850755	,	0.0419388	,	-0.323614	,	-0.973719	,	-0.680238	,\
		0.155451	,	0.442717	,	-0.792786	,	0.925785	,	0.670266	,\
		-0.795609	,	-0.289563	,	0.614236	,	-0.670585	,	0.466877]	,\
		[0.00283691	,	-0.168757	,	-0.134045	,	-0.655235	,	0.172361	,\
		-0.699424	,	0.742285	,	0.0181443	,	0.718971	,	-0.0308272	,\
		-0.931734	,	-0.0327827	,	0.319293	,	0.044473	,	-0.641645	,\
		-0.906465	,	0.496238	,	-0.853211	,	-0.779234	,	-0.979515	,\
		0.772426	,	0.720574	,	-0.873217	,	0.371431	,	-0.826029	,\
		-0.394355	,	0.125097	,	-0.316386	,	-0.701215	,	-0.845742]	,\
		[-0.00332785	,	0.459538	,	-0.518313	,	-0.270738	,	-0.629958	,\
		-0.208143	,	0.7768	,	-0.542743	,	-0.156021	,	0.671736	,\
		-0.999522	,	-0.0525574	,	-0.458319	,	0.587409	,	-0.334639	,\
		-0.468732	,	0.342133	,	-0.838071	,	-0.832362	,	0.658177	,\
		-0.442608	,	-0.158	,	0.313451	,	0.703748	,	-0.755984	,\
		-0.524245	,	-0.187299	,	-0.614524	,	0.429316	,	-0.491171]	,\
		[-0.670967	,	0.935792	,	-0.35605	,	0.175773	,	0.878601	,\
		0.354362	,	-0.1792	,	-0.225034	,	-0.44548	,	0.598865	,\
		-0.63675	,	-0.16696	,	0.240427	,	-0.513443	,	0.812664	,\
		-0.161823	,	-0.120835	,	0.323172	,	0.583739	,	0.732924	,\
		-0.350906	,	-0.612644	,	0.347876	,	0.112573	,	-0.501126	,\
		-0.975015	,	0.869905	,	-0.145285	,	-0.484002	,	-0.475966]	,\
		[0.253969	,	0.168116	,	0.148772	,	0.889593	,	-0.512213	,\
		0.404608	,	-0.0861868	,	-0.879417	,	-0.866462	,	-0.938336	,\
		-0.506162	,	-0.404114	,	0.640099	,	-0.956123	,	-0.576586	,\
		0.060982	,	0.944162	,	0.643442	,	-0.750684	,	-0.639973	,\
		-0.69654	,	0.433098	,	0.615897	,	-0.387919	,	-0.429779	,\
		0.43272	,	0.10301	,	0.358771	,	0.793448	,	-0.0379954]	,\
		[0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423	,\
		-0.270873	,	-0.209617	,	0.968436	,	0.908798	,	0.975851	,\
		-0.865566	,	-0.638281	,	0.333094	,	0.477628	,	0.47261	,\
		0.144597	,	-0.206416	,	0.6937	,	-0.967958	,	-0.0951247]	,\
		[0.998291	,	0.376291	,	-0.962215	,	-0.363174	,	-0.88777	,\
		0.086931	,	0.524476	,	0.956457	,	0.143024	,	0.616481	,\
		0.596118	,	-0.293934	,	-0.63373	,	0.409658	,	0.759892	,\
		0.827175	,	0.228969	,	-0.402829	,	-0.970118	,	0.762559	,\
		0.942716	,	0.70609	,	-0.658158	,	-0.782185	,	-0.806743	,\
		0.2065	,	-0.413743	,	0.406725	,	-0.423813	,	-0.941255]	,\
		[-0.755084	,	-0.721573	,	0.431107	,	-0.221877	,	0.32543	,\
		0.878648	,	-0.419588	,	-0.0752896	,	0.0299447	,	-0.494459	,\
		0.0759643	,	0.0255827	,	0.128944	,	0.17317	,	-0.284309	,\
		-0.565361	,	0.149473	,	0.69331	,	-0.491848	,	0.74916	,\
		-0.249443	,	0.491564	,	0.985068	,	0.678644	,	0.808324	,\
		0.399495	,	-0.333898	,	-0.646636	,	-0.0189709	,	-0.339605]	,\
		[-0.275168	,	-0.932517	,	-0.372497	,	-0.0732907	,	-0.185493	,\
		0.544005	,	-0.478962	,	0.327193	,	-0.525784	,	0.903179	,\
		0.744943	,	0.970612	,	0.00172899	,	-0.726378	,	-0.0985012	,\
		-0.220603	,	-0.594121	,	0.935093	,	-0.216736	,	0.659318	,\
		0.456761	,	-0.109004	,	0.289352	,	-0.566504	,	0.585042	,\
		-0.67704	,	0.996452	,	-0.0685748	,	-0.851985	,	0.416498]	,\
		[-0.165437	,	0.666773	,	-0.976304	,	-0.170024	,	0.905794	,\
		-0.799685	,	0.213464	,	-0.932344	,	-0.668236	,	0.751366	,\
		0.435502	,	-0.470676	,	-0.367062	,	-0.831765	,	-0.294942	,\
		0.217424	,	0.592823	,	0.929094	,	-0.239135	,	-0.41628	,\
		0.650202	,	0.122306	,	-0.237727	,	0.626817	,	-0.227929	,\
		-0.870112	,	0.600442	,	-0.990603	,	0.549151	,	0.512146]	,\
		[0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		0.238261	,	-0.1596	,	-0.827302	,	0.669248	,	0.494475	,\
		0.60624	,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581	,\
		-0.994918	,	-0.0621977	,	0.628171	,	0.761228	,	0.34372	,\
		0.23151	,	0.82132	,	-0.589803	,	0.796275	,	0.57713	,\
		-0.942473	,	-0.610767	,	-0.655472	,	-0.0960986	,	-0.302779]	,\
		[-0.519929	,	-0.560554	,	-0.984415	,	0.601529	,	-0.984103	,\
		0.217909	,	-0.128427	,	-0.262427	,	-0.938208	,	-0.52479	,\
		-0.257078	,	0.939616	,	-0.227661	,	0.115754	,	0.10964	,\
		0.506495	,	0.460303	,	0.897304	,	0.686003	,	0.739986	,\
		-0.627986	,	-0.405551	,	-0.258495	,	-0.796524	,	0.222498	,\
		-0.558804	,	0.312326	,	0.345314	,	0.319143	,	-0.644653]	,\
		[0.163743	,	0.0759916	,	0.695064	,	-0.656856	,	0.074163	,\
		-0.72415	,	0.35978	,	-0.32646	,	-0.96605	,	0.0127605	,\
		0.287161	,	-0.550725	,	-0.433083	,	-0.242821	,	0.878879	,\
		0.526025	,	-0.155339	,	0.0998096	,	0.468761	,	0.324649	,\
		0.81975	,	-0.435823	,	-0.839855	,	0.00282368	,	-0.569165	,\
		-0.798791	,	0.0494081	,	0.367012	,	0.852545	,	0.43557]	,\
		[-0.357004	,	0.314786	,	-0.229239	,	0.530256	,	-0.51327	,\
		-0.899248	,	0.156514	,	0.154329	,	0.499808	,	-0.836327	,\
		0.224232	,	0.16495	,	0.560077	,	-0.813112	,	0.112894	,\
		-0.750417	,	-0.284773	,	-0.271496	,	0.491731	,	-0.712174	,\
		0.584934	,	0.923676	,	0.895312	,	-0.161036	,	-0.995895	,\
		0.791047	,	-0.211323	,	-0.302819	,	0.640735	,	-0.317908]	,\
		[0.473908	,	-0.855725	,	-0.0413591	,	-0.508661	,	0.443453	,\
		-0.22712	,	-0.407783	,	0.657463	,	0.0970092	,	-0.579109	,\
		0.518991	,	0.922338	,	0.337886	,	-0.67474	,	-0.725667	,\
		0.570893	,	-0.0798988	,	-0.917135	,	-0.749545	,	-0.982047	,\
		0.405916	,	0.483328	,	0.282047	,	-0.262206	,	0.784123	,\
		-0.795843	,	0.490091	,	0.372046	,	-0.549437	,	0.0964285]	,\
		[-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086	,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808	,\
		-0.792042	,	-0.144765	,	-0.965748	,	0.0133606	,	-0.0260565	,\
		0.101149	,	0.970191	,	0.532821	,	0.814769	,	-0.0687269	,\
		-0.734976	,	-0.342188	,	-0.315861	,	-0.912834	,	0.24499]	,\
		[-0.228237	,	-0.578066	,	0.307023	,	0.606123	,	0.959635	,\
		0.12919	,	0.721925	,	0.766492	,	0.470845	,	-0.0976466	,\
		-0.240557	,	0.66842	,	0.855535	,	-0.451536	,	0.264961	,\
		0.15731	,	0.281697	,	-0.922955	,	-0.780824	,	0.449716	,\
		0.087545	,	-0.0917108	,	-0.62542	,	-0.110256	,	0.0417576	,\
		-0.0408415	,	0.176461	,	0.740113	,	0.470737	,	-0.914927]	,\
		[0.264319	,	-0.73174	,	0.731548	,	-0.489341	,	0.678946	,\
		0.563174	,	-0.814853	,	-0.949609	,	-0.526794	,	-0.801902	,\
		0.691699	,	-0.660499	,	0.389985	,	0.599856	,	-0.711442	,\
		0.128488	,	0.544144	,	-0.495222	,	0.965229	,	-0.79314	,\
		0.0884339	,	-0.222144	,	0.499412	,	-0.565198	,	0.64824	,\
		0.150039	,	-0.0454542	,	0.604861	,	-0.598288	,	-0.500696]	,\
		[0.44187	,	0.940309	,	-0.240334	,	-0.0276121	,	0.74383	,\
		-0.802627	,	0.378082	,	-0.112673	,	-0.47926	,	-0.3355	,\
		-0.0955366	,	0.0187107	,	0.913887	,	0.123076	,	0.550338	,\
		-0.763681	,	0.0781023	,	0.951666	,	0.734031	,	0.826912	,\
		0.0853141	,	-0.583368	,	-0.157612	,	0.234119	,	0.875043	,\
		-0.116586	,	-0.896382	,	-0.817317	,	-0.948837	,	-0.597427]	,\
		[0.842925	,	-0.144503	,	0.936699	,	-0.443935	,	-0.182996	,\
		-0.868866	,	-0.504041	,	0.926483	,	0.169978	,	-0.00842563	,\
		0.916684	,	0.39175	,	0.759081	,	0.496979	,	-0.200691	,\
		0.0626998	,	-0.977963	,	0.660401	,	0.470569	,	-0.0528868	,\
		0.83125	,	-0.662272	,	0.702768	,	0.875814	,	-0.701221	, \
		0.753047	,	-0.86284	,	-0.589688	,	0.178612	,	-0.720358]
		])

		x0 = self.startp()

		self.name  = 'L3ZDT2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]**2

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += (y[i]**2)

		gx = 1.0 + 9.0/float(self.n-1) * gx

		h = 1.0-(f1/gx)**2
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L3ZDT3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.742377	,	0.426391	,	0.408202	,	0.633885	,	-0.0351053	,\
		0.712758	,	-0.191812	,	-0.390938	,	0.952828	,	0.921519]	,\
		[0.0969326	,	0.089775	,	-0.241157	,	0.0835558	,	-0.420236	,\
		0.00225943	,	0.0101814	,	0.441456	,	0.0633629	,	0.406631	,\
		0.507807	,	0.804148	,	0.963269	,	0.357128	,	-0.832565	,\
		-0.61366	,	-0.204783	,	-0.842476	,	-0.249524	,	-0.0985226	,\
		0.125225	,	0.487649	,	0.147046	,	0.679639	,	0.593707	,\
		0.24476	,	0.941339	,	-0.613783	,	0.402772	,	0.300775]	,\
		[-0.591523	,	-0.606614	,	-0.181873	,	0.692975	,	0.50208	,\
		0.0271269	,	0.804879	,	-0.402973	,	0.800373	,	0.760082	,\
		-0.753397	,	0.617418	,	0.689874	,	0.983384	,	0.668786	,\
		-0.798697	,	-0.244945	,	-0.942649	,	0.402856	,	-0.494672	,\
		-0.545421	,	-0.500243	,	0.154371	,	0.170017	,	-0.259108	,\
		0.956914	,	-0.0620912	,	0.634479	,	0.928617	,	0.464664]	,\
		[0.249008	,	0.370711	,	-0.633174	,	-0.0121906	,	0.42006	,\
		-0.630329	,	-0.763376	,	0.62538	,	0.818945	,	0.891598	,\
		-0.699445	,	0.237731	,	-0.324597	,	-0.800406	,	-0.42585	,\
		0.400334	,	-0.367816	,	0.198455	,	-0.983183	,	0.278976	,\
		0.57488	,	-0.361951	,	-0.0739728	,	0.91438	,	-0.391653	,\
		0.430805	,	0.706102	,	0.423887	,	0.296828	,	-0.265607]	,\
		[0.975863	,	-0.971371	,	-0.124115	,	0.4339	,	-0.254671	,\
		0.803564	,	0.960386	,	-0.0323329	,	0.638181	,	-0.895684	,\
		-0.530324	,	0.282745	,	0.0255867	,	0.287686	,	0.410417	,\
		0.0417966	,	-0.687391	,	0.438773	,	0.287357	,	0.316636	,\
		-0.00138645	,	0.931065	,	-0.748519	,	0.304188	,	-0.266153	,\
		0.553793	,	0.471795	,	0.769147	,	0.059668	,	-0.841617]	,\
		[0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715	,\
		-0.723444	,	-0.577654	,	0.0276004	,	0.0712472	,	-0.622791	,\
		0.923094	,	0.93011	,	-0.945394	,	-0.0934027	,	0.964123]	,\
		[-0.686633	,	-0.711276	,	-0.00325377	,	0.435196	,	-0.710002	,\
		-0.0100638	,	-0.177972	,	-0.491075	,	0.537035	,	-0.924987	,\
		-0.312441	,	0.327779	,	0.184745	,	0.246139	,	-0.936814	,\
		0.0671501	,	-0.527707	,	-0.509489	,	-0.883254	,	0.14851	,\
		-0.311828	,	-0.797099	,	-0.35815	,	0.95808	,	0.907244	,\
		-0.820314	,	-0.894233	,	-0.405896	,	0.0735439	,	0.486645]	,\
		[-0.536704	,	0.359652	,	0.839082	,	0.56817	,	-0.0776788	,\
		-0.878364	,	0.176801	,	-0.548932	,	-0.225601	,	-0.164912	,\
		0.0304653	,	-0.625221	,	-0.13318	,	0.827343	,	-0.101358	,\
		0.439941	,	-0.88216	,	0.170196	,	0.650734	,	-0.0982391	,\
		-0.868862	,	-0.50731	,	-0.848317	,	0.835712	,	0.616391	,\
		0.377022	,	0.63047	,	-0.198619	,	-0.576153	,	0.565373]	,\
		[0.169373	,	-0.975542	,	-0.0297852	,	0.80481	,	0.638317	,\
		0.680494	,	0.471868	,	-0.769787	,	-0.878099	,	-0.973724	,\
		-0.710739	,	-0.144068	,	-0.828545	,	-0.800912	,	0.184654	,\
		0.714817	,	0.307911	,	0.812861	,	-0.403497	,	-0.784382	,\
		0.0193887	,	0.412634	,	-0.169813	,	0.471794	,	0.660792	,\
		0.338806	,	-0.15829	,	0.642516	,	0.355126	,	0.174447]	,\
		[0.298068	,	-0.349803	,	-0.73185	,	0.488106	,	-0.0495073	,\
		-0.360502	,	0.0646258	,	-0.202449	,	-0.717228	,	0.970489	,\
		-0.766576	,	-0.536566	,	-0.628288	,	0.69665	,	0.820713	,\
		-0.262311	,	-0.0755541	,	-0.442313	,	0.621378	,	0.670105	,\
		0.672524	,	-0.105179	,	-0.874749	,	-0.154355	,	-0.774656	,\
		-0.191179	,	-0.972471	,	-0.825361	,	0.779826	,	-0.917201]	,\
		[0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.0850755	,	0.0419388	,	-0.323614	,	-0.973719	,	-0.680238	,\
		0.155451	,	0.442717	,	-0.792786	,	0.925785	,	0.670266	,\
		-0.795609	,	-0.289563	,	0.614236	,	-0.670585	,	0.466877]	,\
		[0.00283691	,	-0.168757	,	-0.134045	,	-0.655235	,	0.172361	,\
		-0.699424	,	0.742285	,	0.0181443	,	0.718971	,	-0.0308272	,\
		-0.931734	,	-0.0327827	,	0.319293	,	0.044473	,	-0.641645	,\
		-0.906465	,	0.496238	,	-0.853211	,	-0.779234	,	-0.979515	,\
		0.772426	,	0.720574	,	-0.873217	,	0.371431	,	-0.826029	,\
		-0.394355	,	0.125097	,	-0.316386	,	-0.701215	,	-0.845742]	,\
		[-0.00332785	,	0.459538	,	-0.518313	,	-0.270738	,	-0.629958	,\
		-0.208143	,	0.7768	,	-0.542743	,	-0.156021	,	0.671736	,\
		-0.999522	,	-0.0525574	,	-0.458319	,	0.587409	,	-0.334639	,\
		-0.468732	,	0.342133	,	-0.838071	,	-0.832362	,	0.658177	,\
		-0.442608	,	-0.158	,	0.313451	,	0.703748	,	-0.755984	,\
		-0.524245	,	-0.187299	,	-0.614524	,	0.429316	,	-0.491171]	,\
		[-0.670967	,	0.935792	,	-0.35605	,	0.175773	,	0.878601	,\
		0.354362	,	-0.1792	,	-0.225034	,	-0.44548	,	0.598865	,\
		-0.63675	,	-0.16696	,	0.240427	,	-0.513443	,	0.812664	,\
		-0.161823	,	-0.120835	,	0.323172	,	0.583739	,	0.732924	,\
		-0.350906	,	-0.612644	,	0.347876	,	0.112573	,	-0.501126	,\
		-0.975015	,	0.869905	,	-0.145285	,	-0.484002	,	-0.475966]	,\
		[0.253969	,	0.168116	,	0.148772	,	0.889593	,	-0.512213	,\
		0.404608	,	-0.0861868	,	-0.879417	,	-0.866462	,	-0.938336	,\
		-0.506162	,	-0.404114	,	0.640099	,	-0.956123	,	-0.576586	,\
		0.060982	,	0.944162	,	0.643442	,	-0.750684	,	-0.639973	,\
		-0.69654	,	0.433098	,	0.615897	,	-0.387919	,	-0.429779	,\
		0.43272	,	0.10301	,	0.358771	,	0.793448	,	-0.0379954]	,\
		[0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423	,\
		-0.270873	,	-0.209617	,	0.968436	,	0.908798	,	0.975851	,\
		-0.865566	,	-0.638281	,	0.333094	,	0.477628	,	0.47261	,\
		0.144597	,	-0.206416	,	0.6937	,	-0.967958	,	-0.0951247]	,\
		[0.998291	,	0.376291	,	-0.962215	,	-0.363174	,	-0.88777	,\
		0.086931	,	0.524476	,	0.956457	,	0.143024	,	0.616481	,\
		0.596118	,	-0.293934	,	-0.63373	,	0.409658	,	0.759892	,\
		0.827175	,	0.228969	,	-0.402829	,	-0.970118	,	0.762559	,\
		0.942716	,	0.70609	,	-0.658158	,	-0.782185	,	-0.806743	,\
		0.2065	,	-0.413743	,	0.406725	,	-0.423813	,	-0.941255]	,\
		[-0.755084	,	-0.721573	,	0.431107	,	-0.221877	,	0.32543	,\
		0.878648	,	-0.419588	,	-0.0752896	,	0.0299447	,	-0.494459	,\
		0.0759643	,	0.0255827	,	0.128944	,	0.17317	,	-0.284309	,\
		-0.565361	,	0.149473	,	0.69331	,	-0.491848	,	0.74916	,\
		-0.249443	,	0.491564	,	0.985068	,	0.678644	,	0.808324	,\
		0.399495	,	-0.333898	,	-0.646636	,	-0.0189709	,	-0.339605]	,\
		[-0.275168	,	-0.932517	,	-0.372497	,	-0.0732907	,	-0.185493	,\
		0.544005	,	-0.478962	,	0.327193	,	-0.525784	,	0.903179	,\
		0.744943	,	0.970612	,	0.00172899	,	-0.726378	,	-0.0985012	,\
		-0.220603	,	-0.594121	,	0.935093	,	-0.216736	,	0.659318	,\
		0.456761	,	-0.109004	,	0.289352	,	-0.566504	,	0.585042	,\
		-0.67704	,	0.996452	,	-0.0685748	,	-0.851985	,	0.416498]	,\
		[-0.165437	,	0.666773	,	-0.976304	,	-0.170024	,	0.905794	,\
		-0.799685	,	0.213464	,	-0.932344	,	-0.668236	,	0.751366	,\
		0.435502	,	-0.470676	,	-0.367062	,	-0.831765	,	-0.294942	,\
		0.217424	,	0.592823	,	0.929094	,	-0.239135	,	-0.41628	,\
		0.650202	,	0.122306	,	-0.237727	,	0.626817	,	-0.227929	,\
		-0.870112	,	0.600442	,	-0.990603	,	0.549151	,	0.512146]	,\
		[0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		0.238261	,	-0.1596	,	-0.827302	,	0.669248	,	0.494475	,\
		0.60624	,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581	,\
		-0.994918	,	-0.0621977	,	0.628171	,	0.761228	,	0.34372	,\
		0.23151	,	0.82132	,	-0.589803	,	0.796275	,	0.57713	,\
		-0.942473	,	-0.610767	,	-0.655472	,	-0.0960986	,	-0.302779]	,\
		[-0.519929	,	-0.560554	,	-0.984415	,	0.601529	,	-0.984103	,\
		0.217909	,	-0.128427	,	-0.262427	,	-0.938208	,	-0.52479	,\
		-0.257078	,	0.939616	,	-0.227661	,	0.115754	,	0.10964	,\
		0.506495	,	0.460303	,	0.897304	,	0.686003	,	0.739986	,\
		-0.627986	,	-0.405551	,	-0.258495	,	-0.796524	,	0.222498	,\
		-0.558804	,	0.312326	,	0.345314	,	0.319143	,	-0.644653]	,\
		[0.163743	,	0.0759916	,	0.695064	,	-0.656856	,	0.074163	,\
		-0.72415	,	0.35978	,	-0.32646	,	-0.96605	,	0.0127605	,\
		0.287161	,	-0.550725	,	-0.433083	,	-0.242821	,	0.878879	,\
		0.526025	,	-0.155339	,	0.0998096	,	0.468761	,	0.324649	,\
		0.81975	,	-0.435823	,	-0.839855	,	0.00282368	,	-0.569165	,\
		-0.798791	,	0.0494081	,	0.367012	,	0.852545	,	0.43557]	,\
		[-0.357004	,	0.314786	,	-0.229239	,	0.530256	,	-0.51327	,\
		-0.899248	,	0.156514	,	0.154329	,	0.499808	,	-0.836327	,\
		0.224232	,	0.16495	,	0.560077	,	-0.813112	,	0.112894	,\
		-0.750417	,	-0.284773	,	-0.271496	,	0.491731	,	-0.712174	,\
		0.584934	,	0.923676	,	0.895312	,	-0.161036	,	-0.995895	,\
		0.791047	,	-0.211323	,	-0.302819	,	0.640735	,	-0.317908]	,\
		[0.473908	,	-0.855725	,	-0.0413591	,	-0.508661	,	0.443453	,\
		-0.22712	,	-0.407783	,	0.657463	,	0.0970092	,	-0.579109	,\
		0.518991	,	0.922338	,	0.337886	,	-0.67474	,	-0.725667	,\
		0.570893	,	-0.0798988	,	-0.917135	,	-0.749545	,	-0.982047	,\
		0.405916	,	0.483328	,	0.282047	,	-0.262206	,	0.784123	,\
		-0.795843	,	0.490091	,	0.372046	,	-0.549437	,	0.0964285]	,\
		[-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086	,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808	,\
		-0.792042	,	-0.144765	,	-0.965748	,	0.0133606	,	-0.0260565	,\
		0.101149	,	0.970191	,	0.532821	,	0.814769	,	-0.0687269	,\
		-0.734976	,	-0.342188	,	-0.315861	,	-0.912834	,	0.24499]	,\
		[-0.228237	,	-0.578066	,	0.307023	,	0.606123	,	0.959635	,\
		0.12919	,	0.721925	,	0.766492	,	0.470845	,	-0.0976466	,\
		-0.240557	,	0.66842	,	0.855535	,	-0.451536	,	0.264961	,\
		0.15731	,	0.281697	,	-0.922955	,	-0.780824	,	0.449716	,\
		0.087545	,	-0.0917108	,	-0.62542	,	-0.110256	,	0.0417576	,\
		-0.0408415	,	0.176461	,	0.740113	,	0.470737	,	-0.914927]	,\
		[0.264319	,	-0.73174	,	0.731548	,	-0.489341	,	0.678946	,\
		0.563174	,	-0.814853	,	-0.949609	,	-0.526794	,	-0.801902	,\
		0.691699	,	-0.660499	,	0.389985	,	0.599856	,	-0.711442	,\
		0.128488	,	0.544144	,	-0.495222	,	0.965229	,	-0.79314	,\
		0.0884339	,	-0.222144	,	0.499412	,	-0.565198	,	0.64824	,\
		0.150039	,	-0.0454542	,	0.604861	,	-0.598288	,	-0.500696]	,\
		[0.44187	,	0.940309	,	-0.240334	,	-0.0276121	,	0.74383	,\
		-0.802627	,	0.378082	,	-0.112673	,	-0.47926	,	-0.3355	,\
		-0.0955366	,	0.0187107	,	0.913887	,	0.123076	,	0.550338	,\
		-0.763681	,	0.0781023	,	0.951666	,	0.734031	,	0.826912	,\
		0.0853141	,	-0.583368	,	-0.157612	,	0.234119	,	0.875043	,\
		-0.116586	,	-0.896382	,	-0.817317	,	-0.948837	,	-0.597427]	,\
		[0.842925	,	-0.144503	,	0.936699	,	-0.443935	,	-0.182996	,\
		-0.868866	,	-0.504041	,	0.926483	,	0.169978	,	-0.00842563	,\
		0.916684	,	0.39175	,	0.759081	,	0.496979	,	-0.200691	,\
		0.0626998	,	-0.977963	,	0.660401	,	0.470569	,	-0.0528868	,\
		0.83125	,	-0.662272	,	0.702768	,	0.875814	,	-0.701221	, \
		0.753047	,	-0.86284	,	-0.589688	,	0.178612	,	-0.720358]
		])

		x0 = self.startp()

		self.name  = 'L3ZDT3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]**2

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += (y[i]**2)

		gx = 1.0 + 9.0/float(self.n-1) * gx

		h = 1.0-np.sqrt(f1/gx) - (f1/gx)*np.sin(10.0*np.pi*f1)
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L3ZDT4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418	,	-0.620254	,	0.843784	,	0.914311	,	-0.788548	,\
		-0.0265389	,	-0.920133	,	0.308861	,	-0.0437502	,	-0.374203	,\
		-0.218632	,	-0.865161	,	-0.715997	,	0.220772	,	0.692356	,\
		0.404396	,	0.449996	,	0.162711	,	0.294454	,	-0.563345	,\
		-0.742377	,	0.426391	,	0.408202	,	0.633885	,	-0.0351053	,\
		0.712758	,	-0.191812	,	-0.390938	,	0.952828	,	0.921519]	,\
		[0.0969326	,	0.089775	,	-0.241157	,	0.0835558	,	-0.420236	,\
		0.00225943	,	0.0101814	,	0.441456	,	0.0633629	,	0.406631	,\
		0.507807	,	0.804148	,	0.963269	,	0.357128	,	-0.832565	,\
		-0.61366	,	-0.204783	,	-0.842476	,	-0.249524	,	-0.0985226	,\
		0.125225	,	0.487649	,	0.147046	,	0.679639	,	0.593707	,\
		0.24476	,	0.941339	,	-0.613783	,	0.402772	,	0.300775]	,\
		[-0.591523	,	-0.606614	,	-0.181873	,	0.692975	,	0.50208	,\
		0.0271269	,	0.804879	,	-0.402973	,	0.800373	,	0.760082	,\
		-0.753397	,	0.617418	,	0.689874	,	0.983384	,	0.668786	,\
		-0.798697	,	-0.244945	,	-0.942649	,	0.402856	,	-0.494672	,\
		-0.545421	,	-0.500243	,	0.154371	,	0.170017	,	-0.259108	,\
		0.956914	,	-0.0620912	,	0.634479	,	0.928617	,	0.464664]	,\
		[0.249008	,	0.370711	,	-0.633174	,	-0.0121906	,	0.42006	,\
		-0.630329	,	-0.763376	,	0.62538	,	0.818945	,	0.891598	,\
		-0.699445	,	0.237731	,	-0.324597	,	-0.800406	,	-0.42585	,\
		0.400334	,	-0.367816	,	0.198455	,	-0.983183	,	0.278976	,\
		0.57488	,	-0.361951	,	-0.0739728	,	0.91438	,	-0.391653	,\
		0.430805	,	0.706102	,	0.423887	,	0.296828	,	-0.265607]	,\
		[0.975863	,	-0.971371	,	-0.124115	,	0.4339	,	-0.254671	,\
		0.803564	,	0.960386	,	-0.0323329	,	0.638181	,	-0.895684	,\
		-0.530324	,	0.282745	,	0.0255867	,	0.287686	,	0.410417	,\
		0.0417966	,	-0.687391	,	0.438773	,	0.287357	,	0.316636	,\
		-0.00138645	,	0.931065	,	-0.748519	,	0.304188	,	-0.266153	,\
		0.553793	,	0.471795	,	0.769147	,	0.059668	,	-0.841617]	,\
		[0.428212	,	0.103064	,	-0.47373	,	-0.300792	,	-0.185507	,\
		0.207359	,	-0.219433	,	0.914104	,	0.184408	,	0.520599	,\
		0.646453	,	-0.401724	,	0.615443	,	-0.0601957	,	-0.748176	,\
		-0.114993	,	0.549589	,	-0.775141	,	0.677726	,	0.610715	,\
		-0.723444	,	-0.577654	,	0.0276004	,	0.0712472	,	-0.622791	,\
		0.923094	,	0.93011	,	-0.945394	,	-0.0934027	,	0.964123]	,\
		[-0.686633	,	-0.711276	,	-0.00325377	,	0.435196	,	-0.710002	,\
		-0.0100638	,	-0.177972	,	-0.491075	,	0.537035	,	-0.924987	,\
		-0.312441	,	0.327779	,	0.184745	,	0.246139	,	-0.936814	,\
		0.0671501	,	-0.527707	,	-0.509489	,	-0.883254	,	0.14851	,\
		-0.311828	,	-0.797099	,	-0.35815	,	0.95808	,	0.907244	,\
		-0.820314	,	-0.894233	,	-0.405896	,	0.0735439	,	0.486645]	,\
		[-0.536704	,	0.359652	,	0.839082	,	0.56817	,	-0.0776788	,\
		-0.878364	,	0.176801	,	-0.548932	,	-0.225601	,	-0.164912	,\
		0.0304653	,	-0.625221	,	-0.13318	,	0.827343	,	-0.101358	,\
		0.439941	,	-0.88216	,	0.170196	,	0.650734	,	-0.0982391	,\
		-0.868862	,	-0.50731	,	-0.848317	,	0.835712	,	0.616391	,\
		0.377022	,	0.63047	,	-0.198619	,	-0.576153	,	0.565373]	,\
		[0.169373	,	-0.975542	,	-0.0297852	,	0.80481	,	0.638317	,\
		0.680494	,	0.471868	,	-0.769787	,	-0.878099	,	-0.973724	,\
		-0.710739	,	-0.144068	,	-0.828545	,	-0.800912	,	0.184654	,\
		0.714817	,	0.307911	,	0.812861	,	-0.403497	,	-0.784382	,\
		0.0193887	,	0.412634	,	-0.169813	,	0.471794	,	0.660792	,\
		0.338806	,	-0.15829	,	0.642516	,	0.355126	,	0.174447]	,\
		[0.298068	,	-0.349803	,	-0.73185	,	0.488106	,	-0.0495073	,\
		-0.360502	,	0.0646258	,	-0.202449	,	-0.717228	,	0.970489	,\
		-0.766576	,	-0.536566	,	-0.628288	,	0.69665	,	0.820713	,\
		-0.262311	,	-0.0755541	,	-0.442313	,	0.621378	,	0.670105	,\
		0.672524	,	-0.105179	,	-0.874749	,	-0.154355	,	-0.774656	,\
		-0.191179	,	-0.972471	,	-0.825361	,	0.779826	,	-0.917201]	,\
		[0.330423	,	0.151614	,	0.884043	,	-0.272951	,	-0.993822	,\
		-0.88565	,	-0.375906	,	-0.708948	,	-0.37902	,	0.576578	,\
		-0.207987	,	-0.865931	,	0.613732	,	-0.525712	,	-0.995728	,\
		0.0850755	,	0.0419388	,	-0.323614	,	-0.973719	,	-0.680238	,\
		0.155451	,	0.442717	,	-0.792786	,	0.925785	,	0.670266	,\
		-0.795609	,	-0.289563	,	0.614236	,	-0.670585	,	0.466877]	,\
		[0.00283691	,	-0.168757	,	-0.134045	,	-0.655235	,	0.172361	,\
		-0.699424	,	0.742285	,	0.0181443	,	0.718971	,	-0.0308272	,\
		-0.931734	,	-0.0327827	,	0.319293	,	0.044473	,	-0.641645	,\
		-0.906465	,	0.496238	,	-0.853211	,	-0.779234	,	-0.979515	,\
		0.772426	,	0.720574	,	-0.873217	,	0.371431	,	-0.826029	,\
		-0.394355	,	0.125097	,	-0.316386	,	-0.701215	,	-0.845742]	,\
		[-0.00332785	,	0.459538	,	-0.518313	,	-0.270738	,	-0.629958	,\
		-0.208143	,	0.7768	,	-0.542743	,	-0.156021	,	0.671736	,\
		-0.999522	,	-0.0525574	,	-0.458319	,	0.587409	,	-0.334639	,\
		-0.468732	,	0.342133	,	-0.838071	,	-0.832362	,	0.658177	,\
		-0.442608	,	-0.158	,	0.313451	,	0.703748	,	-0.755984	,\
		-0.524245	,	-0.187299	,	-0.614524	,	0.429316	,	-0.491171]	,\
		[-0.670967	,	0.935792	,	-0.35605	,	0.175773	,	0.878601	,\
		0.354362	,	-0.1792	,	-0.225034	,	-0.44548	,	0.598865	,\
		-0.63675	,	-0.16696	,	0.240427	,	-0.513443	,	0.812664	,\
		-0.161823	,	-0.120835	,	0.323172	,	0.583739	,	0.732924	,\
		-0.350906	,	-0.612644	,	0.347876	,	0.112573	,	-0.501126	,\
		-0.975015	,	0.869905	,	-0.145285	,	-0.484002	,	-0.475966]	,\
		[0.253969	,	0.168116	,	0.148772	,	0.889593	,	-0.512213	,\
		0.404608	,	-0.0861868	,	-0.879417	,	-0.866462	,	-0.938336	,\
		-0.506162	,	-0.404114	,	0.640099	,	-0.956123	,	-0.576586	,\
		0.060982	,	0.944162	,	0.643442	,	-0.750684	,	-0.639973	,\
		-0.69654	,	0.433098	,	0.615897	,	-0.387919	,	-0.429779	,\
		0.43272	,	0.10301	,	0.358771	,	0.793448	,	-0.0379954]	,\
		[0.511197	,	-0.0997948	,	-0.659756	,	0.575496	,	0.675617	,\
		0.0194674	,	-0.470262	,	0.572576	,	0.351245	,	-0.480477	,\
		0.389633	,	-0.064173	,	0.662131	,	-0.707048	,	-0.340423	,\
		-0.270873	,	-0.209617	,	0.968436	,	0.908798	,	0.975851	,\
		-0.865566	,	-0.638281	,	0.333094	,	0.477628	,	0.47261	,\
		0.144597	,	-0.206416	,	0.6937	,	-0.967958	,	-0.0951247]	,\
		[0.998291	,	0.376291	,	-0.962215	,	-0.363174	,	-0.88777	,\
		0.086931	,	0.524476	,	0.956457	,	0.143024	,	0.616481	,\
		0.596118	,	-0.293934	,	-0.63373	,	0.409658	,	0.759892	,\
		0.827175	,	0.228969	,	-0.402829	,	-0.970118	,	0.762559	,\
		0.942716	,	0.70609	,	-0.658158	,	-0.782185	,	-0.806743	,\
		0.2065	,	-0.413743	,	0.406725	,	-0.423813	,	-0.941255]	,\
		[-0.755084	,	-0.721573	,	0.431107	,	-0.221877	,	0.32543	,\
		0.878648	,	-0.419588	,	-0.0752896	,	0.0299447	,	-0.494459	,\
		0.0759643	,	0.0255827	,	0.128944	,	0.17317	,	-0.284309	,\
		-0.565361	,	0.149473	,	0.69331	,	-0.491848	,	0.74916	,\
		-0.249443	,	0.491564	,	0.985068	,	0.678644	,	0.808324	,\
		0.399495	,	-0.333898	,	-0.646636	,	-0.0189709	,	-0.339605]	,\
		[-0.275168	,	-0.932517	,	-0.372497	,	-0.0732907	,	-0.185493	,\
		0.544005	,	-0.478962	,	0.327193	,	-0.525784	,	0.903179	,\
		0.744943	,	0.970612	,	0.00172899	,	-0.726378	,	-0.0985012	,\
		-0.220603	,	-0.594121	,	0.935093	,	-0.216736	,	0.659318	,\
		0.456761	,	-0.109004	,	0.289352	,	-0.566504	,	0.585042	,\
		-0.67704	,	0.996452	,	-0.0685748	,	-0.851985	,	0.416498]	,\
		[-0.165437	,	0.666773	,	-0.976304	,	-0.170024	,	0.905794	,\
		-0.799685	,	0.213464	,	-0.932344	,	-0.668236	,	0.751366	,\
		0.435502	,	-0.470676	,	-0.367062	,	-0.831765	,	-0.294942	,\
		0.217424	,	0.592823	,	0.929094	,	-0.239135	,	-0.41628	,\
		0.650202	,	0.122306	,	-0.237727	,	0.626817	,	-0.227929	,\
		-0.870112	,	0.600442	,	-0.990603	,	0.549151	,	0.512146]	,\
		[0.180332	,	-0.593814	,	-0.492722	,	0.0646786	,	-0.666503	,\
		0.238261	,	-0.1596	,	-0.827302	,	0.669248	,	0.494475	,\
		0.60624	,	0.0951648	,	-0.160446	,	-0.394585	,	-0.167581	,\
		-0.994918	,	-0.0621977	,	0.628171	,	0.761228	,	0.34372	,\
		0.23151	,	0.82132	,	-0.589803	,	0.796275	,	0.57713	,\
		-0.942473	,	-0.610767	,	-0.655472	,	-0.0960986	,	-0.302779]	,\
		[-0.519929	,	-0.560554	,	-0.984415	,	0.601529	,	-0.984103	,\
		0.217909	,	-0.128427	,	-0.262427	,	-0.938208	,	-0.52479	,\
		-0.257078	,	0.939616	,	-0.227661	,	0.115754	,	0.10964	,\
		0.506495	,	0.460303	,	0.897304	,	0.686003	,	0.739986	,\
		-0.627986	,	-0.405551	,	-0.258495	,	-0.796524	,	0.222498	,\
		-0.558804	,	0.312326	,	0.345314	,	0.319143	,	-0.644653]	,\
		[0.163743	,	0.0759916	,	0.695064	,	-0.656856	,	0.074163	,\
		-0.72415	,	0.35978	,	-0.32646	,	-0.96605	,	0.0127605	,\
		0.287161	,	-0.550725	,	-0.433083	,	-0.242821	,	0.878879	,\
		0.526025	,	-0.155339	,	0.0998096	,	0.468761	,	0.324649	,\
		0.81975	,	-0.435823	,	-0.839855	,	0.00282368	,	-0.569165	,\
		-0.798791	,	0.0494081	,	0.367012	,	0.852545	,	0.43557]	,\
		[-0.357004	,	0.314786	,	-0.229239	,	0.530256	,	-0.51327	,\
		-0.899248	,	0.156514	,	0.154329	,	0.499808	,	-0.836327	,\
		0.224232	,	0.16495	,	0.560077	,	-0.813112	,	0.112894	,\
		-0.750417	,	-0.284773	,	-0.271496	,	0.491731	,	-0.712174	,\
		0.584934	,	0.923676	,	0.895312	,	-0.161036	,	-0.995895	,\
		0.791047	,	-0.211323	,	-0.302819	,	0.640735	,	-0.317908]	,\
		[0.473908	,	-0.855725	,	-0.0413591	,	-0.508661	,	0.443453	,\
		-0.22712	,	-0.407783	,	0.657463	,	0.0970092	,	-0.579109	,\
		0.518991	,	0.922338	,	0.337886	,	-0.67474	,	-0.725667	,\
		0.570893	,	-0.0798988	,	-0.917135	,	-0.749545	,	-0.982047	,\
		0.405916	,	0.483328	,	0.282047	,	-0.262206	,	0.784123	,\
		-0.795843	,	0.490091	,	0.372046	,	-0.549437	,	0.0964285]	,\
		[-0.945716	,	-0.334582	,	0.611894	,	0.281032	,	0.508749	,\
		0.691715	,	-0.198585	,	0.0492812	,	0.959669	,	0.884086	,\
		0.0679849	,	0.449799	,	0.733505	,	-0.00918638	,	0.00446808	,\
		-0.792042	,	-0.144765	,	-0.965748	,	0.0133606	,	-0.0260565	,\
		0.101149	,	0.970191	,	0.532821	,	0.814769	,	-0.0687269	,\
		-0.734976	,	-0.342188	,	-0.315861	,	-0.912834	,	0.24499]	,\
		[-0.228237	,	-0.578066	,	0.307023	,	0.606123	,	0.959635	,\
		0.12919	,	0.721925	,	0.766492	,	0.470845	,	-0.0976466	,\
		-0.240557	,	0.66842	,	0.855535	,	-0.451536	,	0.264961	,\
		0.15731	,	0.281697	,	-0.922955	,	-0.780824	,	0.449716	,\
		0.087545	,	-0.0917108	,	-0.62542	,	-0.110256	,	0.0417576	,\
		-0.0408415	,	0.176461	,	0.740113	,	0.470737	,	-0.914927]	,\
		[0.264319	,	-0.73174	,	0.731548	,	-0.489341	,	0.678946	,\
		0.563174	,	-0.814853	,	-0.949609	,	-0.526794	,	-0.801902	,\
		0.691699	,	-0.660499	,	0.389985	,	0.599856	,	-0.711442	,\
		0.128488	,	0.544144	,	-0.495222	,	0.965229	,	-0.79314	,\
		0.0884339	,	-0.222144	,	0.499412	,	-0.565198	,	0.64824	,\
		0.150039	,	-0.0454542	,	0.604861	,	-0.598288	,	-0.500696]	,\
		[0.44187	,	0.940309	,	-0.240334	,	-0.0276121	,	0.74383	,\
		-0.802627	,	0.378082	,	-0.112673	,	-0.47926	,	-0.3355	,\
		-0.0955366	,	0.0187107	,	0.913887	,	0.123076	,	0.550338	,\
		-0.763681	,	0.0781023	,	0.951666	,	0.734031	,	0.826912	,\
		0.0853141	,	-0.583368	,	-0.157612	,	0.234119	,	0.875043	,\
		-0.116586	,	-0.896382	,	-0.817317	,	-0.948837	,	-0.597427]	,\
		[0.842925	,	-0.144503	,	0.936699	,	-0.443935	,	-0.182996	,\
		-0.868866	,	-0.504041	,	0.926483	,	0.169978	,	-0.00842563	,\
		0.916684	,	0.39175	,	0.759081	,	0.496979	,	-0.200691	,\
		0.0626998	,	-0.977963	,	0.660401	,	0.470569	,	-0.0528868	,\
		0.83125	,	-0.662272	,	0.702768	,	0.875814	,	-0.701221	, \
		0.753047	,	-0.86284	,	-0.589688	,	0.178612	,	-0.720358]
		])

		x0 = self.startp()

		self.name  = 'L3ZDT4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]**2

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += y[i]**2 -10.0*np.cos(4.0*np.pi*y[i])

		gx = 1.0 + 10.0*float(self.n-1) + gx

		h = 1.0-np.sqrt(f1/gx)
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_L3ZDT6():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.MAT   = np.array([ \
		[0.218418,    -0.620254  ,   0.843784  ,  0.914311   , -0.788548  ,\
	0.428212,	0.103064  ,  -0.47373   ,  -0.300792   ,  -0.185507] ,\
	[0.330423,     0.151614  ,   0.884043  , -0.272951   , -0.993822  ,\
	0.511197,	-0.0997948,   -0.659756 ,    0.575496  ,    0.675617],\
	[0.180332 ,   -0.593814  ,  -0.492722  ,  0.0646786  , -0.666503  ,\
	-0.945716,	-0.334582 ,    0.611894 ,    0.281032  ,    0.508749],\
	[-0.0265389,   -0.920133 ,    0.308861 ,  -0.0437502 ,  -0.374203 ,\
	0.207359,	-0.219433 ,    0.914104 ,    0.184408  ,    0.520599],\
	[-0.88565  ,   -0.375906 ,   -0.708948 ,  -0.37902   ,   0.576578 ,\
	0.0194674,	-0.470262 ,    0.572576 ,    0.351245  ,   -0.480477],\
	[0.238261  ,  -0.1596    ,  -0.827302  ,  0.669248   ,  0.494475  ,\
	0.691715,	-0.198585 ,    0.0492812,    0.959669  ,    0.884086],\
	[-0.218632 ,   -0.865161 ,   -0.715997 ,   0.220772  ,   0.692356 ,\
	0.646453,	-0.401724 ,    0.615443 ,   -0.0601957 ,   -0.748176],\
	[-0.207987 ,   -0.865931 ,    0.613732 ,  -0.525712  ,  -0.995728 ,\
	0.389633,	-0.064173 ,    0.662131 ,   -0.707048  ,   -0.340423],\
	[0.60624   ,   0.0951648 ,  -0.160446  , -0.394585   , -0.167581  ,\
	0.0679849,	0.449799  ,   0.733505  ,  -0.00918638 ,   0.00446808],\
	[0.404396  ,   0.449996  ,   0.162711  ,  0.294454   , -0.563345  , \
	-0.114993,	0.549589  ,  -0.775141  ,   0.677726   ,   0.610715] \
		])

		x0 = self.startp()

		self.name  = 'L3ZDT6'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		MAT = self.MAT
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = 0.0
			for j in range(self.n): #j = 1,n
				y[i] += MAT[j,i]*x[j]**2

		f1 = y[0]**2

		gx = 0.0
		for i in range(1,self.n): #i = 2,n
			gx += y[i]**2

		gx = 1.0 + 9.0*(gx/float(self.n-1))**0.25

		h = 1.0-(f1/gx)**2
		f2 = gx*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_LE1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -5.0*np.ones(self.n)
		self.ub    = 10.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'LE1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = (x[0]**2+x[1]**2)**0.125
		f2 = ((x[0]-0.5)**2+(x[1]-0.5)**2)**0.25

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_LOVISON1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = 3.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'LOVISON1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 1.05*x[0]**2+0.98*x[1]**2
		f2 = 0.99*(x[0]-3.0)**2+1.03*(x[1]-2.5)**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_LOVISON2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -0.5*np.ones(self.n)
		self.ub    =  0.5*np.ones(self.n); self.ub[0] = 0.0

		x0 = self.startp()

		self.name  = 'LOVISON2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[1]
		f2 = -(x[1]-x[0]**3)/(x[0]+1.0)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_LOVISON3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n); self.lb[1] = -4.0
		self.ub    = 6.0*np.ones(self.n); self.ub[1] = 4.0

		x0 = self.startp()

		self.name  = 'LOVISON3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]**2+x[1]**2
		f2 = (x[0]-6.0)**2+(x[1]+0.3)**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_LOVISON4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n); self.lb[1] = -1.0
		self.ub    = 6.0*np.ones(self.n); self.ub[1] = 1.0

		x0 = self.startp()

		self.name  = 'LOVISON4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]**2+x[1]**2+4.0*(np.exp(-(x[0]+2.0)**2-x[1]**2)-np.exp(-(x[0]-2.0)**2-x[1]**2))
		f2 = (x[0]-6.0)**2+(x[1]+0.5)**2.0

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_LOVISON5():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 3
		self.n     = n
		self.q     = 3

		self.lb    = -1.0*np.ones(self.n)
		self.ub    =  4.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'LOVISON5'
		self.C = np.array([ \
		[0.218418, -0.620254, 0.843784], \
		[0.914311, -0.788548, 0.428212], \
		[0.103064, -0.47373 , -0.300792] \
		])
		self.alpha = np.array([ \
		[0.407247, 0.665212, 0.575807], \
		[0.942022, 0.363525, 0.00308876], \
		[0.755598, 0.450103, 0.170122] \
		])
		self.beta = np.array([0.575496, 0.675617, 0.180332])
		self.gamma = np.array([-0.593814, -0.492722, 0.0646786])

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f = np.zeros(self.q)
		for j in range(self.q):
			for i in range(self.n):
				f[j] += -self.alpha[i,j]*(x[i]-self.C[j,i])**2

		f[0] = -f[0]
		f[1] = -f[1]-self.beta[1]*np.sin(np.pi*(x[0]+x[1])/self.gamma[1])
		f[2] = -f[2]-self.beta[2]*np.cos(np.pi*(x[0]-x[1])/self.gamma[2])

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_LOVISON6():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 3
		self.n     = n
		self.q     = 3

		self.lb    = -1.0*np.ones(self.n)
		self.ub    =  4.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'LOVISON6'
		self.C = np.array([ \
		[0.218418, -0.620254, 0.843784], \
		[0.914311, -0.788548, 0.428212], \
		[0.103064, -0.47373, -0.300792], \
		[-0.185507, 0.330423, 0.151614]  \
		])
		self.alpha = np.array([ \
		[0.942022, 0.363525, 0.00308876, 0.755598], \
		[0.450103, 0.170122, 0.787748, 0.837808],   \
		[0.590166, 0.203093, 0.253639, 0.532339]    \
		])
		self.beta = np.array([-0.666503, -0.945716, -0.334582, 0.611894])
		self.gamma = np.array([0.281032, 0.508749, -0.0265389, -0.920133])

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f = np.zeros(4)
		for j in range(4):
			for i in range(self.n):
				f[j] += -self.alpha[i,j]*(x[i]-self.C[j,i])**2

		u = np.zeros(self.q)
		u[0] = -f[0]-self.beta[0]*np.exp(f[3]/self.gamma[0])
		u[1] = -f[1]-self.beta[1]*np.sin(np.pi*(x[0]+x[1])/self.gamma[1])
		u[2] = -f[2]-self.beta[2]*np.cos(np.pi*(x[0]-x[1])/self.gamma[2])

		return u

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_LRS1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -50.0*np.ones(self.n)
		self.ub    =  50.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'LRS1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]**2+x[1]**2
		f2 = (x[0]+2.0)**2+x[1]**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MHHM1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 1
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MHHM1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = (x[0]-0.8)**2
		f2 = (x[0]-0.85)**2
		f3 = (x[0]-0.9)**2

		return np.array([f1, f2, f3])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MHHM2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MHHM2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = (x[0]-0.8)**2+(x[1]-0.6)**2
		f2 = (x[0]-0.85)**2+(x[1]-0.7)**2
		f3 = (x[0]-0.9)**2+(x[1]-0.6)**2

		return np.array([f1, f2, f3])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MLF1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 1
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = 20.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MLF1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = (1.0+x[0]/20.0)*np.sin(x[0])
		f2 = (1.0+x[0]/20.0)*np.cos(x[0])

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MLF2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -2.0*np.ones(self.n)
		self.ub    =  2.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MLF2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		y = 2.*x

		f1 = 5.0-((x[0]**2+x[1]-11.0)**2+(x[0]+x[1]**2-7.0)**2)/200.0
		f1 = -f1

		f2 = 5.0-((y[0]**2+y[1]-11.0)**2+(y[0]+y[1]**2-7.0)**2)/200.0
		f2 = -f2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MOP1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 1
		self.n     = n
		self.q     = 2

		self.lb    = -1.e+5*np.ones(self.n)
		self.ub    =  1.e+5*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MOP1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]**2
		f2 = (x[0]-2.0)**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MOP2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 4
		self.n     = n
		self.q     = 2

		self.lb    = -4.0*np.ones(self.n)
		self.ub    =  4.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MOP2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.0
		for i in range(self.n): #i = 1,n
			f1 += (x[i]-1.0/np.sqrt(float(self.n)))**2

		f1 = 1.0 - np.exp(-f1)

		f2 = 0.0
		for i in range(self.n): #i = 1,n
			f2 += (x[i]+1.0/np.sqrt(float(self.n)))**2

		f2 = 1.0 - np.exp(-f2)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MOP3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -np.pi*np.ones(self.n)
		self.ub    =  np.pi*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MOP3'

		self.A1 = 0.50*np.sin(1.0)-2.0*np.cos(1.0)+np.sin(2.0)-1.5*np.cos(2.0)
		self.A2 = 1.50*np.sin(1.0)-np.cos(1.0)+2.0*np.sin(2.0)-0.5*np.cos(2.0)

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		B1 = 0.5*np.sin(x[0])-2.0*np.cos(x[0])+np.sin(x[1])-1.5*np.cos(x[1])
		B2 = 1.5*np.sin(x[0])-np.cos(x[0])+2.0*np.sin(x[1])-0.5*np.cos(x[1])

		f1 = 1.0 + (self.A1-B1)**2 + (self.A2-B2)**2
		f2 = (x[0]+3.0)**2 + (x[1]+1.0)**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MOP4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 3
		self.n     = n
		self.q     = 2

		self.lb    = -5.0*np.ones(self.n)
		self.ub    =  5.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MOP4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.0
		for i in range(self.n-1): #i = 1,n-1
			f1 += (-10.0*np.exp(-0.2*np.sqrt(x[i]**2+x[i+1]**2)))

		f2 = 0.0
		for i in range(self.n): #i = 1,n
			f2 += (np.abs(x[i])**0.8+5.0*np.sin(x[i]**3))

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MOP5():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 3

		self.lb    = -30.0*np.ones(self.n)
		self.ub    =  30.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MOP5'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.5*(x[0]**2+x[1]**2)+np.sin(x[0]**2+x[1]**2)
		f2 = (3.0*x[0]-2.0*x[1]+4.0)**2/8.0+(x[0]-x[1]+1.0)**2/27.0+15.0
		f3 = 1.0/(x[0]**2+x[1]**2+1.0)-1.1*np.exp(-x[0]**2-x[1]**2)

		return np.array([f1, f2, f3])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MOP6():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MOP6'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]
		f2 = (1.0+10.0*x[1])*(1.0-(x[0]/(1.0+10.0*x[1]))**2-x[0]/(1.0+10.0*x[1])*np.sin(8*np.pi*x[0]))

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_MOP7():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 3

		self.lb    = -400.0*np.ones(self.n)
		self.ub    =  400.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'MOP7'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = (x[0]-2.0)**2/2.0+(x[1]+1.0)**2/13.0+3.0
		f2 = (x[0]+x[1]-3.0)**2/36.0+(-x[0]+x[1]+2.0)**2/8.0-17.0
		f3 = (x[0]+2.0*x[1]-1.0)**2/175.0+(-x[0]+2.0*x[1])**2/17.0-13.0

		return np.array([f1, f2, f3])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_OKA1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)

		self.lb[0] = 6.0*np.sin(np.pi/12.0); 	    self.ub[0] = 6.0*np.sin(np.pi/12.0)+2.0*np.pi*np.cos(np.pi/12.0)
		self.lb[1] = -2.0*np.pi*np.sin(np.pi/12.0);	self.ub[1] = 6.0*np.cos(np.pi/12.0)

		x0 = self.startp()

		self.name  = 'OKA1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return (self.lb+self.lb)/2.0

	def functs(self,x):

		y = np.zeros(self.n)
		y[0] = np.cos(np.pi/12.0)*x[0]-np.sin(np.pi/12.0)*x[1]
		y[1] = np.sin(np.pi/12.0)*x[0]+np.cos(np.pi/12.0)*x[1]

		f1 = y[0]
		f2 = np.sqrt(2.0*np.pi)-np.sqrt(np.abs(y[0]))+2.0*np.abs(y[1]-3.0*np.cos(y[0])-3.0)**(1.0/3.0)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_OKA2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 3
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)

		self.lb[0] = -np.pi; self.ub[0] = np.pi
		self.lb[1] = -5.0;	self.ub[1] = 5.0
		self.lb[2] = -5.0;	self.ub[2] = 5.0

		x0 = self.startp()

		self.name  = 'OKA2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]
		f2 = 1.0-(x[0]+np.pi)**2.0/(4.0*np.pi**2)+np.abs(x[1]-5.0*np.cos(x[0]))**(1.0/3.0)+ \
			 np.abs(x[2]-5.0*np.sin(x[0]))**(1.0/3.0)

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_QV1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 2

		self.lb    = -5.12*np.ones(self.n)
		self.ub    =  5.12*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'QV1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 0.0
		for i in range(self.n): #i = 1,n
			f1 += (x[i]**2-10.0*np.cos(2.0*np.pi*x[i])+10.0)
		f1 = (f1/float(self.n))**0.25

		f2 = 0.0
		for i in range(self.n): #i = 1,n
			f2 += ((x[i]-1.5)**2-10.0*np.cos(2.0*np.pi*(x[i]-1.5))+10.0)
		f2 = (f2/float(self.n))**0.25

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_SCH1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 1
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = 5.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'SCH1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		if(x[0] <= 1.0):
			f1 = -x[0]
		elif (x[0] <= 3.0):
			f1 = -2.0 + x[0]
		elif (x[0] <= 4.0):
			f1 = 4.0 - x[0]
		else:
			f1 = -4.0 + x[0]

		f2 = (x[0]-5.0)**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_SK1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 1
		self.n     = n
		self.q     = 2

		self.lb    = -10.0*np.ones(self.n)
		self.ub    =  10.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'SK1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = -x[0]**4-3.0*x[0]**3+10.0*x[0]**2+10.0*x[0]+10.0
		f2 = 0.5*x[0]**4+2.0*x[0]**3+10.0*x[0]**2-10.0*x[0]+5.0

		f1 = -f1
		f2 = -f2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_SK2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 4
		self.n     = n
		self.q     = 2

		self.lb    = -10.0*np.ones(self.n)
		self.ub    =  10.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'SK2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = -(x[0]-2.0)**2-(x[1]+3.0)**2-(x[2]-5.0)**2-(x[3]-4.0)**2+5.0

		a = 0.0
		b = 0.0
		for i in range(self.n): #i = 1,n
			a += np.sin(x[i])
			b += x[i]**2
		f2 = a / (1.0 + b/100.0)

		f1 = -f1
		f2 = -f2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_SP1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -1.0*np.ones(self.n)
		self.ub    =  5.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'SP1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = (x[0]-1.0)**2+(x[0]-x[1])**2
		f2 = (x[1]-3.0)**2+(x[0]-x[1])**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_SSFYY1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -100.0*np.ones(self.n)
		self.ub    =  100.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'SSFYY1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]**2+x[1]**2
		f2 = (x[0]-1.0)**2+(x[1]-2.0)**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_SSFYY2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 1
		self.n     = n
		self.q     = 2

		self.lb    = -100.0*np.ones(self.n)
		self.ub    =  100.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'SSFYY2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 10.0+x[0]**2-10.0*np.cos(x[0]*np.pi/2.0)
		f2 = (x[0]-4.0)**2

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_TKLY1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 4
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n); self.lb[0] = 0.1
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'TKLY1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return (self.lb+self.ub)/2.0

	def functs(self,x):

		f1 = x[0]
		f2 = 1.0
		for i in range(1,self.n): #i = 2,n
			f2 *= (2.0-np.exp(-((x[i]-0.1)/0.004)**2)-0.8*np.exp(-((x[i]-0.9)/0.4)**2))
		f2 = f2 / x[0]

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_VFM1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 3

		self.lb    = -2.0*np.ones(self.n)
		self.ub    =  2.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'VFM1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]**2+(x[1]-1.0)**2
		f2 = x[0]**2+(x[1]+1.0)**2 + 1.0
		f3 = (x[0]-1.0)**2+x[1]**2 + 2.0

		return np.array([f1, f2, f3])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_VU1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -3.0*np.ones(self.n)
		self.ub    =  3.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'VU1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 1.0/(x[0]**2+x[1]**2+1.0)
		f2 = x[0]**2+3.0*x[1]**2+1.0

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_VU2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 2
		self.n     = n
		self.q     = 2

		self.lb    = -3.0*np.ones(self.n)
		self.ub    =  3.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'VU2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]+x[1]+1.0
		f2 = x[0]**2+2.0*x[1]-1.0

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = 2.0*float(i+1)

		x0 = self.startp()

		self.name  = 'WFG1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 0.8
		BB = 0.75
		CC = 0.85
		AAA = 0.02
		alpha = 1.0
		AAAA = 5.0
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		A = np.ones(self.q-1)
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		for i in range(k): #i = 1,k
			t1[i] = y[i]
		for i in range(k,self.n): #i = k+1,n
			t1[i] = np.abs(y[i]-0.35)/np.abs(np.floor(0.35-y[i])+0.35)
		t2 = np.zeros(self.n)
		for i in range(k): #i = 1,k
			t2[i] = t1[i]
		for i in range(k,self.n): #i = k+1,n
			t2[i] = AA + np.minimum(0.0,float(np.floor(t1[i]-BB)))*(AA*(BB-t1[i])) / \
				    BB - np.minimum(0.0,float(np.floor(CC-t1[i])))*(1.0-AA)*(t1[i]-CC)/(1.0-CC)
		t3 = np.zeros(self.n)
		w = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			t3[i] = t2[i]**AAA
			w[i] = 2.0*float(i+1)
		t4 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			gsum1 = 0.0
			gsum2 = 0.0
			for j in range(i*2,(i+1)*2): #j = ((i-1)*k/(M-1)+1), (i*k/(M-1))
				gsum1 += (w[j]*t3[j])
			for j in range(i*2,(i+1)*2): #j = ((i-1)*k/(M-1)+1), (i*k/(M-1))
				gsum2 += w[j]
			t4[i] = gsum1/gsum2
		gsum1 = 0.0
		gsum2 = 0.0
		for j in range(k,self.n): #j = k+1,n
			gsum1 += (w[j]*t3[j])
			gsum2 += w[j]
		t4[self.q-1] = gsum1/gsum2
		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t4[self.q-1],A[i])*(t4[i]-0.5)+0.5
		x[self.q-1] = t4[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*(1.0-np.cos(x[i]*np.pi/2.0))
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*(1.0-np.cos(x[i]*np.pi/2.0))
			h[j] = h[j]*(1.0-np.sin(x[self.q-j]*np.pi/2.0))
		h[self.q-1] = (1.0-x[0]-(np.cos(2.0*AAAA*np.pi*x[0]+np.pi/2.0))/(2.0*AAAA*np.pi))**alpha

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = 2.0*float(i+1)

		x0 = self.startp()

		self.name  = 'WFG2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 2.0
		alpha = 1.0
		beta = 1.0
		AAAA = 5.0
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		A = np.ones(self.q-1)
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		for i in range(k): #i = 1,k
			t1[i] = y[i]
		for i in range(k,self.n): #i = k+1,n
			t1[i] = np.abs(y[i]-0.35)/np.abs(np.floor(0.35-y[i])+0.35)
		t2 = np.zeros(k+2)
		for i in range(k): #i = 1,k
			t2[i] = t1[i]
		for i in range(k,k+2): #i = k+1,k+l/2
			gsum2 = 0.0
			for ii in range(k+2*(i-k),k+2*(i-k)+2): #ii = (k+2*(i-k)-1), (k+2*(i-k))
				gsum1 = 0.0
				for jj in range(int(AA-1)): #jj = 0, AA-2
					gsum1 += np.abs(t1[ii]-t1[(k+2*(i-k))+(((ii+jj-(k+2*(i-k))+1)) % (k+2*(i-k)+1-(k+2*(i-k))+1))])
				gsum2 += t1[ii] + gsum1
			t2[i] = gsum2 / (float(2)/float(AA)*(AA/2.0)*(1.0+2.0*float(AA)-2.0*(AA/2.0)))

		t3 = np.zeros(self.q)
		w = np.ones(self.n)
		for i in range(self.q-1): #i = 1,M-1
			gsum1 = 0.0
			gsum2 = 0.0
			for j in range(2*i,2*i+2): #j = ((i-1)*2+1), (i*2)
				gsum1 += w[j]*t2[j]
				gsum2 += w[j]
			t3[i] = gsum1 / gsum2
		gsum1 = 0.0
		gsum2 = 0.0
		for j in range(k,k+2): #j = k+1,k+l/2
			gsum1 += w[j]*t2[j]
			gsum2 += w[j]
		t3[self.q-1] = gsum1 / gsum2
		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],A[i])*(t3[i]-0.5)+0.5
		x[self.q-1] = t3[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*(1.0-np.cos(x[i]*np.pi/2.0))
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*(1.0-np.cos(x[i]*np.pi/2.0))
			h[j] = h[j]*(1.0-np.sin(x[self.q-j]*np.pi/2.0))
		h[self.q-1] = (1.0-x[0]-(np.cos(2.0*AAAA*np.pi*x[0]+np.pi/2.0))/(2.0*AAAA*np.pi))**alpha

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = 2.0*float(i+1)

		x0 = self.startp()

		self.name  = 'WFG3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 2.0
		alpha = 1.0
		beta = 1.0
		AAAA = 5.0
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		A = np.zeros(self.q-1)
		A[0] = 1.0
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		for i in range(k): #i = 1,k
			t1[i] = y[i]
		for i in range(k,self.n): #i = k+1,n
			t1[i] = np.abs(y[i]-0.35)/np.abs(np.floor(0.35-y[i])+0.35)
		t2 = np.zeros(k+2)
		for i in range(k): #i = 1,k
			t2[i] = t1[i]
		for i in range(k,k+2): #i = k+1,k+l/2
			gsum2 = 0.0
			for ii in range(k+2*(i-k),k+2*(i-k)+2): #ii = (k+2*(i-k)-1), (k+2*(i-k))
				gsum1 = 0.0
				for jj in range(int(AA-1)): #jj = 0, AA-2
					gsum1 += np.abs(t1[ii]-t1[(k+2*(i-k))+(((ii+jj-(k+2*(i-k))+1)) % (k+2*(i-k)+1-(k+2*(i-k))+1))])
				gsum2 += t1[ii] + gsum1
			t2[i] = gsum2 / (float(2)/float(AA)*(AA/2.0)*(1.0+2.0*float(AA)-2.0*(AA/2.0)))

		t3 = np.zeros(self.q)
		w = np.ones(self.n)
		for i in range(self.q-1): #i = 1,M-1
			gsum1 = 0.0
			gsum2 = 0.0
			for j in range(2*i,2*i+2): #j = ((i-1)*2+1), (i*2)
				gsum1 += w[j]*t2[j]
				gsum2 += w[j]
			t3[i] = gsum1 / gsum2
		gsum1 = 0.0
		gsum2 = 0.0
		for j in range(k,k+2): #j = k+1,k+l/2
			gsum1 += w[j]*t2[j]
			gsum2 += w[j]
		t3[self.q-1] = gsum1 / gsum2
		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],A[i])*(t3[i]-0.5)+0.5
		x[self.q-1] = t3[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*x[i]
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*x[i]
			h[j] = h[j]*(1.0-x[self.q-j])
		h[self.q-1] = 1.0-x[0]

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = float(i+1)

		x0 = self.startp()

		self.name  = 'WFG4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 30
		BB = 10
		CC = 0.35
		alpha = 1.0
		beta = 1.0
		AAAA = 5.0
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		A = np.ones(self.q-1)
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		for i in range(self.n): #i = 1,n
			t1[i] = (1.0+np.cos((4.0*AA+2.0)*np.pi*(0.5-np.abs(y[i]-CC)/(2.0*(np.floor(CC-y[i])+CC)))) \
			      + 4.0*BB*(np.abs(y[i]-CC)/(2.0*(np.floor(CC-y[i])+CC)))**2)/(BB+2.0)
		t2 = np.zeros(self.q)
		w = np.ones(self.n)

		for i in range(self.q-1): #i=1,M-1
			gsum1=0.0
			gsum2=0.0
			for ii in range(i*2,i*2+2): #ii=((i-1)*2+1),(i*2)
				gsum1 += (w[ii]*t1[ii])
				gsum2 += w[ii]
			t2[i] = gsum1/gsum2
		gsum1=0.0
		gsum2=0.0
		for ii in range(k,self.n): #ii=k+1,n
			gsum1 += (w[ii]*t1[ii])
			gsum2 += w[ii]
		t2[self.q-1] = gsum1/gsum2

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t2[self.q-1],A[i])*(t2[i]-0.5)+0.5
		x[self.q-1] = t2[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*np.sin(x[i]*np.pi/2.0)
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*np.sin(x[i]*np.pi/2.0)
			h[j] = h[j]*np.cos(x[self.q-j]*np.pi/2.0)
		h[self.q-1] = np.cos(x[0]*np.pi/2.0)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG5():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = 2.0*float(i+1)

		x0 = self.startp()

		self.name  = 'WFG5'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 0.35
		BB = 0.001
		CC = 0.05
		AAA = 0.02
		alpha = 1.0
		AAAA = 5.0

		A = np.ones(self.q-1)
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		for i in range(self.n): #i = 1,n
			t1[i] = 1.0+(np.abs(y[i]-AA)-BB)*((np.floor(y[i]-AA+BB)*(1.0-CC+(AA-BB)/BB))/(AA-BB)+	\
			       (np.floor(AA+BB-y[i])*(1.0-CC+(1.0-AA-BB)/BB))/(1.0-AA-BB)+1.0/BB)

		t2 = np.zeros(self.q)
		w = np.ones(self.n)

		for i in range(self.q-1): #i=1,M-1
			gsum1=0.0
			gsum2=0.0
			for ii in range(i*2,i*2+2): #ii=((i-1)*2+1),(i*2)
				gsum1 += (w[ii]*t1[ii])
				gsum2 += w[ii]
			t2[i] = gsum1/gsum2
		gsum1=0.0
		gsum2=0.0
		for ii in range(k,self.n): #ii=k+1,n
			gsum1 += (w[ii]*t1[ii])
			gsum2 += w[ii]
		t2[self.q-1] = gsum1/gsum2

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t2[self.q-1],A[i])*(t2[i]-0.5)+0.5
		x[self.q-1] = t2[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*np.sin(x[i]*np.pi/2.0)
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*np.sin(x[i]*np.pi/2.0)
			h[j] = h[j]*np.cos(x[self.q-j]*np.pi/2.0)
		h[self.q-1] = np.cos(x[0]*np.pi/2.0)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG6():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = 2.0*float(i+1)

		x0 = self.startp()

		self.name  = 'WFG6'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 0.35
		BB = 0.001
		CC = 0.05
		AAA = 0.02
		alpha = 1.0
		AAAA = 5.0

		A = np.ones(self.q-1)
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		for i in range(k): #i = 1,k
			t1[i] = y[i]
		for i in range(k,self.n): #i = k+1,n
			t1[i] = np.abs(y[i]-0.35)/np.abs(np.floor(0.35-y[i])+0.35)

		t2 = np.zeros(self.q)
		w = np.ones(self.n)

		for i in range(self.q-1): #i = 1,M-1
			gsum2 = 0.0
			for ii in range(i*2,i*2+2): #ii = ((i-1)*2+1), (i*2)
				gsum1 = 0.0
				for jj in range(1): #jj = 0, (k/(M-1)-2)
					gsum1 += np.abs(t1[ii]-t1[(i*2)+((ii+jj-(i*2)+1) % 2)])
				gsum2 += t1[ii] + gsum1
			#t2[i] = gsum2 / (((float(i+1)*float(k)/2.0)-(float(i+1)*float(k)/2.0+1.0)+1.0)/  \
			#		 (float(k)/2.0)*np.ceil(float(k)/2.0/2.0)*(1.0+2.0*float(k)/  \
			#			2.0-2*np.ceil(float(k)/2.0/2.0)))
			t2[i] = gsum2 / (((float(i+1)*2.0)-(float(i)*2.0+1.0)+1.0)/(2.0)*6.0)

		t2[self.q-1] = 0.0
		for ii in range(k,self.n): #ii = k+1,n
			gsum1 = 0.0
			for jj in range(3): #jj = 0, l-2
				gsum1 += np.abs(t1[ii]-t1[k+((ii+jj-(k+1)+1) % (self.n-k))])
			t2[self.q-1] += t1[ii] + gsum1
		t2[self.q-1] = t2[self.q-1] / ((float(self.n-k)/4)*2.0*(5))

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t2[self.q-1],A[i])*(t2[i]-0.5)+0.5
		x[self.q-1] = t2[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*np.sin(x[i]*np.pi/2.0)
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*np.sin(x[i]*np.pi/2.0)
			h[j] = h[j]*np.cos(x[self.q-j]*np.pi/2.0)
		h[self.q-1] = np.cos(x[0]*np.pi/2.0)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG7():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = 2.0*float(i+1)

		x0 = self.startp()

		self.name  = 'WFG7'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 0.98/49.98
		BB = 0.02
		CC = 50.0
		AAA = 0.02
		alpha = 1.0
		AAAA = 5.0

		w = np.ones(self.n)
		A = np.ones(self.q-1)
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		rsum = np.zeros(k)
		for i in range(k): #i = 1,k
			gsum1 = 0.0
			gsum2 = 0.0
			for j in range(i+1,self.n): #j = i+1,n
				gsum1 += w[j]*y[j]
				gsum2 += w[j]
			rsum[i] = gsum1 / gsum2
		for i in range(k): #i = 1,k
			t1[i] = y[i]**(BB+(CC-BB)*(AA-(1.0-2.0*rsum[i])*np.abs(np.floor(0.5-rsum[i])+AA)))
		for i in range(k,self.n): #i = k+1,n
			t1[i] = y[i]

		t2 = np.zeros(self.n)

		for i in range(k): #i = 1,k
			t2[i] = t1[i]
		for i in range(k,self.n): #i = k+1,n
			t2[i] = np.abs(t1[i]-0.35)/np.abs(np.floor(0.35-t1[i])+0.35)

		t3 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			gsum1 = 0.0
			gsum2 = 0.0
			for j in range(2*i,2*i+2): #j = ((i-1)*2+1),(i*2)
				gsum1 += w[j]*t2[j]
				gsum2 += w[j]
			t3[i] = gsum1 / gsum2
		gsum1 = 0.0
		gsum2 = 0.0
		for j in range(k,self.n): #j = k+1,n
			gsum1 += w[j]*t2[j]
			gsum2 += w[j]
		t3[self.q-1] = gsum1 / gsum2

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],A[i])*(t2[i]-0.5)+0.5
		x[self.q-1] = t3[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*np.sin(x[i]*np.pi/2.0)
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*np.sin(x[i]*np.pi/2.0)
			h[j] = h[j]*np.cos(x[self.q-j]*np.pi/2.0)
		h[self.q-1] = np.cos(x[0]*np.pi/2.0)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG8():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = 2.0*float(i+1)

		x0 = self.startp()

		self.name  = 'WFG8'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 0.98/49.98
		BB = 0.02
		CC = 50.0
		AAA = 0.02
		alpha = 1.0
		AAAA = 5.0

		w = np.ones(self.n)
		A = np.ones(self.q-1)
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		rsum = np.zeros(k)
		for i in range(k,self.n): #i = k+1,n
			gsum1 = 0.0
			gsum2 = 0.0
			for j in range(i): #j = 1,i-1
				gsum1 += w[j]*y[j]
				gsum2 += w[j]
			rsum[i-k] = gsum1 / gsum2
		for i in range(k): #i = 1,k
			t1[i] = y[i]
		for i in range(k,self.n): #i = k+1,n
			t1[i] = y[i]**(BB+(CC-BB)*(AA-(1.0-2.0*rsum[i-k])*np.abs(np.floor(0.5-rsum[i-k])+AA)))
		t2 = np.zeros(self.n)

		for i in range(k): #i = 1,k
			t2[i] = t1[i]
		for i in range(k,self.n): #i = k+1,n
			t2[i] = np.abs(t1[i]-0.35)/np.abs(np.floor(0.35-t1[i])+0.35)

		t3 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			gsum1 = 0.0
			gsum2 = 0.0
			for j in range(2*i,2*i+2): #j = ((i-1)*2+1),(i*2)
				gsum1 += w[j]*t2[j]
				gsum2 += w[j]
			t3[i] = gsum1 / gsum2
		gsum1 = 0.0
		gsum2 = 0.0
		for j in range(k,self.n): #j = k+1,n
			gsum1 += w[j]*t2[j]
			gsum2 += w[j]
		t3[self.q-1] = gsum1 / gsum2

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],A[i])*(t2[i]-0.5)+0.5
		x[self.q-1] = t3[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*np.sin(x[i]*np.pi/2.0)
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*np.sin(x[i]*np.pi/2.0)
			h[j] = h[j]*np.cos(x[self.q-j]*np.pi/2.0)
		h[self.q-1] = np.cos(x[0]*np.pi/2.0)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_WFG9():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 8
		self.n     = n
		self.q     = 3

		self.lb    = np.zeros(self.n)
		self.ub    = np.zeros(self.n)
		for i in range(self.n):
			self.ub[i] = 2.0*float(i+1)

		x0 = self.startp()

		self.name  = 'WFG9'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,z):

		AA = 0.98/49.98
		BB = 0.02
		CC = 50.0
		AAA = 0.35
		BBB = 0.001
		CCC = 0.05
		AAAA = 30.0
		BBBB = 95.0
		CCCC = 0.35

		w = np.ones(self.n)
		A = np.ones(self.q-1)
		S = 2.0*np.ones(self.n)
		for i in range(self.q): #i = 1,M
			S[i] *= float(i+1)
		zmax = 2.0*np.ones(self.n)
		for i in range(self.n): #i = 1,n
			zmax[i] *= float(i+1)
		y = np.zeros(self.n)
		for i in range(self.n): #i = 1,n
			y[i] = z[i]/zmax[i]
		t1 = np.zeros(self.n)
		k = 4
		rsum = np.zeros(k)
		for i in range(k): #i = 1,k
			gsum1 = 0.0
			gsum2 = 0.0
			for j in range(i+1,self.n): #j = i+1,n
				gsum1 += w[j]*y[j]
				gsum2 += w[j]
			rsum[i-k] = gsum1 / gsum2
		for i in range(k): #i = 1,k
			t1[i] = y[i]**(BB+(CC-BB)*(AA-(1.0-2.0*rsum[i])*np.abs(np.floor(0.50-rsum[i])+AA)))
		for i in range(k,self.n): #i = k+1,n
			t1[i] = y[i]
		t2 = np.zeros(self.n)

		for i in range(k): #i = 1,k
			t2[i] = 1+(np.abs(t1[i]-AAA)-BBB)*((np.floor(t1[i]-AAA+BBB)*(1.0-CCC+(AAA-BBB)/BBB))/(AAA-BBB))
			t2[i] += (np.abs(t1[i]-AAA)-BBB)*((np.floor(AAA+BBB-t1[i])*(1.0-CCC+(1.0-AAA-BBB)/BBB))/(1.0-AAA-BBB)+1.0/BBB)
		for i in range(k,self.n): #i = k+1,n
			t2[i] = 1.0+np.cos((4.0*AAAA+2.0)*np.pi*(0.5-np.abs(t1[i]-CCCC)/(2.0*(np.floor(CCCC-t1[i])+CCCC))))
			t2[i] = t2[i]+4.0*BBBB*(np.abs(t1[i]-CCCC)/(2.0*(np.floor(CCCC-t1[i])+CCCC)))**2
			t2[i] = t2[i]/(BBBB+2.0)

		t3 = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			gsum1 = 0.0
			for ii in range(2*i,2*i+2): #ii = ((i-1)*k/(M-1)+1),(i*k/(M-1))
				gsum1 += t2[ii]
				for jj in range(1): #jj = 0, (k/(M-1)-2)
					gsum1 += np.abs(t2[ii]-t2[(i*2)+((ii+jj-(i*2)+1) % 2)])

			gsum1=gsum1/3.0
			t3[i] = gsum1

		gsum1 = 0.0
		for ii in range(k,self.n): #ii = k+1,n
			gsum1 += t2[ii]
			for jj in range(3): #jj = 0, 2
				gsum1 += np.abs(t2[ii]-t2[k+((ii+jj-(k+1)+1) % (self.n-k))] )

		gsum1 = gsum1/10.0
		t3[self.q-1] = gsum1

		x = np.zeros(self.q)
		for i in range(self.q-1): #i = 1,M-1
			x[i] = np.maximum(t3[self.q-1],A[i])*(t2[i]-0.5)+0.5
		x[self.q-1] = t3[self.q-1]
		h = np.zeros(self.q)
		h[0] = 1.0
		for i in range(self.q-1): #i = 1,M-1
			h[0] = h[0]*np.sin(x[i]*np.pi/2.0)
		for j in range(1,self.q-1): #j = 2,M-1
			h[j] = 1.0
			for i in range(self.q-j): #i = 1,M-j
				h[j] = h[j]*np.sin(x[i]*np.pi/2.0)
			h[j] = h[j]*np.cos(x[self.q-j]*np.pi/2.0)
		h[self.q-1] = np.cos(x[0]*np.pi/2.0)

		f = np.zeros(self.q)
		for i in range(self.q): #i = 1,M
			f[i] = x[self.q-1]+S[i]*h[i]

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_ZDT1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'ZDT1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]
		f2 = 0.0

		for i in range(1,self.n): #i=2, n
			f2 += x[i]

		f2 = 1.0 + 9.0/float(self.n-1) * f2

		h = 1.0 - np.sqrt(f1/f2)

		f2 = f2*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_ZDT2():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'ZDT2'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]
		f2 = 0.0

		for i in range(1,self.n): #i=2, n
			f2 += x[i]

		f2 = 1.0 + 9.0/float(self.n-1) * f2

		h = 1.0 - (f1/f2)**2

		f2 = f2*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_ZDT3():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 30
		self.n     = n
		self.q     = 2

		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		x0 = self.startp()

		self.name  = 'ZDT3'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]
		f2 = 0.0

		for i in range(1,self.n): #i=2, n
			f2 += x[i]

		f2 = 1.0 + 9.0/float(self.n-1) * f2

		h = 1.0 - np.sqrt(f1/f2) - (f1/f2) * np.sin(10.0*np.pi*f1)

		f2 = f2*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_ZDT4():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 2

		self.lb    = -5.0*np.ones(self.n); self.lb[0] = 0.0
		self.ub    =  5.0*np.ones(self.n); self.ub[0] = 1.0

		x0 = self.startp()

		self.name  = 'ZDT4'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = x[0]
		f2 = 1.0+10.0*float(self.n-1)

		for i in range(1,self.n): #i=2, n
			f2 += x[i]**2-10.0*np.cos(4.0*np.pi*x[i])

		h = 1.0 - np.sqrt(f1/f2)

		f2 = f2*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_ZDT6():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 2

		self.lb    = -5.0*np.ones(self.n); self.lb[0] = 0.0
		self.ub    =  5.0*np.ones(self.n); self.ub[0] = 1.0

		x0 = self.startp()

		self.name  = 'ZDT6'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f1 = 1.0 - np.exp(-4.0*x[0])*np.sin(6.0*np.pi*x[0])**6
		f2 = 0.0

		for i in range(1,self.n): #i=2, n
			f2 += x[i]

		f2 = ( f2/float(self.n-1) )**0.25
		f2 = 9.0 * f2
		f2 = f2 + 1.0

		h = 1.0 - (f1/f2)**2

		f2 = f2*h

		return np.array([f1, f2])

	def get_bounds(self):

		return (self.lb,self.ub)

class problem_ZLT1():
	'''
	Problemi della collezione di DMS
	'''
	n     = 0
	q     = 0
	def __init__(self):
		n = 10
		self.n     = n
		self.q     = 3

		self.lb    = -1000.0*np.ones(self.n)
		self.ub    =  1000.0*np.ones(self.n)

		x0 = self.startp()

		self.name  = 'ZLT1'

	def getdim(self):
		'''
		ritorna la tuple (n,q)
		'''
		return (self.n,self.q)

	def startp(self):
		return np.zeros(self.n)

	def functs(self,x):

		f = np.zeros(self.q)
		for k in range(self.q): #k=1,M
			f[k] = (x[k]-1.0)**2

			for i in range(self.n): #i=1, n
				if not(i==k):
					f[k] += x[i]**2

		return f

	def get_bounds(self):

		return (self.lb,self.ub)

class DMSproblems_list():
	def __init__(self):
		self.list = [ \
		'bk1',	'cl1', 'deb41',	'deb53',	'deb512a',	'deb512b',	'deb512c', \
		'deb513',	'deb521a',	'deb521b',	'dg01',	'dpam1',	'dtlz1',	'dtlz1n2', \
		'dtlz2',	'dtlz2n2',	'dtlz3',	'dtlz3n2',	'dtlz4',	'dtlz4n2',	'dtlz5', \
		'dtlz5n2',	'dtlz6',	'dtlz6n2',	'ex005',	'far1',	'fes1',	'fes2',	'fes3', \
		'fonseca',	'i1',	'i2',	'i3',	'i4',	'i5',	'ikk1',	'im1',	'jin1', \
		'jin2',	'jin3',	'jin4',	'kursawe',	'l1zdt4',	'l2zdt1',	'l2zdt2',	'l2zdt3', \
		'l2zdt4',	'l2zdt6',	'l3zdt1',	'l3zdt2',	'l3zdt3',	'l3zdt4',	'l3zdt6', \
		'le1',	'lovison1',	'lovison2',	'lovison3',	'lovison4',	'lovison5',	'lovison6', \
		'lrs1',	'mhhm1',	'mhhm2',	'mlf1',	'mlf2',	'mop1',	'mop2',	'mop3',	'mop4', \
		'mop5',	'mop6',	'mop7',	'oka1',	'oka2',	'qv1',	'sch1',	'sk1',	'sk2', 'sp1',		\
		'ssfyy1', 'ssfyy2', 'tkly1', 'vfm1', 'vu1', 'vu2', 'wfg1', 'wfg2', 'wfg3', 'wfg4', \
		'wfg5', 'wfg6', 'wfg7', 'wfg8', 'wfg9', 'zdt1', 'zdt2', 'zdt3', 'zdt4', 'zdt6', 'zlt1' \
		]

#==============================================================

'''
subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 8
	m = 0
	q = 3

	return
end subroutine setdim

subroutine startp(n,x)
	implicit none
	integer	:: n
	real*8		:: x(n), l(n), u(n)

	!call setbounds(n,l,u)

	x = 0.d0

	return
end subroutine startp

subroutine functs(n,z,M,f)
	implicit none
	integer	:: n, M, i, j
	real*8		:: z(n), f(M)
	integer, parameter :: k = 4
	integer, parameter :: l = 4
	real*8, parameter :: pi = 4.d0*atan(1.d0)
	real*8, parameter :: pi2 = 2.d0*atan(1.d0)
	real*8		:: S(3), zmax(8), A(2), y(8), h(3)
	real*8		:: t1(8), t2(8), t3(8), w(8), t4(3), x(3)
!	real*8		:: S(M), zmax(n), A(M-1), y(n), h(M)
!	real*8		:: t1(n), t2(n), t3(n), w(n), t4(M), x(M)
	real*8		:: gsum1, gsum2
	real*8, parameter :: AA = 0.8d0;
	real*8, parameter :: BB = 0.75d0;
	real*8, parameter :: CC = 0.85d0;
	real*8, parameter :: AAA = 0.02d0;
	real*8, parameter :: alpha = 1.d0
	real*8, parameter :: AAAA = 5.d0

!param S {m in 1..M} := 2*m;
	do i = 1,M
		S(i) = 2.d0*dble(i)
	enddo
!param zmax {i in 1..n} := 2*i;
	do i = 1,n
		zmax(i) = 2.d0*dble(i)
	enddo
!param A {i in 1..M-1} := 1;
	do i = 1,M-1
		A(i) = 1.d0
	enddo

!# problem variables
!# transform z into [0,1] set
!var y{i in 1..n} = z[i]/zmax[i];
	do i = 1,n
		y(i) = z(i)/zmax(i)
	enddo

!# first level mapping
!var t1{i in 1..n} = if i <= k then y[i]
!      else abs(y[i]-0.35)/abs(floor(0.35-y[i])+0.35);
	do i = 1,k
		t1(i) = y(i)
	enddo
	do i = k+1,n
		t1(i) = abs(y(i)-0.35d0)/abs(floor(0.35d0-y(i))+0.35d0)
	enddo

!# second level mapping
!var t2{i in 1..n} = if i<=k then t1[i]
!    else AA+ min(0,floor(t1[i]-BB))*(AA*(BB-t1[i]))/BB-min(0,floor(CC-t1[i]))*(1-AA)*(t1[i]-CC)/(1-CC);
	do i = 1,k
		t2(i) = t1(i)
	enddo
	do i = k+1,n
		t2(i) = AA + min(0.d0,dble(floor(t1(i)-BB)))*(AA*(BB-t1(i))) / &
			BB - min(0.d0,dble(floor(CC-t1(i))))*(1.d0-AA)*(t1(i)-CC)/(1.d0-CC)
	enddo

!# third level mapping
!param w{i in 1..n} := 2*i;
!var t3{i in 1..n} = t2[i]^AAA;
	do i = 1,n
		t3(i) = t2(i)**AAA
		w(i) = 2.d0*dble(i)
	enddo

!# forth level mapping
!var t4{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t3[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
!    else (sum {j in k+1..n} (w[j]*t3[j]))/(sum {j in k+1..n} w[j]);
	do i = 1,M-1
		gsum1 = 0.d0
		gsum2 = 0.d0
		do j = ((i-1)*k/(M-1)+1), (i*k/(M-1))
			gsum1 = gsum1 + (w(j)*t3(j))
		enddo
		do j = ((i-1)*k/(M-1)+1), (i*k/(M-1))
			gsum2 = gsum2 + w(j)
		enddo
		t4(i) = gsum1/gsum2
	enddo
	gsum1 = 0.d0
	gsum2 = 0.d0
	do j = k+1,n
		gsum1 = gsum1 + (w(j)*t3(j))
		gsum2 = gsum2 + w(j)
	enddo
	t4(M) = gsum1/gsum2

!# Define objective function variables
!var x{i in 1..M} = if i<=M-1 then max(t4[M],A[i])*(t4[i]-0.5)+0.5
!    else t4[M];
	do i = 1,M-1
		x(i) = max(t4(M),A(i))*(t4(i)-0.5d0)+0.5d0
	enddo
	x(M) = t4(M)

!# Define objective function function h
!var h{m in 1..M} = if m==1 then prod {i in 1..M-1} (1-cos(x[i]*pi2))
!    else if m<=M-1 then (prod {i in 1..M-m} (1-cos(x[i]*pi2)))*(1-sin(x[M-m+1]*pi2))
!        else (1-x[1]-(cos(2*AAAA*pi*x[1]+pi2))/(2*AAAA*pi))^alpha;
	h(1) = 1.d0
	do i = 1,M-1
		h(1) = h(1)*(1.d0-cos(x(i)*pi2))
	enddo
	do j = 2,M-1
		h(j) = 1.d0
		do i = 1,M-j
			h(j) = h(j)*(1.d0-cos(x(i)*pi2))
		enddo
		h(j) = h(j)*(1.d0-sin(x(M-j+1)*pi2))
	enddo
	h(M) = (1.d0-x(1)-(cos(2.d0*AAAA*pi*x(1)+pi2))/(2.d0*AAAA*pi))**alpha

!# The objective functions
!minimize fobj {m in 1..M}:
!    x[M]+S[m]*h[m];
	do i = 1,M
		f(i) = x(M)+S(i)*h(i)
	enddo

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n, i
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb = 0.d0
	do i = 1,n
		ub(i) = 2.d0*dble(i)
	enddo

	return
end subroutine setbounds
'''
