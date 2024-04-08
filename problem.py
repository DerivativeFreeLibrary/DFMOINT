import numpy as np

class problem():
	'''
	Nel problema assumiamo che le prime <ncont> variabili sono continue
	e le successive <nint> sono discrete.
	Le continue variano tra lb e ub
	Le discrete sono intere e variano tra lbint e ubint quindi assumono valori
	lbint lbint+1 lbint+2 ... ubint-1 ubint
	'''
	n     = 0
	ncont = 0
	nint  = 0
	m     = 0
	q     = 0
	def __init__(self):
		self.n     = 10
		self.ncont = 5
		self.nint  = 5
		self.m     = 9
		self.q     = 2
		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont])
		J = np.arange(len(x)-1)
		return x[J]**2 + x[J+1]**2 + x[J]*x[J+1] - 1.0

	def getdim(self):
		'''
		ritorna la tuple (n,ncont,nint,m,q)
		'''
		return (self.n,self.ncont,self.nint,self.m,self.q)

	def startp(self):
		(lb,ub) = self.get_bounds()
		return np.hstack((lb[:self.ncont],ub[self.ncont:]))
		#return np.zeros(self.n)

	def functs(self,xmix):
		n = len(xmix)
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont])
		f = np.zeros(self.q)
		for i in range(n):
			f[0] += np.sqrt(np.abs(x[i]-np.exp(((i+1)/(n))**2.0)/3.0))
			f[1] += ((x[i]-0.5*np.cos(10.0*np.pi*((i+1)/(n)))-0.5)**2.0)
		return f

	def get_bounds(self):
		lbcont = np.zeros(self.ncont)
		ubcont = np.ones(self.ncont)
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

'''
	f = 0.d0
	do i = 1,n
		f(1) = f(1) + (abs(x(i)-exp((dble(i)/dble(n))**2.d0)/3.d0)**0.5d0)
		f(2) = f(2) + ((x(i)-0.5d0*cos(10.d0*pi*(dble(i)/dble(n)))-0.5d0)**2.d0)
	enddo
'''
