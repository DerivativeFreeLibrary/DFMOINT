#from problem import problem
import numpy as np
import pygmo as pg
import subprocess
import os
#from problem_cec09 import *
import math
from filter_factory import filter_factory
from filter_factory import front_factory
#import pyomo.environ as pe
#from pyomo.opt import SolverFactory
from DMSprob_factory import *

def transform_point(y,ideal,nadir,q):

	yout = np.zeros(q)
	for i in range(q):
		if ideal[i] == nadir[i]:
			yout[i] = (y[i] - ideal[i]) #0.0 #y[i]
		else:
			yout[i] = (y[i] - ideal[i]) / (nadir[i] - ideal[i])
	return yout
def transform_front(F,ideal,nadir,q):
	Fout = np.empty((0, q))
	for i in range(len(F)):
		Fout = np.append(Fout, [transform_point(F[i,:],ideal,nadir,q)], axis=0)
	return Fout

class problem_urgenze():
	'''
	problem defined in README.md
	'''
	n     = 0
	ncont = 0
	nint  = 0
	m     = 0
	q     = 0
	def __init__(self,n=16,ctype='z',nint=16):
#		F =  10.0
#		sigma = 10.0
		assert n >= 16
		assert nint <= n
		self.n     = n
		self.ncont = n-nint
		self.nint  = nint
		assert self.ncont >= 0
		assert self.nint  >= 2

		self.q     = 2
		self.lb    =   np.array([ 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0,0])
		self.ub    =   np.array([20,20,20,20,20,20,20,20,20,20,20,20,20,20,14,6])

		self.ctype = ctype
#		assert ctype='z'

		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Urgenze'

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'z':
			return np.array([x[0]-x[1],x[2]-x[3],x[4]-x[5],x[6]-x[7],x[8]-x[9],x[10]-x[11],x[12]-x[13],\
							x[0]-x[1]+x[2]-x[3]+x[4]-x[5]+x[6]-x[7]+x[8]-x[9]+x[10]-x[11]+x[12]-x[13]+21])

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		f = np.zeros(self.q)

		vin = self.fconstriq(x)
		#print(vin)
		#input()
		if np.max(vin) <= 0.0:
			''' write variables into file variables.txt '''
			fid = open("Input/variables.txt","w")
			for i in range(n):
				print("%2d" % x[i],file=fid)
#			print("%f" % x[0],file=fid)
			fid.close()
			print('scritte le variabili')
			#input()

			# commentato da Stefano il 26-07-2022 per prova
			# ''' launch simulator '''
			# subprocess.call('lancio_simulatore.cmd')

			# ''' read outputs from file '''
			# outs = np.zeros(9)
			# i = 0
			# with open("Output/responses.txt","r") as fid:
				# line = fid.readline()
				# while line != '' and i < 9:
					# outs[i] = float(line)
					# print("%2d" % outs[i])
					# i += 1
					# line = fid.readline()
			# #input()
			# f1 = np.sum(outs)
			# fine commentato da Stefano il 26-07-2022 per prova
			# f1 di appoggio Stefano 26-07-2022
			f1 = -x[1]*x[1]-x[2]*x[2]
			# fine f1 di appoggio Stefano 26-07-2022
			f2 =  x[1]- x[0]
			f2+=  x[3]- x[2]
			f2+=  x[5]- x[4]
			f2+=  x[7]- x[6]
			f2+= x[9]- x[8]
			f2+= x[11]-x[10]
			f2+= x[13]-x[12]
		else:
			outs = 1.e+10*np.ones(9)
			f1 = 1.e+10
			f2 = 1.e+10

		print(x)
		print(outs)
		#print('hit return to continue...')
		#input()

		return np.array([f1,f2])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = self.ub[self.ncont:] - self.lb[self.ncont:] #100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return self.m
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_2():
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
	def __init__(self,nint):
		self.n     = 4
		assert nint < self.n
		if nint < 0:
			nint = 2
		self.ncont = self.n-nint #2
		self.nint  = nint #2
		self.m     = 0 #9
		self.q     = 2
		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.B     = 100
		self.p1    = np.random.permutation(self.B+1)
		self.p2    = np.random.permutation(self.B+1)
		F          = 10.0
		sigma      = 10.0
		self.lb[0] = F/sigma
		self.ub[0] = 3*F/sigma
		self.lb[1] = np.sqrt(2.0)*F/sigma
		self.ub[1] = 3*F/sigma
		self.lb[2] = np.sqrt(2.0)*F/sigma
		self.ub[2] = 3*F/sigma
		self.lb[3] = F/sigma
		self.ub[3] = 3*F/sigma
		self.name  = 'Problem CL1'
		return

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]
		ubcont = self.ub[:self.ncont]
		#lbcont = self.lb[self.ncont:]
		#ubcont = self.ub[self.ncont:]
		lbint = np.zeros(self.nint)
		ubint = self.B*np.ones(self.nint)
		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		return np.empty(0)

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
		x[2] = self.p1[x[2].round().astype(int)]
		x[3] = self.p2[x[3].round().astype(int)]
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + x[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])
		f = np.zeros(self.q)
		F = 10.0
		E = 2.e+5
		L = 200.0
		sigma = 10.0
		#f[0] = L*(2.0*x[0] + np.sqrt(2.0)*x[1] + np.sqrt(x[2]) + x[3])
		#f[1] = F*L/E*(2.0/x[0] + (2.0*np.sqrt(2.0))/x[1] - (2.0*np.sqrt(2.0))/x[2] + 2.0/x[3])
		f[0] = L*(2.0*x[2] + np.sqrt(2.0)*x[3] + np.sqrt(x[0]) + x[1])
		f[1] = F*L/E*(2.0/x[2] + (2.0*np.sqrt(2.0))/x[3] - (2.0*np.sqrt(2.0))/x[0] + 2.0/x[1])
		return f

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		return objs

class problem_solar8():
	'''
	wrapper per il problema SOLAR 8

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
	Ic    = 0
	Iz    = 0
	def __init__(self):
		self.n     = 13
		self.ncont = 11 #2
		self.nint  = 2 #2
		self.m     = 9 #9
		self.q     = 2
		#X0               ( 11.0 11.0 200.0 10.0 10.0 2650 89.0  0.5  8.0 36 0.30 0.020 0.0216 )
		self.lbvero = np.array([ 1.0,  1.0,  20.0,  1.0,  1.0,      1,  1.0,  0.0,  1.0,      1, 0.01, 0.005, 0.0060] )
		self.ubvero = np.array([40.0, 40.0, 250.0, 30.0, 30.0,  1.e+4, 89.0, 20.0, 20.0,  1.e+4, 5.00, 0.100, 0.1000] )
		self.Ic = [0,1,2,3,4,6,7,8,10,11,12]
		self.Iz = [5,9]
		self.lb = np.zeros(self.n)
		self.ub = np.ones(self.n)
		self.lb[5] = 1;    self.lb[9] = 1;
		self.ub[5] = 10000;self.ub[9] = 10000;
		self.scale = np.array([1e+12, 1.e+8, 1.e+4, 1.e+2, 1.0, 1.0, 1.e+8, 1.e-3, 1.e+1, 1.e+12, 1.e-1])

		x0 = self.startp()

		lb,ub = self.get_bounds()

		#for i in range(self.n):
		#	print(lb[i],x0[i],ub[i])

		self.name  = 'Problem SOLAR 8'

		self.eps = 0.1*np.ones((self.q,self.m))

		fob,ciq = self.run_bb(x0)
		#fob = self.functs(x0)
		#ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3*self.scale[self.q+i]
				else:
					self.eps[k,i] = 1.e-1*self.scale[self.q+i]

		return

	def get_bounds(self):
		Ic = self.Ic
		Iz = self.Iz
		lbcont = self.lb[Ic]
		ubcont = self.ub[Ic]
		#lbcont = self.lb[self.ncont:]
		#ubcont = self.ub[self.ncont:]
		lbint = self.lb[Iz]
		ubint = self.ub[Iz]
		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		n = len(xmix)
		Ic = self.Ic
		Iz = self.Iz
		x = np.zeros(n)
		x[Ic] = self.lbvero[Ic] + xmix[:self.ncont]*(self.ubvero[Ic]-self.lbvero[Ic])
		x[Iz] = xmix[self.ncont:]
		g = np.zeros(self.m)

		with open("x_solar8.txt","w") as finputx:
			for i in range(n):
				if i in Iz:
					print("%11d " % np.round(x[i]),end='',file=finputx)
				else:
					print("%18.9e " % x[i],end='',file=finputx)

		#os.system("/home/giampo/solar/bin/solar 8 x_solar8.txt > out.txt")
		os.system("solar_WINDOWS.exe 8 x_solar8.txt -prec=0.1 > out.txt")
		#os.system("./solar 8 -fid=0.1 x_solar8.txt > out.txt")
		
		with open("out.txt","r") as fo:
			line = fo.readline()
			arr = line.split(" ")

			for i in range(self.m): #elem in arr:
				g[i] = float(arr[self.q+i])/self.scale[self.q+i]

			#	print(g[i],end='')
			#	print(" ",end='')
			#print()

		#print("g = ",g)
		#print("Hit RETURN to continue")
		#input()
		return g

	def getdim(self):
		'''
		ritorna la tuple (n,ncont,nint,m,q)
		'''
		return (self.n,self.ncont,self.nint,self.m,self.q)

	def startp(self):
		(lb,ub) = self.get_bounds()
		Ic = self.Ic
		x0 = np.array([11.0, 11.0, 200.0, 10.0, 10.0, 89.0,  0.5,  8.0, 0.30, 0.020, 0.0216, 2650, 36 ])
		x = np.copy(x0)
		x[:self.ncont] = (x0[:self.ncont] - self.lbvero[Ic])/(self.ubvero[Ic]-self.lbvero[Ic])
        
		return x    
		return x0
		return np.hstack((lb[:self.ncont],ub[self.ncont:]))
		#X0 ( 11.0 11.0 200.0 10.0 10.0 2650 89.0  0.5  8.0 36 0.30 0.020 0.0216 )
		#return np.zeros(self.n)

	def functs(self,xmix):
		n = len(xmix)
		Ic = self.Ic
		Iz = self.Iz
		x = np.zeros(n)
		x[Ic] = self.lbvero[Ic] + xmix[:self.ncont]*(self.ubvero[Ic]-self.lbvero[Ic])
		#print(self.lbvero[Ic])
		#print(self.ubvero[Ic]-self.lbvero[Ic])
		#print(xmix[:self.ncont])
		#input()
        
		x[Iz] = xmix[self.ncont:]
		f = np.zeros(self.q)

		with open("x_solar8.txt","w") as finputx:
			for i in range(n):
				if i in Iz:
					print("%11d " % np.round(x[i]),end='',file=finputx)
				else:
					print("%18.9e " % x[i],end='',file=finputx)

		#os.system("/home/giampo/solar/bin/solar 8 x_solar8.txt > out.txt")
		os.system("solar_WINDOWS.exe 8 x_solar8.txt -prec=0.1 > out.txt")
		#os.system("./solar 8 -fid=0.1 x_solar8.txt > out.txt")
		
		with open("out.txt","r") as fo:
			line = fo.readline()
			arr = line.split()

			for i in range(self.q): #elem in arr:
				f[i] = float(arr[i])/self.scale[i]
                
			#	print(f[i],end='')
			#	print(" ",end='')
			#print()

		#print("f = ",f)
		#print("Hit RETURN to continue")
		#input()
		return f
	
	def run_bb(self,xmix):
		n = len(xmix)
		Ic = self.Ic
		Iz = self.Iz
		x = np.zeros(n)
		x[Ic] = self.lbvero[Ic] + xmix[:self.ncont]*(self.ubvero[Ic]-self.lbvero[Ic])
		x[Iz] = xmix[self.ncont:]

		f = np.zeros(self.q)
		g = np.zeros(self.m)

		with open("x_solar8.txt","w") as finputx:
			for i in range(n):
				if i in Iz:
					print("%11d " % np.round(x[i]),end='',file=finputx)
				else:
					print("%18.9e " % x[i],end='',file=finputx)

		#os.system("/home/giampo/solar/bin/solar 8 x_solar8.txt > out.txt")
		os.system("solar_WINDOWS.exe 8 x_solar8.txt -prec=0.1 > out.txt")
		#os.system("./solar 8 -fid=0.1 x_solar8.txt > out.txt")

		with open("out.txt","r") as fo:
			line = fo.readline()
			arr = line.split()

			for i in range(self.q): #elem in arr:
				f[i] = float(arr[i])/self.scale[i]
				
			for i in range(self.m): #elem in arr:
				g[i] = float(arr[self.q+i])#/self.scale[self.q+i]
                
			#	print(f[i],end='')
			#	print(" ",end='')
			#print()

		#print("f = ",f)
		#print("Hit RETURN to continue")
		#input()
		return f,g

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_solar9():
	'''
	wrapper per il problema SOLAR 9

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
	Ic    = 0
	Iz    = 0
	def __init__(self):
		self.n     = 29
		self.ncont = 22 #2
		self.nint  = 7 #2
		self.m     = 17 #9
		self.q     = 2
		self.lb    = np.array([ 1.0,  1.0,  20.0,  1.0,  1.0,    1,  1.0,  0.0,  1.0, 793.0,  1.0,  1.0, 0.01, 0.01, 495.0,    1, 0.01, 0.0050, 0.006, 0.007,  0.5, 0.0050, 0.006, 0.15,     2,     1,  1, 1, 1] )
		self.ub    = np.array([40.0, 40.0, 250.0, 30.0, 30.0,1.e+4, 89.0, 20.0, 20.0, 995.0, 50.0, 30.0, 5.00, 5.00, 650.0,1.e+4, 5.00, 0.1000, 0.100, 0.200, 10.0, 0.1000, 0.100, 0.40, 1.e+3, 1.e+6, 10, 9, 8] )
		self.Ic = [0,1,2,3,4,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23]
		self.Iz = [5,15,24,25,26,27,28]

		x0 = self.startp()

		lb,ub = self.get_bounds()

		#for i in range(self.n):
		#	print(lb[i],x0[i],ub[i])

		self.name  = 'Problem SOLAR 9'

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def get_bounds(self):
		Ic = self.Ic
		Iz = self.Iz
		lbcont = self.lb[Ic]
		ubcont = self.ub[Ic]
		#lbcont = self.lb[self.ncont:]
		#ubcont = self.ub[self.ncont:]
		lbint = self.lb[Iz]
		ubint = self.ub[Iz]
		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		n = len(xmix)
		Ic = self.Ic
		Iz = self.Iz
		x = np.zeros(n)
		x[Ic] = xmix[:self.ncont]
		x[Iz] = xmix[self.ncont:]
		g = np.zeros(self.m)

		with open("x_solar9.txt","w") as finputx:
			for i in range(n):
				if i in Iz:
					print("%11d " % np.round(x[i]),end='',file=finputx)
				else:
					print("%18.9e " % x[i],end='',file=finputx)

		os.system("/home/giampo/solar/bin/solar 9 x_solar9.txt > out.txt")

		with open("out.txt","r") as fo:
			line = fo.readline()
			arr = line.split(" ")

			for i in range(self.m): #elem in arr:
				g[i] = float(arr[self.q+i])

				#print(float(elem),end='')
				#print(" ",end='')
			#print()

		#print(g)
		#print("Hit RETURN to continue")
		#input()
		return g

	def getdim(self):
		'''
		ritorna la tuple (n,ncont,nint,m,q)
		'''
		return (self.n,self.ncont,self.nint,self.m,self.q)

	def startp(self):
		(lb,ub) = self.get_bounds()

		return np.array([9.0,  9.0, 150.0,  6.0,  8.0, 45.0,  0.5,  5.0, 900.0,  9.0,  9.0, 0.30, 0.20, 560.0, 0.30, 0.0165, 0.018, 0.017, 10.0, 0.0155, 0.016, 0.20, 1000, 500, 3, 12000,  1, 2, 2 ])
		return np.hstack((lb[:self.ncont],ub[self.ncont:]))
		#X0 ( 11.0 11.0 200.0 10.0 10.0 2650 89.0  0.5  8.0 36 0.30 0.020 0.0216 )
		#return np.zeros(self.n)

	def functs(self,xmix):
		n = len(xmix)
		Ic = self.Ic
		Iz = self.Iz
		x = np.zeros(n)
		x[Ic] = xmix[:self.ncont]
		x[Iz] = xmix[self.ncont:]
		f = np.zeros(self.q)

		with open("x_solar9.txt","w") as finputx:
			for i in range(n):
				if i in Iz:
					print("%11d " % np.round(x[i]),end='',file=finputx)
				else:
					print("%18.9e " % x[i],end='',file=finputx)

		os.system("/home/giampo/solar/bin/solar 9 x_solar9.txt > out.txt")

		with open("out.txt","r") as fo:
			line = fo.readline()
			arr = line.split()

			for i in range(self.q): #elem in arr:
				f[i] = float(arr[i])

				#print(float(elem),end='')
				#print(" ",end='')
			#print()

		#print(f)
		#print("Hit RETURN to continue")
		#input()
		return f

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_2_discretizzato():
	'''
	è il problem_2 in cui le variabili sono tutte continue anche le ultime 2
	che dovrebbero essere discrete
	'''
	n     = 0
	ncont = 0
	nint  = 0
	m     = 0
	q     = 0
	def __init__(self):
		self.n     = 4
		self.ncont = 4
		self.nint  = 0
		self.m     = 0 #9
		self.q     = 2
		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.B     = 100
		self.p1    = np.random.permutation(self.B+1)
		self.p2    = np.random.permutation(self.B+1)
		F          = 10.0
		sigma      = 10.0
		self.lb[0] = F/sigma
		self.ub[0] = 3*F/sigma
		self.lb[1] = np.sqrt(2.0)*F/sigma
		self.ub[1] = 3*F/sigma
		self.lb[2] = np.sqrt(2.0)*F/sigma
		self.ub[2] = 3*F/sigma
		self.lb[3] = F/sigma
		self.ub[3] = 3*F/sigma
		self.name  = 'Problem CL1 (discretizzato)'
		return

	def get_bounds(self):
		#lbcont = self.lb[:2]
		#ubcont = self.ub[:2]
		lbcont = self.lb[2:]
		ubcont = self.ub[2:]
		lbint = np.zeros(2)
		ubint = self.B*np.ones(2)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		return np.empty(0)

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

		#permuta le variabili intere per distruggere l'ordinamento naturale
		x[2] = self.p1[x[2].round().astype(int)]
		x[3] = self.p2[x[3].round().astype(int)]

		#riporta le ultime due variabili nei box veri
		(lb,ub) = self.get_bounds()
		x[2:] = self.lb[:2] + x[:2]*(self.ub[:2]-self.lb[:2])/(ub[2:]-lb[2:])
		f = np.zeros(self.q)
		F = 10.0
		E = 2.e+5
		L = 200.0
		sigma = 10.0
		#f[0] = L*(2.0*x[0] + np.sqrt(2.0)*x[1] + np.sqrt(x[2]) + x[3])
		#f[1] = F*L/E*(2.0/x[0] + (2.0*np.sqrt(2.0))/x[1] - (2.0*np.sqrt(2.0))/x[2] + 2.0/x[3])
		f[0] = L*(2.0*x[2] + np.sqrt(2.0)*x[3] + np.sqrt(x[0]) + x[1])
		f[1] = F*L/E*(2.0/x[2] + (2.0*np.sqrt(2.0))/x[3] - (2.0*np.sqrt(2.0))/x[0] + 2.0/x[1])
		return f

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		return objs

class problem_1():
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
		self.m     = 0 #9
		self.q     = 2
		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)
		self.name  = 'Test problem 1'

		x0 = self.startp()

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1
		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])
		J = np.arange(len(x)-1)
		return np.empty(0)
		#return x[J]**2 + x[J+1]**2 + x[J]*x[J+1] - 1.0

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
		'''
			f = 0.d0
			do i = 1,n
				f(1) = f(1) + (abs(x(i)-exp((dble(i)/dble(n))**2.d0)/3.d0)**0.5d0)
				f(2) = f(2) + ((x(i)-0.5d0*cos(10.d0*pi*(dble(i)/dble(n)))-0.5d0)**2.d0)
			enddo
		'''
		n = len(xmix)
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])
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

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

def fconstr_a(x):
    J = np.arange(len(x)-2)
    return  (3-2*x[J+1])*x[J+1] - x[J] - 2*x[J+2] + 1

def fconstr_b(x):
    J = np.arange(len(x)-2)
    return  (3-2*x[J+1])*x[J+1] - x[J] - 2*x[J+2] + 2.5

def fconstr_c(x):
    J = np.arange(len(x)-1)
    return x[J]**2 + x[J+1]**2 + x[J]*x[J+1] - 2*x[J] - 2*x[J+1] + 1

def fconstr_d(x):
    J = np.arange(len(x)-1)
    return x[J]**2 + x[J+1]**2 + x[J]*x[J+1] - 1

def fconstr_e(x):
    J = np.arange(len(x)-2)
    return (3-0.5*x[J+1])*x[J+1] - x[J] -2*x[J+2] +1

def fconstr_f(x):
    J = np.arange(len(x)-2)
    return np.array([np.sum((3-0.5*x[J+1])*x[J+1] - x[J] -2*x[J+2] +1)])

def fconstr_z(x):
    return np.array([-1.0]) #np.array([],dtype=np.float64) #np.array([-1.0])

class problem_cec09_1():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 2
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)

		self.lb    = -1.0*np.ones(self.n)
		self.lb[0] =  0.0
		self.ub    =      np.ones(self.n)

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(1) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		'''
			f = 0.d0
			do i = 1,n
				f(1) = f(1) + (abs(x(i)-exp((dble(i)/dble(n))**2.d0)/3.d0)**0.5d0)
				f(2) = f(2) + ((x(i)-0.5d0*cos(10.d0*pi*(dble(i)/dble(n)))-0.5d0)**2.d0)
			enddo
		'''
		n = len(xmix)
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(2, n, 2)
		J2 = np.arange(1, n, 2)
		obj_1 = x[0] + (2 / len(J1)) * sum(
			[(x[i] - np.sin(6 * np.pi * x[0] + (i + 1) * np.pi / n)) ** 2 for i in J1])
		obj_2 = 1 - np.sqrt(x[0]) + (2 / len(J2)) * sum(
			[(x[i] - np.sin(6 * np.pi * x[0] + (i + 1) * np.pi / n)) ** 2 for i in J2])

		self.FRONT.update_point([obj_1,obj_2])
		#print('F:',self.FRONT.F)
		#input()
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
			#print('In try')
			#print(self.ref_point,self.idl_point)
			#print(FT)
			#print(nadir)
		except:
			#print('In except')
			self.HV.append(0.0)
		#print(self.HV[-1])
		#input()

		return np.array([obj_1, obj_2])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_2():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 2
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = -1.0*np.ones(self.n)
		self.lb[0] =  0.0
		self.ub    =      np.ones(self.n)

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(2) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(2, n, 2)
		J2 = np.arange(1, n, 2)

		y_odd = 2 * sum([(x[j] -
			(0.3 * x[0] ** 2 * np.cos(
			24 * np.pi * x[0] + 4 * (j + 1) * np.pi / n) + 0.6 * x[0])
			* np.cos(6 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J1]) / len(J1)

		y_even = 2 * sum([(x[j] -
			(0.3 * x[0] ** 2 * np.cos(
			24 * np.pi * x[0] + 4 * (j + 1) * np.pi / n) + 0.6 * x[0])
			* np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J2]) / len(J2)

		obj_1 = x[0] + y_odd
		obj_2 = 1 - np.sqrt(x[0]) + y_even

		self.FRONT.update_point([obj_1,obj_2])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_3():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 2
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = np.zeros(self.n)
		self.ub    = np.ones(self.n)

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(3) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(2, n, 2)
		J2 = np.arange(1, n, 2)

		y = [x[j] - x[0] ** (0.5 * (1 + 3 * (j - 1) / (n - 2))) for j in np.arange(n)]
		y_1 = sum([y[j] ** 2 for j in J1])

		y_2 = np.prod([np.cos(20 * y[j] * np.pi / np.sqrt(j + 1.0)) for j in J1])

		y_odd = 2 * (4 * y_1 - 2 * y_2 + 2) / len(J1)

		obj_1 = x[0] + y_odd

		y_even = 2 * (
			4 * sum([y[j] ** 2 for j in J2])
			- 2 * np.prod([np.cos(20 * y[j] * np.pi / np.sqrt(j + 1.0)) for j in J2])
			+ 2) / len(J2)

		obj_2 = 1 - np.sqrt(x[0]) + y_even

		self.FRONT.update_point([obj_1,obj_2])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_4():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 2
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = -2.0*np.ones(self.n)
		self.ub    =  2.0*np.ones(self.n)
		self.lb[0] = 0.0
		self.ub[0] = 1.0

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(4) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(2, n, 2)
		J2 = np.arange(1, n, 2)

		y = [x[j] - np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / n) for j in np.arange(n)]

		obj_1 = x[0] + 2 * sum([np.abs(y[j] / (1 + np.exp(2 * np.abs(y[j])))) for j in J1]) / len(J1)
		obj_2 = 1 - x[0] ** 2 + 2 * sum(
			[np.abs(y[j] / (1 + np.exp(2 * np.abs(y[j])))) for j in J2]) / len(J2)

		self.FRONT.update_point([obj_1,obj_2])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_5():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 2
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = -1.0*np.ones(self.n)
		self.ub    =      np.ones(self.n)
		self.lb[0] = 0.0

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(5) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(2, n, 2)
		J2 = np.arange(1, n, 2)

		N = 10
		eps = 0.1

		def h(t):
			return 2 * t ** 2 - np.cos(4 * np.pi * t) + 1

		y = [x[j] - np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / n) for j in np.arange(n)]

		obj_1 = x[0] + (1 / (2 * N) + eps) * np.abs(np.sin(2 * N * np.pi * x[0])) + 2 * sum(
			[h(y[j]) for j in J1]) / len(J1)
		obj_2 = 1 - x[0] + (1 / (2 * N) + eps) * np.abs(
			np.sin(2 * N * np.pi * x[0])) + 2 * sum([h(y[j]) for j in J2]) / len(J2)

		self.FRONT.update_point([obj_1,obj_2])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_6():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 2
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = -1.0*np.ones(self.n)
		self.ub    =      np.ones(self.n)
		self.lb[0] = 0.0

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(6) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(2, n, 2)
		J2 = np.arange(1, n, 2)

		N = 2
		eps = 0.1
		y = [x[j] - np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / self.n) for j in np.arange(n)]

		obj_1 = x[0] + (
			np.maximum(0.0, 2 * (1. / (2 * N) + eps) * np.sin(2 * N * np.pi * x[0]))
			+ 2 * (
			4 * sum([y[j] ** 2 for j in J1]) - 2 * np.prod(
			[np.cos(20 * y[j] * np.pi / np.sqrt(j + 1)) for j in J1])
			+ 2) / len(J1))

		obj_2 = 1 - x[0] + (
			np.maximum(0.0, 2 * (1. / (2 * N) + eps) * np.sin(2 * N * np.pi * x[0]))
			+ 2 * (
			4 * sum([y[j] ** 2 for j in J2]) - 2 * np.prod(
			[np.cos(20 * y[j] * np.pi / np.sqrt(j + 1)) for j in J2])
			+ 2) / len(J2))

		self.FRONT.update_point([obj_1,obj_2])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_7():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 2
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = -1.0*np.ones(self.n)
		self.ub    =      np.ones(self.n)
		self.lb[0] = 0.0

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(7) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(2, n, 2)
		J2 = np.arange(1, n, 2)

		y = [x[j] - np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / self.n) for j in np.arange(n)]
		obj_1 = x[0] ** (1 / 5) + (2 / len(J1)) * sum([y[j] ** 2 for j in J1])
		obj_2 = 1 - x[0] ** (1 / 5) + (2 / len(J2)) * sum([y[j] ** 2 for j in J2])

		self.FRONT.update_point([obj_1,obj_2])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_8():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 3
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = -2.0*np.ones(self.n)
		self.ub    =  2.0*np.ones(self.n)
		self.lb[0] = 0.0
		self.lb[1] = 0.0
		self.ub[0] = 1.0
		self.ub[1] = 1.0

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(8) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(3, n, 3)
		J2 = np.arange(4, n, 3)
		J3 = np.arange(2, n, 3)

		obj_1 = np.cos(0.5 * np.pi * x[0]) * np.cos(0.5 * np.pi * x[1]) + (2 / len(J1)) * sum(
			[(x[j] - 2 * x[1] * np.sin(2 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J1]
		)
		obj_2 = np.cos(0.5 * np.pi * x[0]) * np.sin(0.5 * np.pi * x[1]) + (2 / len(J2)) * sum(
			[(x[j] - 2 * x[1] * np.sin(2 * np.pi * x[0] + (j + 1) * np.pi / self.n)) ** 2 for j in J2]
		)
		obj_3 = np.sin(0.5 * np.pi * x[0]) + (2 / len(J3)) * sum(
			[(x[j] - 2 * x[1] * np.sin(2 * np.pi * x[0] + (j + 1) * np.pi / self.n)) ** 2 for j in J3]
		)

		self.FRONT.update_point([obj_1,obj_2,obj_3])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2, obj_3])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_9():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 3
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = -2.0*np.ones(self.n)
		self.ub    =  2.0*np.ones(self.n)
		self.lb[0] = 0.0
		self.lb[1] = 0.0
		self.ub[0] = 1.0
		self.ub[1] = 1.0

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(9) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1

		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(3, n, 3)
		J2 = np.arange(4, n, 3)
		J3 = np.arange(2, n, 3)
		eps = 0.1

		obj_1 = 0.5 * (np.maximum(0.0, (1 + eps) * (1 - 4 * (2 * x[0] - 1) ** 2))
			+ 2 * x[0]) * x[1] + (2 / len(J1)) * sum([(x[j] - 2 * x[1] * np.sin(
			2 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J1])

		obj_2 = 0.5 * (np.maximum(0.0, (1 + eps) * (1 - 4 * (2 * x[0] - 1) ** 2))
			- 2 * x[0] + 2) * x[1] + (2 / len(J2)) * sum([(x[j] - 2 * x[1] * np.sin(
			2 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J2])

		obj_3 = 1 - x[1] + (2 / len(J3)) * sum([(x[j] - 2 * x[1] * np.sin(
			2 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J3])

		self.FRONT.update_point([obj_1,obj_2,obj_3])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2, obj_3])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_cec09_10():
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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,n=10,ctype='a',nint=-1):
		assert n >= 4
		assert nint < n
		self.n     = n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			assert self.ncont >= 2
			assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		self.q     = 3
		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)
		self.lb    = -2.0*np.ones(self.n)
		self.ub    =  2.0*np.ones(self.n)
		self.lb[0] = 0.0
		self.lb[1] = 0.0
		self.ub[0] = 1.0
		self.ub[1] = 1.0

		self.ctype = ctype
		x0 = self.startp()
		self.m     = self.fconstriq(x0).shape[0]

		self.name  = 'Problem CEC09(10) '+str(n)+' '+ctype

		self.eps = 0.1*np.ones((self.q,self.m))

		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1
		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		J1 = np.arange(3, n, 3)
		J2 = np.arange(4, n, 3)
		J3 = np.arange(2, n, 3)
		y = [x[j] - 2 * x[1] * np.sin(2 * np.pi * x[0] + (j + 1) * np.pi / self.n) for j in np.arange(n)]

		obj_1 = np.cos(0.5 * x[0] * np.pi) * np.cos(0.5 * x[1] * np.pi) + (2 / len(J1)) * sum(
			[4 * y[j] ** 2 - np.cos(8 * np.pi * y[j]) + 1 for j in J1])

		obj_2 = np.cos(0.5 * x[0] * np.pi) * np.sin(0.5 * x[1] * np.pi) + (2 / len(J2)) * sum(
			[4 * y[j] ** 2 - np.cos(8 * np.pi * y[j]) + 1 for j in J2])

		obj_3 = np.sin(0.5 * x[0] * np.pi) + (2 / len(J3)) * sum(
			[4 * y[j] ** 2 - np.cos(8 * np.pi * y[j]) + 1 for j in J3])

		self.FRONT.update_point([obj_1,obj_2,obj_3])
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])

		return np.array([obj_1, obj_2, obj_3])

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_luis():
	'''
	Problemi della collezione di DMS

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
	HV    = []
	ref_point = []
	idl_point = []
	def __init__(self,prbsel='bk1',ctype='a',nint=-1):
		#print(prbsel)
		if prbsel=='bk1':
			prb = problem_BK1()
		elif prbsel == 'cl1':
			prb = problem_CL1()
		elif prbsel == 'deb41':
			prb = problem_DEB41()
		elif prbsel == 'deb53':
			prb = problem_DEB53()
		elif prbsel == 'deb512a':
			prb = problem_DEB512a()
		elif prbsel == 'deb512b':
			prb = problem_DEB512b()
		elif prbsel == 'deb512c':
			prb = problem_DEB512c()
		elif prbsel == 'deb513':
			prb = problem_DEB513()
		elif prbsel == 'deb521a':
			prb = problem_DEB521a()
		elif prbsel == 'deb521b':
			prb = problem_DEB521b()
		elif prbsel == 'dg01':
			prb = problem_DG01()
		elif prbsel == 'dpam1':
			prb = problem_DPAM1()
		elif prbsel == 'dtlz1':
			prb = problem_DTLZ1()
		elif prbsel == 'dtlz1n2':
			prb = problem_DTLZ1n2()
		elif prbsel == 'dtlz2':
			prb = problem_DTLZ2()
		elif prbsel == 'dtlz2n2':
			prb = problem_DTLZ2n2()
		elif prbsel == 'dtlz3':
			prb = problem_DTLZ3()
		elif prbsel == 'dtlz3n2':
			prb = problem_DTLZ3n2()
		elif prbsel == 'dtlz4':
			prb = problem_DTLZ4()
		elif prbsel == 'dtlz4n2':
			prb = problem_DTLZ4n2()
		elif prbsel == 'dtlz5':
			prb = problem_DTLZ5()
		elif prbsel == 'dtlz5n2':
			prb = problem_DTLZ5n2()
		elif prbsel == 'dtlz6':
			prb = problem_DTLZ6()
		elif prbsel == 'dtlz6n2':
			prb = problem_DTLZ6n2()
		elif prbsel == 'ex005':
			prb = problem_EX005()
		elif prbsel == 'far1':
			prb = problem_FAR1()
		elif prbsel == 'fes1':
			prb = problem_FES1()
		elif prbsel == 'fes2':
			prb = problem_FES2()
		elif prbsel == 'fes3':
			prb = problem_FES3()
		elif prbsel == 'fonseca':
			prb = problem_FONSECA()
		elif prbsel == 'i1':
			prb = problem_I1()
		elif prbsel == 'i2':
			prb = problem_I2()
		elif prbsel == 'i3':
			prb = problem_I3()
		elif prbsel == 'i4':
			prb = problem_I4()
		elif prbsel == 'i5':
			prb = problem_I5()
		elif prbsel == 'ikk1':
			prb = problem_IKK1()
		elif prbsel == 'im1':
			prb = problem_IM1()
		elif prbsel == 'jin1':
			prb = problem_JIN1()
		elif prbsel == 'jin2':
			prb = problem_JIN2()
		elif prbsel == 'jin3':
			prb = problem_JIN3()
		elif prbsel == 'jin4':
			prb = problem_JIN4()
		elif prbsel == 'kursawe':
			prb = problem_KURSAWE()
		elif prbsel == 'l1zdt4':
			prb = problem_L1ZDT4()
		elif prbsel == 'l2zdt1':
			prb = problem_L2ZDT1()
		elif prbsel == 'l2zdt2':
			prb = problem_L2ZDT2()
		elif prbsel == 'l2zdt3':
			prb = problem_L2ZDT3()
		elif prbsel == 'l2zdt4':
			prb = problem_L2ZDT4()
		elif prbsel == 'l2zdt6':
			prb = problem_L2ZDT6()
		elif prbsel == 'l3zdt1':
			prb = problem_L3ZDT1()
		elif prbsel == 'l3zdt2':
			prb = problem_L3ZDT2()
		elif prbsel == 'l3zdt3':
			prb = problem_L3ZDT3()
		elif prbsel == 'l3zdt4':
			prb = problem_L3ZDT4()
		elif prbsel == 'l3zdt6':
			prb = problem_L3ZDT6()
		elif prbsel == 'le1':
			prb = problem_LE1()
		elif prbsel == 'lovison1':
			prb = problem_LOVISON1()
		elif prbsel == 'lovison2':
			prb = problem_LOVISON2()
		elif prbsel == 'lovison3':
			prb = problem_LOVISON3()
		elif prbsel == 'lovison4':
			prb = problem_LOVISON4()
		elif prbsel == 'lovison5':
			prb = problem_LOVISON5()
		elif prbsel == 'lovison6':
			prb = problem_LOVISON6()
		elif prbsel == 'lrs1':
			prb = problem_LRS1()
		elif prbsel == 'mhhm1':
			prb = problem_MHHM1()
		elif prbsel == 'mhhm2':
			prb = problem_MHHM2()
		elif prbsel == 'mlf1':
			prb = problem_MLF1()
		elif prbsel == 'mlf2':
			prb = problem_MLF2()
		elif prbsel == 'mop1':
			prb = problem_MOP1()
		elif prbsel == 'mop2':
			prb = problem_MOP2()
		elif prbsel == 'mop3':
			prb = problem_MOP3()
		elif prbsel == 'mop4':
			prb = problem_MOP4()
		elif prbsel == 'mop5':
			prb = problem_MOP5()
		elif prbsel == 'mop6':
			prb = problem_MOP6()
		elif prbsel == 'mop7':
			prb = problem_MOP7()
		elif prbsel == 'oka1':
			prb = problem_OKA1()
		elif prbsel == 'oka2':
			prb = problem_OKA2()
		elif prbsel == 'qv1':
			prb = problem_QV1()
		elif prbsel == 'sch1':
			prb = problem_SCH1()
		elif prbsel == 'sk1':
			prb = problem_SK1()
		elif prbsel == 'sk2':
			prb = problem_SK2()
		elif prbsel == 'sp1':
			prb = problem_SP1()
		elif prbsel == 'ssfyy1':
			prb = problem_SSFYY1()
		elif prbsel == 'ssfyy2':
			prb = problem_SSFYY2()
		elif prbsel == 'tkly1':
			prb = problem_TKLY1()
		elif prbsel == 'vfm1':
			prb = problem_VFM1()
		elif prbsel == 'vu1':
			prb = problem_VU1()
		elif prbsel == 'vu2':
			prb = problem_VU2()
		elif prbsel == 'wfg1':
			prb = problem_WFG1()
		elif prbsel == 'wfg2':
			prb = problem_WFG2()
		elif prbsel == 'wfg3':
			prb = problem_WFG3()
		elif prbsel == 'wfg4':
			prb = problem_WFG4()
		elif prbsel == 'wfg5':
			prb = problem_WFG5()
		elif prbsel == 'wfg6':
			prb = problem_WFG6()
		elif prbsel == 'wfg7':
			prb = problem_WFG7()
		elif prbsel == 'wfg8':
			prb = problem_WFG8()
		elif prbsel == 'wfg9':
			prb = problem_WFG9()
		elif prbsel == 'zdt1':
			prb = problem_ZDT1()
		elif prbsel == 'zdt2':
			prb = problem_ZDT2()
		elif prbsel == 'zdt3':
			prb = problem_ZDT3()
		elif prbsel == 'zdt4':
			prb = problem_ZDT4()
		elif prbsel == 'zdt6':
			prb = problem_ZDT6()
		elif prbsel == 'zlt1':
			prb = problem_ZLT1()
		else:
			prb = problem_BK1()
		self.prb = prb
		self.list = DMSproblems_list()

		self.n, self.q = self.prb.getdim()
		n = self.n
		#assert n >= 4
		assert nint < n
		if nint < 0:
			self.ncont = math.floor(n/2)
			self.nint  = n - self.ncont
			#assert self.ncont >= 2
			#assert self.nint  >= 2
		else:
			self.ncont = n-nint
			self.nint  = nint
		#print(self.list.list)
		#print(self.prb.startp())
		#print(self.prb.functs(self.prb.startp()))

		self.FRONT = front_factory(self.q)
		self.HV.clear()
		self.ref_point = np.inf*np.ones(self.q)
		self.idl_point = np.inf*np.ones(self.q)

		self.lb,self.ub = self.prb.get_bounds()

		self.ctype = ctype
		x0 = self.startp()
		if ctype == 'z':
			self.m = 0
			self.m = self.fconstriq(x0).shape[0]
		else:
		    self.m = self.fconstriq(x0).shape[0]

		self.name  = 'Problem DMS '+str(n)+' '+ctype+' '+self.prb.name
		#print('Hit RETURN to continue')
		#input()

		self.eps = 0.1*np.ones((self.q,self.m))

		if n < 4:
			return
		fob = self.functs(x0)
		ciq = self.fconstriq(x0)
		for k in range(self.q):
			for i in range(self.m):
				if np.maximum(0.0,ciq[i]) < 1.0:
					self.eps[k,i] = 1.e-3
				else:
					self.eps[k,i] = 1.e-1
		return

	def fconstreq(self,x):
		return np.empty(0)

	def fconstriq(self,xmix):
		x = np.copy(xmix)
		(lb,ub) = self.get_bounds()
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		if self.ctype == 'a':
			return fconstr_a(x)
		if self.ctype == 'b':
			return fconstr_b(x)
		if self.ctype == 'c':
			return fconstr_c(x)
		if self.ctype == 'd':
			return fconstr_d(x)
		if self.ctype == 'e':
			return fconstr_e(x)
		if self.ctype == 'f':
			return fconstr_f(x)
		if self.ctype == 'z':
			return fconstr_z(x)

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
		if self.nint > 0:
			x[self.ncont:] = self.lb[self.ncont:] + xmix[self.ncont:]*(self.ub[self.ncont:]-self.lb[self.ncont:])/(ub[self.ncont:]-lb[self.ncont:])

		f = self.prb.functs(x)
		self.FRONT.update_point(f)
		#print('F:',self.FRONT.F)
		FF = self.FRONT.get_clean_front(self.ref_point)
		ind = np.where(np.isnan(FF))
		try:
			FF = np.delete(FF,ind,axis=0)
			FT = transform_front(FF,self.idl_point,self.ref_point,self.q)
			nadir = transform_point(self.ref_point,self.idl_point,self.ref_point,self.q)
			if len(FT) > 0:
				hv = pg.hypervolume(FT)
				#print('R:',self.ref_point)
				self.HV.append(hv.compute(nadir))
			else:
				self.HV.append(0.0)
		except:
			self.HV.append(0.0)
		#print(self.HV[-1])
		#input()

		return f

	def get_bounds(self):
		lbcont = self.lb[:self.ncont]   #(self.ncont)
		ubcont = self.ub[:self.ncont]
		lbint = np.zeros(self.nint)
		ubint = 100*np.ones(self.nint)

		return (np.hstack((lbcont,lbint)),np.hstack((ubcont,ubint)))

	def get_nic(self):
		return 0 #self.m
	def get_nec(self):
		return 0
	def get_nobj(self):
		return self.q
	def get_nix(self):
		return self.nint

	def fitness(self,x):
		objs = self.functs(x)
		ci   = self.fconstriq(x)
		#return np.hstack((objs,ci))
		fmax = np.zeros(self.q)
		for i in range(self.m):
			for j in range(self.q):
				fmax[j] += np.maximum(0.0,ci[i]/self.eps[j,i])

		fpen = objs + fmax
		return fpen

class problem_factory():
	def __init__(self,prob_id=1,n=10,ctype='z',nint=-1,prbsel='deb512c'):
		#print(prbsel)
		if prob_id == 1:
			self.prob = problem_1()
		elif prob_id == 2:
			self.prob = problem_2(nint=nint)
		elif prob_id == 3:
			self.prob = problem_2_discretizzato()
		elif prob_id == 99:
			self.prob = problem_urgenze()
		elif prob_id == 108:
			self.prob = problem_solar8()
		elif prob_id == 109:
			self.prob = problem_solar9()
		elif prob_id == 1001:
			self.prob = problem_cec09_1(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1002:
			self.prob = problem_cec09_2(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1003:
			self.prob = problem_cec09_3(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1004:
			self.prob = problem_cec09_4(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1005:
			self.prob = problem_cec09_5(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1006:
			self.prob = problem_cec09_6(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1007:
			self.prob = problem_cec09_7(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1008:
			self.prob = problem_cec09_8(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1009:
			self.prob = problem_cec09_9(n=n,ctype=ctype,nint=nint)
		elif prob_id == 1010:
			self.prob = problem_cec09_10(n=n,ctype=ctype,nint=nint)
		elif prob_id == 2000:
			#print(prbsel)
			self.prob = problem_luis(prbsel=prbsel,ctype=ctype,nint=nint)

	def prob_ids(self):
		return [1,2,3,99,108,109, \
			    1001,1002,1003,1004,1005,1006,1007,1008,1009,1010, \
				2000]

	def describe(self):
		ids = self.prob_ids()

		for prob_id in ids:
			print('prob.id = ',prob_id,': ',end='')
			if prob_id == 1:
				print(problem_1().name)
			elif prob_id == 2:
				print(problem_2().name)
			elif prob_id == 3:
				print(problem_2_discretizzato().name)
			elif prob_id == 99:
				print(problem_urgenze().name)
			elif prob_id == 108: # SOLAR 8 problem
				print(problem_solar8().name)
			elif prob_id == 109: # SOLAR 9 problem
				print(problem_solar9().name)
			elif prob_id == 1001:
				print(problem_cec09_1().name)
			elif prob_id == 1002:
				print(problem_cec09_2().name)
			elif prob_id == 1003:
				print(problem_cec09_3().name)
			elif prob_id == 1004:
				print(problem_cec09_4().name)
			elif prob_id == 1005:
				print(problem_cec09_5().name)
			elif prob_id == 1006:
				print(problem_cec09_6().name)
			elif prob_id == 1007:
				print(problem_cec09_7().name)
			elif prob_id == 1008:
				print(problem_cec09_8().name)
			elif prob_id == 1009:
				print(problem_cec09_9().name)
			elif prob_id == 1010:
				print(problem_cec09_10().name)
			elif prob_id == 2000:
				print(problem_luis().name)
