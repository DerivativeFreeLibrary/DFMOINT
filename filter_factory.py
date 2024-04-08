import numpy as np
from utility import minore_uguale, domina
from filter_type import filter_elem

class front_factory():
	def __init__(self,q,Fin=None):
		self.q = q
		self.F = np.empty((0,q))
		if isinstance(Fin, np.ndarray):
			self.F = np.copy(Fin)
	def get_clean_front(self,ref_point):
		Ltemp = np.empty((0,self.q))
		for i in range(len(self.F)):
			if domina(self.F[i,:],ref_point):
				Ltemp = np.append(Ltemp,[self.F[i,:]],axis=0)
		return Ltemp
	def update_point(self,f):
		'''
		possono succedere TRE cose:
		1) Lkappa[i] < della tuple (x,f)
		2) (x,f) < Lkappa
		3) (x,f) non domina e non è dominato da Lkappa[i]
		'''
		flag = 0
		Ltemp = np.zeros((0,self.q))
		while flag < 3:
			if len(self.F) == 0:
				flag = 3
				continue
			for i in range(len(self.F)):
				#print(i,' ',f,' ',e.fobs)
				#print('A',self.F)
				#print('B',self.F[i,:])
				if minore_uguale(self.F[i,:],f) and minore_uguale(f,self.F[i,:]):
					flag = 1
					return False
				if(domina(self.F[i,:],f)):
					flag = 1
					return False
				elif(domina(f,self.F[i,:])):
					flag = 2
					self.F = np.delete(self.F,i,axis=0)
					break
				else:
					flag = 3

		self.F = np.append(self.F,[f],axis=0)
		return True

class filter_factory():
	def __init__(self):
		self.Lkappa = np.empty(0,dtype=filter_elem)
		self.Lnew   = np.empty(0,dtype=filter_elem)
		self.Ltilde = np.empty(0,dtype=filter_elem)
		self.iprint = 0

	def sizes(self):
		return self.Lkappa.shape[0], self.Ltilde.shape[0], self.Lnew.shape[0]

	def init_Lkappa(self):
		self.Lkappa = np.empty(0,dtype=filter_elem)

	def init_Lnew(self):
		self.Lnew = np.empty(0,dtype=filter_elem)

	def init_Ltilde(self):
		self.Ltilde = np.empty(0,dtype=filter_elem)

	def intorno(self,x,f,delta):
		theta = 1.0

		''' iterate through Lkappa '''
		for i, e in enumerate(self.Lkappa):
			if np.all(np.abs(e.fobs-f) < np.minimum(1.0,theta*delta)):
				return False

		''' iterate through Lnew '''
		for i, e in enumerate(self.Lnew):
			if np.all(np.abs(e.fobs-f) < np.minimum(1.0,theta*delta)):
				return False

		''' iterate through Ltilde '''
		for i, e in enumerate(self.Ltilde):
			if np.all(np.abs(e.fobs-f) < np.minimum(1.0,theta*delta)):
				return False

		return True

	def add_point(self,x,f,alfa,alfac,alfa_tilde,xi,allones):
		''' add tuple dict_elem = (x,f,alfa,alfac) to Lnew'''

		self.Lnew = np.append(self.Lnew,filter_elem(x,f,alfa,alfac,alfa_tilde,xi,allones))
		if self.iprint > 0:
			''' print out Lnew '''
			for i, e in enumerate(self.Lnew):
				print('Lnw:',i,end='')
				e.print(flagx=False)

	def add_point1(self,e):
		''' add tuple dict_elem = (x,f,alfa,alfac) to Lnew'''

		self.Lnew = np.append(self.Lnew,e)

		if self.iprint > 0:
			''' print out Lnew '''
			for i, e in enumerate(self.Lnew):
				print('Lnw:',i,end='')
				e.print(flagx=False)


	def update_point(self,x,f,alfa,alfa_c,alfa_tilde,xi,allones):
		'''
		possono succedere TRE cose:
		1) Lkappa[i] < della tuple (x,f)
		2) (x,f) < Lkappa
		3) (x,f) non domina e non è dominato da Lkappa[i]
		'''
		flag = 0
		Ltemp = np.empty(0,dtype=filter_elem)
		while flag < 3:
			if len(self.Lkappa) == 0:
				flag = 3
				continue
			for i, e in enumerate(self.Lkappa):
				#print(i,' ',f,' ',e.fobs)
				if minore_uguale(e.fobs,f) and minore_uguale(f,e.fobs):
					flag = 1
					return False
				if(domina(e.fobs,f)):
					flag = 1
					return False
				elif(domina(f,e.fobs)):
					flag = 2
					self.Lkappa = np.delete(self.Lkappa,i)
					break
				else:
					flag = 3

		self.Lkappa = np.append(self.Lkappa,filter_elem(x,f,alfa,alfa_c,alfa_tilde,xi,allones))
		return True

	def merge_sets(self):
		'''
		questa subroutine fa il merge delle liste Lnew e Lkappa
		'''
		flag = False

		for i,e in enumerate(self.Lnew):
			flag = flag or update_point(self,e.x,e.fobs,e.alfa,e.alfa_c,e.alfa_tilde,e.xi,e.flag_allones)
		return flag

	def domina(self,f,rho):
		for i,e in enumerate(self.Lkappa):
			if domina(e.fobs-rho,f):
				return True
		return False

	def print(self,prob):
		print('------------- print filter Lnew -----------------')
		totfeas = 0
		for i,e in enumerate(self.Lnew):
			print(' %7d' % i,end='')
			#print(' x=',end='')
			#for j,x in enumerate(e['x']):
			#	print(' %20.13e' % x,end='')
			print(' f=',end='')
			for j,f in enumerate(e.fobs):
				print(' %20.13e' % f,end='')
			print(' alfa= %13.6e' % e.alfa,end='')
			for j,alfa in enumerate(e.alfa_c):
				print(' %13.6e' % alfa,end='')

			ciq = prob.fconstriq(e.x)
			viol = 0.0
			if len(ciq) > 0:
				viol = np.maximum(0.0,np.max(ciq))
			if viol < 1.e-3:
				totfeas += 1
			print(' viol= %20.13e' % viol)
		print('totfeas = ',totfeas)
		print('------------- END print filter Lnew -------------')

	def print_filter(self,prob):
		'''
		prints out the filter in Lkappa
		'''
		print('------------- print filter Lkappa -----------------')
		totfeas = 0
		for i,e in enumerate(self.Lkappa):
			print(' %7d' % i,end='')
			#print(' x=',end='')
			#for j,x in enumerate(e['x']):
			#	print(' %20.13e' % x,end='')
			print(' f=',end='')
			for j,f in enumerate(e.fobs):
				print(' %20.13e' % f,end='')
			#print(' alfa= %20.13e' % e['alfa'],end='')
			#for j,alfa in enumerate(e['alfa_c']):
			#	print(' %20.13e' % alfa,end='')

			ciq = prob.fconstriq(e.x)
			viol = 0.0
			if len(ciq) > 0:
				viol = np.maximum(0.0,np.max(ciq))
			if viol < 1.e-3:
				totfeas += 1
			print(' viol= %20.13e' % viol)
		print('totfeas = ',totfeas)
		print('------------- END print filter Lkappa -------------')
