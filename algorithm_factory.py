import numpy as np
from filter_type import filter_elem
from filter_factory import filter_factory
from cache_factory import cache_factory
from utility import maggiore_stretto, minore_uguale, merge
import ghalton
import sobol_seq


class dfmoint():
	'''
	Classe del metodo di ottimizzazione per problemi multi obiettivo
	misti interi.
	I metodi esportati dalla classe sono:
	* __init__
	* run   (è il metodo main)
	* gen_base
	* gram_schmidt
	* spread_ord
	* functs_pen
	* diagonale
	* linesearchbox_cont
	* linesearchbox_dense
	* discrete_linesearch
	* stop
	* halton
	'''

	def __init__(self,problem=[],n=1,ctype='z'):
		self.prob = problem
		self.FILTER = filter_factory()
		(n,ncont,nint,m,q) = self.prob.getdim()
		print('solving problem ',self.prob.name)
		print('algo init: ',n,' ',ncont,' ',nint,' ',m,' ',q)
		self.CACHE = cache_factory(n,m,q)

	def gen_base(self,d):
	    n = len(d)
	    ind = np.argmax(np.abs(d))
	    H = np.zeros((n,n))
	    H[:,0] = d
	    #print(ind)
	    for i in range(1,ind+1):
	        H[i-1,i] = np.double(1.0)
	    for i in range(ind+1,n):
	        H[i,i] = np.double(1.0)

	    #print('================================================')
	    #print('ind = ',ind,'d_dense = ')
	    #for j in range(n):
	    #	print('%15.6e ' % d[j],end='')
	    #print()
	    #for i in range(n):
	    #	for j in range(n):
	    #		print('%15.6e ' % H[i,j],end='')
	    #	print()
	    #input()
	    return H

	def gram_schmidt(self,H):
		(n,n) = H.shape
		for i in range(1,n):
		    proj = np.double(0.0)
		    for j in range(i):
		        proj += (np.dot(H[:,i],H[:,j])/np.dot(H[:,j],H[:,j]))*H[:,j]
		    H[:,i] -= proj

		for i in range(n):
		    H[:,i] = H[:,i]/np.linalg.norm(H[:,i])

		for i in range(n):
			for j in range(n):
				if(np.abs(H[i,j]) < 1.e-4):
					H[i,j] = 0.0

		for i in range(n):
		    H[:,i] = H[:,i]/np.linalg.norm(H[:,i])

		return H

	def spread_ord(self):
		#print('---------------------------------------------------')
		#print('------------------  in SPREAD  --------------------')
		#print('---------------------------------------------------')

		(n,ncont,nint,m,q) = self.prob.getdim()
		ncomp  = 0
		Lktemp = np.empty(0,dtype=filter_elem)

		ndim = len(self.FILTER.Lkappa)
		delta = np.zeros(ndim)
		fobj = np.zeros((ndim,q))
		indf = [i for i in range(ndim)]
		for i,e in enumerate(self.FILTER.Lkappa):
			ncomp += 1
			fobj[i,:] = e.fobs

		if ncomp <= 1:
			#print('---------------------------------------------------')
			#print('------------------ out SPREAD  --------------------')
			#print('---------------------------------------------------')
			return self.FILTER.Lkappa,False,delta

		app  = np.zeros(ndim)
		iapp = np.zeros(ndim,dtype=int)
		for j in range(q):
			indf = [i for i in range(ndim)]
			alfas = fobj[:,j]
			ipt = np.argsort(alfas)
			for i in range(ndim):
				app[ndim-(i+1)] = alfas[ipt[i]]
				iapp[ndim-(i+1)]= indf[ipt[i]]


			#print('spread_ord: alfas=',alfas)
			#print('spread_ord: app=',app)
			alfas = app.copy()
			indf = iapp.copy()

			#for i in range(ndim):
			#	print(alfas[i],self.FILTER.Lkappa[indf[i]].fobs[j])

			#print('spread_ord: ipt=',ipt)
			#print('spread_ord: indf=',indf)
			#print('spread_ord: alfas=',alfas)

			delta[indf[0]] += 1.e+10 + np.abs(alfas[0])
			delta[indf[-1]] += 1.e+9 + np.abs(alfas[-1])
			for i in range(1,ndim-1):
			#	print('i=',i,' indf[i]=',indf[i],' delta[indf[i]]=',delta[indf[i]])
			#	print('a1=',(alfas[i-1]-alfas[i+1]),' a2=',(alfas[0]-alfas[-1]))
				delta[indf[i]] += np.abs((alfas[i-1]-alfas[i+1])/(alfas[0]-alfas[-1]))

			#for i in range(ndim):
			#	print('spread_ord: delta[indf[i]]=',delta[indf[i]])

		delta /= ndim
		#print('spread_ord: delta=',delta)
		ipt = np.argsort(delta)
		for i in range(ndim):
			app[ndim-(i+1)] = delta[ipt[i]]
			Lktemp = np.insert(Lktemp,0,self.FILTER.Lkappa[ipt[i]])

		#print('spread_ord:')
		#for i,e in enumerate(Lktemp):
		#	print(i,Lktemp[i].fobs)
		#print(app)
		#print('---------------------------------------------------')
		#print('------------------ out SPREAD  --------------------')
		#print('---------------------------------------------------')

		#input()
		return Lktemp,True,app

	def functs_pen(self,fob,ciq,eps):
		(n,ncont,nint,m,q) = self.prob.getdim()
		fmax = np.zeros(q)
		viol = 0.0
		for i in range(m):
			for j in range(q):
				fmax[j] += np.maximum(0.0,ciq[i]/eps[j,i])
			viol += np.maximum(0.0,ciq[i])

		fpen = fob + fmax

		return fpen, viol

	def diagonale(self):
		(bl,bu) = self.prob.get_bounds()
		(n,ncont,nint,m,q) = self.prob.getdim()

		punti = np.zeros((n,n))

		if n > 1:
			for i in range(n):
				punti[:ncont,i] =         bl[:ncont] + (bu[:ncont]-bl[:ncont])*i/(n-1)
				punti[ncont:,i] = np.ceil(bl[ncont:] + (bu[ncont:]-bl[ncont:])*i/(n-1))
		else:
			punti[0,0] = bl[0]

		return punti

	def linesearchbox_cont(self,x,f,d,alfa,alfa_d,alfa_dense,alfa_tilde,z,fz,i_corr,alfa_max,iprint,bl,bu,nf,pn,ifront,eps,xi):
		(n,ncont,nint,m,q) = self.prob.getdim()
		gamma  = 1.e-6
		delta  = 0.5
		delta1 = 0.5
		alfac  = np.copy(alfa_d)
		j      = i_corr
		zdelta = np.copy(x)
		if iprint > 1:
			print('variabile continua  j =',j,'    d(j) =',d[j],' alfa=',alfa_d[j])

		while True:
			if ifront == 1:
				if iprint > 1:
					print(' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa)
				alfa_d[j] = alfa
				self.FILTER.add_point(z,fz,alfa_dense,alfa_d,alfa_tilde,xi,False)
				return alfa, z, nf, ifront

			if d[j] > 0.0:
				if alfa/delta1 - (bu[j]-x[j]) < - 1.e-6:
					alfaex = alfa/delta1
				else:
					alfaex = bu[j]-x[j]
					ifront = 1
					if iprint > 1:
						print(' punto espan. sulla front.')
			else:
				if alfa/delta1 - (x[j]-bl[j]) < - 1.e-6:
					alfaex = alfa/delta1
				else:
					alfaex = x[j]-bl[j]
					ifront = 1
					if iprint > 1:
						print(' punto espan. sulla front.')

			zdelta[j] = x[j]+alfaex*d[j]

			found_incache,ind,fob,ciq = self.CACHE.find(zdelta)
			if not found_incache:
				if self.prob.name == "Problem SOLAR 8":
					fob,ciq = self.prob.run_bb(zdelta)
				else:
					fob = self.prob.functs(zdelta)
					ciq = self.prob.fconstriq(zdelta)
				#fob = self.prob.functs(zdelta)
				#ciq = self.prob.fconstriq(zdelta)
				self.CACHE.insert(zdelta,fob,ciq)
				nf += 1
			fzdelta,violzdelta = self.functs_pen(fob,ciq,eps)

			if iprint > 1:
				print('n,q,m=',n,q,m)
				print('  fob=',fob)
				print('  ciq=',ciq)
				print(' fzex=',fzdelta,'  alfaex=',alfaex)
			if iprint >= 2:
				for i in range(n):
					print(' z(',i,')=',z[i],zdelta[i])

			fpar = f - gamma*alfaex*alfaex
			fpar2= fz-gamma*(alfaex-alfa)*(alfaex-alfa)
			fpar2= fz-gamma*alfaex*alfaex+gamma*alfa*alfa
			if not minore_uguale(fzdelta,fpar2):
				if iprint > 1:
					print(' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa)

				alfa_d[j] = alfa
				self.FILTER.add_point(z,fz,alfa_dense,alfa_d,alfa_tilde,xi,False)

			#self.FILTER.print_filter(self.prob)
			#input()
			if self.FILTER.domina(fzdelta,gamma*alfaex*alfaex):
				if iprint > 1:
					print('Lkappa_domina')
				return alfa, z, nf, ifront

			else:
				if iprint > 1:
					print('fzdelta NOT > fpar')
				if not (ifront == 1):
					fz = fzdelta
					z  = np.copy(zdelta)
					violz = violzdelta
					alfa = alfaex
				else:
					if iprint > 1:
						print(' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa)

					alfa_d[j] = alfaex
					self.FILTER.add_point(zdelta,fzdelta,alfa_dense,alfa_d,alfa_tilde,xi,False)

		return alfa, z, nf, ifront

	def linesearchbox_dense(self,x,f,d,      alfa,alfac, alfa_dense,alfa_tilde,z,fz,alfa_max,iprint,bl,bu,nf,pn,ifront,eps,xi):
		(n,ncont,nint,m,q) = self.prob.getdim()
		gamma  = 1.e-6
		delta  = 0.5
		delta1 = 0.5
		ifront = 0

		if iprint > 0:
			print('direzione halton, alfa=',alfa_dense)

		alfa = alfa_dense
		alfaex = alfa
		#############################################
		# PROJECTED STRONG EXPANSION
		#############################################

		while True:
			alfaex = alfa/delta1
			zdelta = x + alfaex*d
			zz     = np.maximum(bl,np.minimum(bu,x+alfa*d))
			zdelta = np.maximum(bl,np.minimum(bu,zdelta))
			if np.linalg.norm(zdelta-zz) <= 1.e-8:
				d *= -1
				alfa = alfa/2.0
				#print('zdelta vicino zz')
				break

			found_incache,ind,fob,ciq = self.CACHE.find(zdelta)
			if not found_incache:
				if self.prob.name == "Problem SOLAR 8":
					fob,ciq = self.prob.run_bb(zdelta)
				else:
					fob = self.prob.functs(zdelta)
					ciq = self.prob.fconstriq(zdelta)
				#fob = self.prob.functs(zdelta)
				#ciq = self.prob.fconstriq(zdelta)
				nf += 1
				self.CACHE.insert(zdelta,fob,ciq)
			fzdelta, violzdelta = self.functs_pen(fob,ciq,eps)
			if iprint > 0:
				print(' alfaex = ',alfaex)
				print(' d = ',d)
				print(' x = ',x)
				print(' z = ',zdelta)
				print(' fzex=',fzdelta,'  alfaex=',alfaex)
			if iprint >= 2:
				for i in range(n):
					print(' z(',i,')=',z[i])

			fpar  = f  - gamma*alfaex*alfaex
			fpar2 = fz - gamma*(alfaex-alfa)*(alfaex-alfa)
			fpar2 = fz - gamma*alfaex*alfaex + gamma*alfa*alfa

			if iprint > 0:
				print(' fzdelta =',fzdelta,'   fpar2 =',fpar2)

			if not minore_uguale(fzdelta,fpar2):
				if iprint > 0:
					print(' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa)
				self.FILTER.add_point(z,fz,alfa,alfac,alfa_tilde,xi,False)

			if self.FILTER.domina(fzdelta,gamma*alfaex*alfaex):
				if iprint > 0:
					print('Lkappa_domina')
				return alfa, z, nf,ifront
			else:
				if iprint > 0:
					print('fzdelta NOT > fpar')
				fz = fzdelta
				z  = np.copy(zdelta)
				violz = violzdelta
				alfa = alfaex

		return alfa, z, nf, ifront

	def discrete_linesearch(self,y,d,alfa_d,alfa_dense,alpha_tilde,ii,  lb,ub,f_ref,eps,xi,m,ncont,nf,outlev):
		#	 			   (x,d,alfa_d,alfa_dense,alpha_tilde,idir,lb,ub,f,    eps,xi,m,ncont,nf,iprint)
		#
		# Function discrete_linesearch
		#
		# Purpose:
		#
		# This function performs a discrete linesearch along a primitive direction
		# d (d \in Z^n). It possibly updates the list Lnew.
		# ATTENZIONE: fallimento vuol dire che non è stato aggiunto alcun punto
		#             a Lnew
		#
		# Inputs:
		#
		# y            : starting point for the linesearch
		#
		# d            : search direction
		#
		# alpha_tilde  : starting stepsize
		#
		# lb, ub       : lower and upper bounds
		#
		# f_ref        : reference o.f. value
		#
		# Output:
		#
		#
		# alpha        : 1) alpha > 0 if linesearch finds a point guaranteeing
		#                simple decrease: f(y+alpha d)<f_ref
		#                2) alpha = 0 failure
		#
		# nf           : number of function avaluations
		#

		# calculate dimension of the problem
		n = len(d)

		# initialize vector alpha_max
		alpha_max = np.inf * np.ones(n)

		# calculate max alpha
		indices = ( d > 0 )

		alpha_max[indices]=np.divide(ub[indices] - y[indices],d[indices])

		indices = ( d < 0 )

		alpha_max[indices]=np.divide(lb[indices] - y[indices],d[indices])

		#compute starting alpha
		alpha_bar  = np.floor( min(alpha_max) )
		alpha_init = min(alpha_tilde[ii], alpha_bar)

		if outlev >= 2:
		    print('discrete_search: alpha_init = ',alpha_init)

		#Build first point for starting linesearch
		if (alpha_init > 0):
			y_trial = y + alpha_init * d
			found_incache,ind,fob,ciq = self.CACHE.find(y_trial)
			if not found_incache:
				if self.prob.name == "Problem SOLAR 8":
					fob,ciq = self.prob.run_bb(y_trial)
				else:
					fob = self.prob.functs(y_trial)
					ciq = self.prob.fconstriq(y_trial)
				#fob = self.prob.functs(y_trial)
				#ciq = self.prob.fconstriq(y_trial)
				nf += 1
				self.CACHE.insert(y_trial,fob,ciq)
			f_trial, viol_trial = self.functs_pen(fob,ciq,eps)
		else:
			f_trial = np.inf


	    # cicle for updating alpha
		if outlev >= 1:
			print('discrete_search: ftrial = ',f_trial)
			print('               : f_ref  = ',f_ref)
			print('               : alpha_init = ',alpha_init)

		# if (alpha_init > 0) and (f_trial <= f_ref - xi):

		#if (alpha_init > 0) and not maggiore_stretto(f_trial, f_ref - xi):
		if (alpha_init > 0) and not self.FILTER.domina(f_trial,xi):

			# initialize alpha and best point
			if outlev >= 1:
				print('discrete_search: f_trial not > f_ref - xi')
			alpha=alpha_init
			x = y_trial
			f = f_trial

			#calculate trial point
			if outlev >= 1:
				print('discrete_search: alpha=',alpha,' alphabar=',alpha_bar)
			if alpha < alpha_bar:
				y_trial = y + min(alpha_bar,2*alpha)* d
				found_incache,ind,fob,ciq = self.CACHE.find(y_trial)
				if not found_incache:
					if self.prob.name == "Problem SOLAR 8":
						fob,ciq = self.prob.run_bb(y_trial)
					else:
						fob = self.prob.functs(y_trial)
						ciq = self.prob.fconstriq(y_trial)
					#fob = self.prob.functs(y_trial)
					#ciq = self.prob.fconstriq(y_trial)
					nf += 1
					self.CACHE.insert(y_trial,fob,ciq)
				f_trial, viol_trial = self.functs_pen(fob,ciq,eps)
			else:
				f_trial = np.inf


			# expansion step (increase stepsize)
			#while (alpha<alpha_bar) and (f_trial <= f_ref - xi):
			while (alpha<alpha_bar) and not self.FILTER.domina(f_trial,xi):
			#while (alpha<alpha_bar) and not maggiore_stretto(f_trial, f - xi):
				# se alpha < alpha_bar e
				# NON ESISTE un punto nel filtro che DOMINA sufficientemente (con xi) f_trial
				# allora ...
				if outlev >= 1:
					print('discrete_search: NOT EXISTS e in Lkappa: e.fobs - xi "domina" f_trial (while loop)')
				# alpha calulation and best point updatingd

				if not minore_uguale(f_trial,f - xi):
					if outlev >= 1:
						print(' discrete_search: f_trial NOT <= f - xi')
						print(' discrete_search: AGGIUNGE PUNTo Lnew fz =',f,'   alfa =',alpha)

					alpha_tilde[ii] = alpha
					self.FILTER.add_point(x,f,alfa_dense,alfa_d,alpha_tilde,xi,False)

				alpha=min(alpha_bar, 2*alpha)

				# best point updating
				x = y_trial
				f = f_trial

				#next point to be tested
				if outlev >= 1:
					print('discrete_search(1): alpha=',alpha,' alphabar=',alpha_bar)
				if(alpha < alpha_bar):
					y_trial = y + min(alpha_bar, 2* alpha) * d
					found_incache,ind,fob,ciq = self.CACHE.find(y_trial)
					if not found_incache:
						if self.prob.name == "Problem SOLAR 8":
							fob,ciq = self.prob.run_bb(y_trial)
						else:
							fob = self.prob.functs(y_trial)
							ciq = self.prob.fconstriq(y_trial)
						#fob = self.prob.functs(y_trial)
						#ciq = self.prob.fconstriq(y_trial)
						nf += 1
						self.CACHE.insert(y_trial,fob,ciq)
					f_trial, viol_trial = self.functs_pen(fob,ciq,eps)
				else:
					f_trial = np.inf

			if not self.FILTER.domina(f,xi):
				if outlev >= 1:
					print('discrete_search: NOT EXISTS e in Lkappa: e.fobs - xi "domina" f (after while loop)')
					print(' AGGIUNGE PUNTo Lnew fz =',f,'   alfa =',alpha)

				alpha_tilde[ii] = alpha
				self.FILTER.add_point(x,f,alfa_dense,alfa_d,alpha_tilde,xi,False)

		else:
			###########################################
			alpha = 0

		return alpha, nf
##########################################################################
# END OF CODE discrete_linesearch
##########################################################################

	def stop(self,alfa_d,alfa_dense,alfa_tilde,alfa_max,nf,f,alfa_stop,nf_max):
		(n,ncont,nint,m,q) = self.prob.getdim()
		istop = 0
		alfa_max = 0.0

		if n > 1:
			alfa_max = np.maximum(alfa_max,alfa_dense)

		if alfa_max <= alfa_stop:
			istop = 1

		if nf > nf_max:
			istop = 2

		return istop

	def halton(self,index):
		(n,ncont,nint,m,q) = self.prob.getdim()
		n = ncont
		primes = np.load('PRIMES.sav.npy')
		x = np.zeros(n)
		if n==1:
			x[0] = 1.0
			return x
		base = np.zeros(n)
		for i in range(n):
			base[i] = primes[i]
			x[i] = 0.0
			f = 1.0/(base[i])
			j = index
			while j > 0:
				x[i] += f*np.mod(j,base[i])
				#print('j = ',j)
				j //= np.int(base[i])
				#print('j = ',j)
				f /= base[i]
				#print('f = ',f)
		#print('halton: ',index)
		#print('halton: ',x)
		#input()
		return x

	def prime_vector(self,d):
	    n = len(d)
	    flag = 0
	    if(n==1):
	        flag = True
	        return flag
	    temp = np.gcd(np.array(abs(d[0]),dtype=int),np.array(abs(d[1]),dtype=int))
	    if(n==2):
	        flag = (temp == 1)
	        return flag
	    for i in np.arange(2,n,1):
	        temp = np.gcd(temp,np.array(abs(d[i]),dtype=int))
	        #temp = numpy_gcd(temp,abs(d[i]))
	        if temp == 1:
	            flag = True
	            return flag
	    if temp != 1:
	        flag = False
	        return flag

	##########################################################################
	# END OF CODE prime_vector
	##########################################################################

	def generate_dirs(self,ncont,n,D,eta,betaLS,Phalton,ihalton):
	    #
	    # Function generate_dirs
	    #
	    # Purpose:
	    #
	    # This function generate new integer directions which are added to set D
	    #
	    # Inputs:
	    #
	    # n            : dimension of the problem
	    #
	    # D            : matrix of current directions (one per each column)
	    #
	    # alpha_tilde  : array of stepsizes along direction already in D
	    #
	    # Output:
	    #
	    # Dout         : [new_direction D]
	    #
	    # succout      : [0 succ]
	    #
	    # alpha        : array of stepsizes along the directions in Dout
	    #                alpha = [new_step_sizes alpha_tilde]
	    #

	    mD = np.shape(D)[1]

	    for j in range(1000):
	        #keyboard
		    if True:
		        v = 2*np.asarray(Phalton[ihalton-1], dtype = np.float64) - np.ones(n)
		        ihalton += 1
		        v = eta*(v/np.linalg.norm(v))

		        if (np.linalg.norm(v) < 1e-16):
		            break
		    else:

		        if(n > 1):
		            v = 2*np.asarray(Phalton[ihalton-1], dtype = np.float64) - np.ones(n)
		            v = 10.0*np.random.rand()*v
		            ihalton += 1
		        else:
		            v = np.array([1.0])

	        #d = abs(round(v)) good if H=norm(d)^2*eye(n,n) - 2*d*d' used
		    d = np.round(v)

	        #now check whether d is a prime vector
		    #if True:
		    if self.prime_vector(d) == True:
		        trovato = False
		        #check whether d is already in D
		        d = np.reshape(d,(len(d),1))
		        d = np.concatenate((np.zeros((ncont,1)),d),axis=0)
		        DIFF1 = D - np.tile(d,(1,mD))
		        DIFF2 = D + np.tile(d,(1,mD))
		        if( min ( np.sum(abs(DIFF1),axis=0)) == 0 ) or ( min ( np.sum(abs(DIFF2),axis=0)) == 0 ):
		            trovato = True

		        if trovato == False:
		            #print(d)
		            #input()
		            H       = d.copy() #norm(d)^2*eye(n,n) - 2*d*d'
		            Dout    = np.hstack((H,D))
		            iexit   = 1
		            return Dout, iexit, ihalton

	    Dout    = D
	    iexit   = 0

	    return Dout, iexit, ihalton

	##########################################################################
	# END OF CODE generate_dirs
	##########################################################################

	def run(self,alfa_stop=1.e-6,iprint=0,nf_max=2000,print_all_filters=False,versione_FORTRAN=False):
		(bl,bu) = self.prob.get_bounds()
		#print(bl,bu)
		#input()
		(n,ncont,nint,m,q) = self.prob.getdim()

		eps = 0.1*np.ones((q,m))

		alfainiz = np.minimum(10.0,np.max(bu-bl)/10.0)
		alfaciniz = np.zeros(n)
		for i in range(ncont):
			alfaciniz[i] = np.minimum(10.0,(bu[i]-bl[i])/10.0)
		xi = 1.e-3 #0.5 #1.0

		print_line = 0
		nf = 0
		coef_delta = 1.0
		eta = 1.5
		gamma=1.e-6
		num_gen = 0
		num_merge = 0
		num_used = 0
		istop = 0
		cambio_eps = False
		alfa_dense = 0.0
		index_halton = 1000 + 2*n
		index_sobol  = 10000
		index_primitive = 7
		hschoice = 2
		''' d_coord è la matrice identità di dimensione n x n'''
		d_coord = np.eye(n)

		'''
		D è la matrice che contiene (inizialmente) le direzioni coord delle
		discrete. Poi, viene arricchita con le direzioni primitive
		'''
		if nint > 0:
			D = np.concatenate((np.zeros((ncont,nint)),np.identity(nint)),axis=0)
			ndir          = np.shape(D)[1]
			alfa_tilde    = np.round((bu[ncont:]-bl[ncont:])/2.0)
			for i in range(nint):
				alfa_tilde[i] = np.maximum(1.0,np.round(np.minimum(10.0,(bu[ncont+i]-bl[ncont+i])/10.0)))
		else:
			alfa_tilde    = np.array([])
			ndir = 0
		old_maxalpha  = np.inf

		print(bu)
		print(bl)
		print(alfa_tilde)
		#input()
		alfa_d = np.zeros(n)
		sequencer   = ghalton.Halton(ncont)       ## (input) la dimensione dello spazio
		Phalton     = sequencer.get(1000000)   ## (input) il numero di punti da generare nello spazio
		sequencer_int = ghalton.Halton(nint)
		Phalton_int   = sequencer_int.get(1000000)

		if ncont > 1:
			if hschoice == 1:
				d_dense = np.asarray(Phalton[index_halton-1], dtype = np.double)
				#d_dense = self.halton(index_halton)
			else:
				d_dense, index_sobol = sobol_seq.i4_sobol(ncont, index_sobol)
		else:
			d_dense = np.array([1.0])

		''' H è una matrice a blocchi. Il blocco in alto a sinistra ncont x ncont
		    è una matrice ortonormale di direzioni dense. Il blocco in basso
			a destra l'identità di dimensione nint x nint
		'''
		if ncont > 0:
			H = self.gen_base(d_dense)
			H = self.gram_schmidt(H)
			if nint > 0:
				H = np.vstack((H,np.zeros((nint,ncont))))
				H = np.hstack((H,np.vstack((np.zeros((ncont,nint)),np.eye(nint)))))
			if iprint >= 1:
				print(H)
				print('Hit return to continue...')
				input()
		if ncont == 0 and nint > 0:
			H = np.eye(nint)

		if iprint > 1:
			print('call startp ...')
		x = self.prob.startp()

		if iprint > 1:
			print('bl=',bl)
			print(' x=',x)
			print('bu=',bu)
			print('check box constraints ...')
			input()

		for i in range(n):
			if (x[i] < bl[i]) or (x[i] > bu[i]):
				print('The starting point violates bound constraints, STOP')
				return False

		if iprint > 1:
			print('compute fob + constr onto initial point')
		found_incache,ind,fob,ciq = self.CACHE.find(x)
		if not found_incache:
			if self.prob.name == "Problem SOLAR 8":
				fob,ciq = self.prob.run_bb(x)
			else:
				fob = self.prob.functs(x)
				ciq = self.prob.fconstriq(x)
			nf += 1
			self.CACHE.insert(x,fob,ciq)

		if iprint > 1:
			print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
			print('fob = ',fob)
			print('ciq = ',ciq)
			print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
			print('--------------------------------------')
			print('--------- EPSILON INIZIALI -----------')
			print('--------------------------------------')
		if self.prob.name == "Problem SOLAR 8":
			for k in range(q):
				for i in range(m):
					if np.maximum(0.0,ciq[i]) < 1.0:
						eps[k,i] = 1.e-3*self.prob.scale[q+i]
					else:
						eps[k,i] = 1.e-1*self.prob.scale[q+i]
		else:
			for k in range(q):
				for i in range(m):
					if np.maximum(0.0,ciq[i]) < 1.0:
						eps[k,i] = 1.e-3
					else:
						eps[k,i] = 1.e-1

		finiz, violiniz = self.functs_pen(fob,ciq,eps)
		self.FILTER.update_point(x,finiz,alfainiz,alfaciniz,alfa_tilde,xi,False)

		punti = self.diagonale()
		for i in range(n): #[n-1]: #range(n):
			x = punti[:,i]
			#for t in range(n):
			#	print(' %16.9e' % x[t],end='')
			#print()
			#input()
			found_incache,ind,fob,ciq = self.CACHE.find(x)
			if not found_incache:
				if self.prob.name == "Problem SOLAR 8":
					fob,ciq = self.prob.run_bb(x)
				else:
					fob = self.prob.functs(x)
					ciq = self.prob.fconstriq(x)
				#fob = self.prob.functs(x)
				#ciq = self.prob.fconstriq(x)
				nf += 1
				self.CACHE.insert(x,fob,ciq)
			finiz, violiniz = self.functs_pen(fob,ciq,eps)
			#if self.FILTER.intorno(x,finiz,1.e-2):
			self.FILTER.update_point(x,finiz,alfainiz, alfaciniz,alfa_tilde,xi,False)

		#self.FILTER.print_filter(self.prob)
		#input()
		########################################
		# INIZIO CICLO PRINCIPALE
		########################################
		#######   85     0     0  1506  4229 Lkappa
		print(' Lkap  Ltil  Lnew    nf ncach   type      ndir')
		print('-----------------------------------------------')
		while True:
			if print_all_filters:
				fid_Lkappa = open('Lkappa_'+str(num_merge)+'.txt','w')
				for i,e in enumerate(self.FILTER.Lkappa):
					if m >= 1:
						ciq = self.prob.fconstriq(e.x)
						viol = np.maximum(0.0,np.max(ciq))
					else:
						viol = 0.0
					if viol <= 1.e-3:
						print(' %6d' % i,end='',file=fid_Lkappa)
						for fob in e.fobs:
							print(' %13.6e' % fob,end='',file=fid_Lkappa)
						for i in range(ncont):
							print(' %10.3e' % e.x[i],end='',file=fid_Lkappa)
						for i in range(nint):
							print(' %3d' % e.x[ncont+i],end='',file=fid_Lkappa)
						print(file=fid_Lkappa)
				fid_Lkappa.close()

			for i in range(q):
				fbest = 1.e+100
				ibest = -1
				for pn,e in enumerate(self.FILTER.Lkappa):
					#print(i,' ',pn,' ',e)
					if fbest > e.fobs[i]:
						fbest = e.fobs[i]
						ibest = pn

				if ibest >= 0:
					e = self.FILTER.Lkappa[ibest]
					self.FILTER.Lkappa = np.delete(self.FILTER.Lkappa,ibest)
					self.FILTER.Lkappa = np.insert(self.FILTER.Lkappa,0,e)

			self.FILTER.Lkappa, flag_spread, deltaf = self.spread_ord()

			'''
			computes Lmaxalfa and Lmaxalfa2
			'''
			Lmaxalfa = 0.0
			Lmaxalfa2= 0.0
			#if num_used < 100:
			#	num_used = 0

			for i,e in enumerate(self.FILTER.Lkappa):
				num_used += 1
				if Lmaxalfa < e.alfa:
					Lmaxalfa = e.alfa
				if Lmaxalfa2 < e.alfa:
					Lmaxalfa2 = e.alfa
				Lmaxalfa2 = np.maximum(Lmaxalfa2,np.max(e.alfa_c))
			if iprint > 1:
				print('Lmaxalfa  = ',Lmaxalfa)
				print('Lmaxalfa2 = ',Lmaxalfa2)

			if versione_FORTRAN and (Lmaxalfa2 < alfa_stop):
				break

			if (Lmaxalfa2 < alfa_stop) and ( \
			    (nint == 0) or \
			    ((nint > 0) and ((index_primitive > 300000) or (ndir > 2000)))):
				break

			trovato = False
			#if num_used < 100:
			#	num_gen = 0
			if num_used > 0:
				num_gen += 1

			######################################################
			# Refine self.FILTER.Lkappa
			# exploring just the continuous variables
			######################################################
			for pnum,e in enumerate(self.FILTER.Lkappa):
				pn = pnum
				evolve_pop = False
				viol = 0.0
				if m >= 1:
					x = e.x
					found_incache,ind,fob,ciq = self.CACHE.find(x)
					if not found_incache:
						ciq = self.prob.fconstriq(x)
					viol = np.maximum(viol,np.max(ciq))

				alfa_max = np.maximum(e.alfa,np.max(e.alfa_c))
				if (deltaf[pn] == 0.0) and flag_spread:
					if iprint > 0:
						print(pn)

				#self.FILTER.print_filter(self.prob)
				if versione_FORTRAN:
					condizione = ((e.alfa > alfa_stop) and \
					              ((deltaf[pn] >= coef_delta) or (not flag_spread)) and \
								  ((e.alfa > 0.1*Lmaxalfa) or (pn <= q)) )
				else:
					condizione = ((e.alfa > alfa_stop) and \
					              ((deltaf[pn] >= 0.01*coef_delta) or (not flag_spread)) and \
								  ((e.alfa > 0.1*Lmaxalfa) or (pn <= q)) )
				#condizione = True

				#print('deltaf=',deltaf)
				if iprint > 0:
					print('pn=',pn)
					print('alfa=',e.alfa)
					print('deltaf[pn]=',deltaf[pn],' coef_delta=',coef_delta,' flagspread=',flag_spread)
					print('0.1*Lmaxalfa=',0.1*Lmaxalfa)
					print('condizione = ',condizione)

				#input()

				if condizione:
					if iprint >= 0:
						print('EVOLVE POP ?',evolve_pop)
						print('ANALIZZA PUNTO ',pn)
					if iprint > 1:
						print(' nf=',nf,
							  ' size(Lkappa)=',len(self.FILTER.Lkappa),
							  ' size(Ltilde)=',len(self.FILTER.Ltilde),
							  ' size(Lnew)=',len(self.FILTER.Lnew))
						print('d=',[d_coord[i,i] for i in range(n)])
					###################################
					# ASSEGNA PUNTO e ALFA
					###################################
					x = e.x
					f = e.fobs
					xi = e.xi
					if iprint >= 0:
						print('         f= ', f)
					for i in range(n):
						alfa_d[i] = e.alfa_c[i]
					alfa_dense = e.alfa
					alfa_tilde = e.alfa_tilde
					if ndir > len(alfa_tilde):
						alfaint = np.max(alfa_tilde)
						alfa_tilde = np.hstack((alfaint*np.ones(ndir-len(alfa_tilde)),alfa_tilde))
						self.FILTER.Lkappa[pn].alfa_tilde = np.copy(alfa_tilde)

					trovato = True
					if iprint > 1:
						if n < 4:
							print('         x= ', x)
						print('alfa_dense= ', alfa_dense)
						print('alfa_coord= ', alfa_d)
						print('alfa_tilde= ', alfa_tilde)
						print('         f= ', f)

					z = np.copy(x)
					#########################################
					# ANALISI SINGOLO PUNTO
					#########################################
					if iprint > 1:
						print('PROVA', alfa_max)
						print('----------------------------------------------')
						print('nf=%5d   alfamax=%12.5e',nf,alfa_max)
					if iprint >= 2:
						for i in range(n):
							print(' x[',i,']=',x[i])

					''' esplorazione delle direzioni coordinate (continue)'''
					for ii in range(ncont):

						d = np.copy(d_coord[ii,:])
						if iprint > 1:
							print('uso la ',ii,' direzione coord')

						ifront = 0
						for ielle in [1, 2]:
							if ifront == 1:
								break
							if iprint > 1:
								print('ENTRO NELLA LINESEARCH COORD CON alfa=',alfa_d[ii],nf)
							if d[ii] > 0.0:
								if alfa_d[ii] - (bu[ii]-x[ii]) < - 1.e-6:
									alfa = np.maximum(1.e-24,alfa_d[ii])
								else:
									alfa = bu[ii] - x[ii]
									ifront = 1
									if iprint > 1:
										print('a', x[ii],  bl[ii], bu[ii], alfa_d[ii], d[ii], alfa)
										print(' punto espan. sulla front. *')
							else:
								if alfa_d[ii] - (x[ii]-bl[ii]) < -1.e-6:
									alfa = np.maximum(1.e-24,alfa_d[ii])
								else:
									alfa = x[ii] - bl[ii]
									ifront = 1
									if iprint > 1:
										print('b', x[ii],  bl[ii], bu[ii], alfa_d[ii], d[ii], alfa)
										print(' punto espan. sulla front. *')

							if iprint > 1:
								print('alfa_max==================>',alfa_max)
							if np.abs(alfa) <= 1.e-3*np.minimum(1.0,alfa_max):
								d[ii] *= -1
								ifront = 0
								if iprint > 1:
									print(' direzione opposta per alfa piccolo')
									print(' j =',ii,'    d(j) =',d[ii])
									print(' alfa=',alfa,'    alfamax=',alfa_max)
								if (alfa == 0.0) and (ielle == 2):
									self.FILTER.Lkappa[pn].alfa_c[ii] /= 2.0
								alfa = 0.0
								continue
							z[ii] = x[ii] + alfa*d[ii]
							#print('(1) z=',z,' di=',d[ii],' alfa=',alfa)
							#self.CACHE.print()
							found_incache,ind,fob,ciq = self.CACHE.find(z)
							if not found_incache:
								if iprint > 1:
									print('NON trovato punto in cache')
									#print('... z=',z)
								if self.prob.name == "Problem SOLAR 8":
									fob,ciq = self.prob.run_bb(z)
								else:
									fob = self.prob.functs(z)
									ciq = self.prob.fconstriq(z)
								#fob = self.prob.functs(z)
								#ciq = self.prob.fconstriq(z)
								self.CACHE.insert(z,fob,ciq)
								nf += 1
							else:
								if iprint > 1:
									print('trovato punto in cache')
									print(fob,self.prob.functs(z))
							#self.CACHE.print()
							#input()

							fz, violz = self.functs_pen(fob,ciq,eps)
							if iprint > 1:
								print(' fz =',fz,'   alfa =',alfa)
							if iprint >= 2:
								for i in range(n):
									print(' z(',i,')=',z[i])

							fpar = f - gamma*alfa*alfa
							#print('fz=',fz,d[ii])
							#print('fpar=',fpar)
							#if maggiore_stretto(fz,fpar):
							#if self.FILTER.domina(fz,gamma*alfa*alfa):
							if (versione_FORTRAN and maggiore_stretto(fz,fpar)) or \
							   ((not versione_FORTRAN) and self.FILTER.domina(fz,gamma*alfa*alfa)):
								if iprint > 1:
									print('FAILURE STEP')
								if ielle == 2:
									self.FILTER.Lkappa[pn].alfa_c[ii] /= 2.0
								else:
									d[ii] *= -1.0
								z[ii] = x[ii]
							else:
								if iprint > 1:
									print('(Do a projected strong expansion)')
								alfa, z, nf, ifront = self.linesearchbox_cont(x,f,d,alfa,alfa_d,alfa_dense,alfa_tilde,z,fz,ii,alfa_max,iprint,bl,bu,nf,pn,ifront,eps,xi)
								if iprint > 1:
									print('ESCO  DALLA LINESEARCH COORD CON alfa=',alfa_d[ii],nf)

								if np.abs(alfa) >= 1.e-12:
									if iprint > 1:
										print('coord con successo ============')

						#enddo for ielle in [1,2]

					#enddo for ii in range(n)

					if iprint > 0:
						print('-x-x-x-x-x-x-x TRA COORD e DENSE x-x-x-x-x-x-x-x-x-x')
						self.FILTER.print(self.prob)

					''' esplorazione delle direzioni dense, solo rispetto alle
					    continue quindi il ciclo va da 1 a ncont, e solo se
						ncont > 1
					'''
					if ncont > 1:
						for ii in range(ncont):
							#self.FILTER.print(self.prob)
							#input()

							d_dense = np.copy(H[:,ii])
							if iprint > 0:
								print('uso la',ii,' direzione densa',d_dense)

							for ielle in [1,2]:
								if iprint > 1:
									print(' ielle =',ielle)
								if iprint > 1:
									print('ENTRO NELLA LINESEARCH DENSE CON alfa=',alfa_dense,nf)

								z = x+alfa_dense*d_dense
								z = np.maximum(bl,np.minimum(bu,z))
								if np.linalg.norm(z-x) <= 1.e-8:
									d_dense *= -1
									#print('z vicino x')
									continue

								found_incache,ind,fob,ciq = self.CACHE.find(z)
								if not found_incache:
									if self.prob.name == "Problem SOLAR 8":
										fob,ciq = self.prob.run_bb(z)
									else:
										fob = self.prob.functs(z)
										ciq = self.prob.fconstriq(z)
									#fob = self.prob.functs(z)
									#ciq = self.prob.fconstriq(z)
									nf += 1
									self.CACHE.insert(z,fob,ciq)
								fz, viol = self.functs_pen(fob,ciq,eps)

								fpar = f - gamma*alfa_dense*alfa_dense

								if iprint > 1:
									print(' z=',z)
									print(' x=',x)
									print(' d_dense=',d_dense)
									print(' alfa_dense=',alfa_dense)
									print('fz=',fz,' fpar=',fpar,' |z-x|=',np.linalg.norm(z-x))
								#if maggiore_stretto(fz,fpar) or (np.linalg.norm(z-x) <= 1.e-8):
								#if self.FILTER.domina(fz,gamma*alfa_dense*alfa_dense) or (np.linalg.norm(z-x) <= 1.e-8):
								if (versione_FORTRAN and (maggiore_stretto(fz,fpar) or (np.linalg.norm(z-x) <= 1.e-8))) or \
								   ((not versione_FORTRAN) and (self.FILTER.domina(fz,gamma*alfa_dense*alfa_dense) or (np.linalg.norm(z-x) <= 1.e-8))):
									if iprint > 1:
										print('   fz > fpar')
									#############################
									# (Failure step)
									#############################
									if ielle == 2:
										if iprint > 1:
											print('   ielle = 2, riduco alfa')
										self.FILTER.Lkappa[pn].alfa = alfa_dense/2.0
									else:
										if iprint > 1:
											print('   ielle = 1, giro d')
										d_dense *= -1
								else:
									#print('fz NOT > fpar')
									####################################
									# (Do a projected strong expansion)
									####################################
									if iprint > 1:
										print('=====> d_dense = ',d_dense)
									alfa, z, nf, ifront = self.linesearchbox_dense(x,f,d_dense,alfa,alfa_d,alfa_dense,alfa_tilde,z,fz,alfa_max,iprint,bl,bu,nf,pn,ifront,eps,xi)
									if iprint > 1:
										print('=====> d_dense = ',d_dense)
									if iprint > 1:
										print('ESCO  DALLA LINESEARCH DENSE CON alfa=',alfa,nf)

									if np.abs(alfa) >= 1.e-12:
										if iprint > 1:
											print('denso con successo ============')
								continue
							#enddo for ielle in [1,2]
						#enddo for ii in range(n)

						if iprint > 1:
							if hschoice == 1:
								print('index_halton =',index_halton,alfa)
							else:
								print('index_sobol =',index_sobol,alfa)
					'''
					REDUCE THE STEPSIZES: for the dense, coord, primitive directions
					'''
					if self.FILTER.Lkappa[pn].alfa > alfa_stop:
						self.FILTER.Lkappa[pn].alfa /= 2.0
					for i in range(ncont):
						if self.FILTER.Lkappa[pn].alfa_c[i] > alfa_stop:
							self.FILTER.Lkappa[pn].alfa_c[i] /= 2.0

					if iprint > 0:
						print('merge Lnew Ltilde')
					self.FILTER.Ltilde, flag_change = merge(self.FILTER.Lnew,self.FILTER.Ltilde,iprint)

					if iprint > 0:
						self.FILTER.print(self.prob)

					if iprint >= 0:
						print('%5d %5d %5d' % (self.FILTER.sizes()),end='')
						print(' %5d %5d %6s %8d' % (nf, self.CACHE.ncache_hits,'Ltilde',ndir))
						print_line += 1
						if np.mod(print_line,30) == 0:
							print(' Lkap  Ltil  Lnew    nf ncach   type      ndir')
							print('-----------------------------------------------')
					#self.FILTER.init_Lnew()

					istop = self.stop(alfa_d,alfa_dense,alfa_tilde,alfa_max,nf,f,alfa_stop,nf_max)
					if istop >= 2:
						''' stops if nf >= nf_max '''
						break
					############################################
					# FINE ANALISI SINGOLO PUNTO
					############################################

				#endif if condizione
			#endfor pnum,e in enumerate(self.FILTER.Lkappa):

			if iprint >= 0:
				print('------------------------------------------------')
				print('  FINE esplorazione CONTINUE: trovato = ',trovato)
				print('------------------------------------------------')
				#self.FILTER.print_filter(self.prob)
				#input()

			'''
			L'esplorazione delle continue è terminata!
			Generate the new matrix H with a new dense direction
			H ha la seguente struttura:

			      / M  0 \   M matrice ncont x ncont
			  H = |      |   I identità nint x nint
				  \ 0  I /
			'''
			index_halton += 2*n

			if(ncont > 1):
				if hschoice == 1:
					d_dense = np.asarray(Phalton[index_halton-1], dtype = np.double)
					#d_dense = self.halton(index_halton)
				else:
					d_dense, index_sobol = sobol_seq.i4_sobol(ncont, index_sobol)
			else:
				d_dense = np.array([1.0])

			if ncont > 0:
				H = self.gen_base(d_dense)
				H = self.gram_schmidt(H)

				if nint > 0:
					H = np.vstack((H,np.zeros((nint,ncont))))
					H = np.hstack((H,np.vstack((np.zeros((ncont,nint)),np.eye(nint)))))
				if iprint > 0:
					print('================================================')
					print('d_dense = ')
					for j in range(n):
						print('%15.6e ' % d_dense[j],end='')
					print()
					for i in range(n):
						for j in range(n):
							print('%15.6e ' % H[i,j],end='')
						print()
					#print(H)
					print('================================================')
					print('Hit return to continue...')
					input()
			if ncont == 0 and nint > 0:
				H = np.eye(nint)

			flag_successo_discrete_GLOB = False
			if nint > 0:
				'''
				l'esplorazione lungo le direzioni primitive avviene
				solo se ci sono variabili discrete !
				'''
				if iprint > 1:
					print('merge Ltilde Lkappa')
				self.FILTER.Lkappa, flag_change = merge(self.FILTER.Ltilde,self.FILTER.Lkappa,iprint)
				num_merge += 1
				self.FILTER.init_Lnew()
				self.FILTER.init_Ltilde()

				self.FILTER.Lkappa, flag_spread, deltaf = self.spread_ord()
				######################################################
				# Refine self.FILTER.Lkappa
				# exploring just the discrete variables
				######################################################
				for pnum,e in enumerate(self.FILTER.Lkappa):
					if istop >= 2:
						''' stops if nf >= nf_max '''
						break

					pn = pnum
					evolve_pop = False
					viol = 0.0
					if m >= 1:
						x = e.x
						found_incache,ind,fob,ciq = self.CACHE.find(x)
						if not found_incache:
							ciq = self.prob.fconstriq(x)
						viol = np.maximum(viol,np.max(ciq))

					alfa_max = np.maximum(e.alfa,np.max(e.alfa_c))
					if (deltaf[pn] == 0.0) and flag_spread:
						if iprint > 0:
							print(pn)

					condizione = True

					if iprint > 1:
						print('pn=',pn)
						print('alfa=',e.alfa)
						print('deltaf[pn]=',deltaf[pn],' coef_delta=',coef_delta,' flagspread=',flag_spread)
						print('0.1*Lmaxalfa=',0.1*Lmaxalfa)
						print('condizione = ',condizione)

					if condizione:
						if iprint > 1:
							print('EVOLVE POP ?',evolve_pop)
							print('ANALIZZA PUNTO ',pn)
						if iprint > 1:
							print(' nf=',nf,
								  ' size(Lkappa)=',len(self.FILTER.Lkappa),
								  ' size(Ltilde)=',len(self.FILTER.Ltilde),
								  ' size(Lnew)=',len(self.FILTER.Lnew))
							print('d=',[d_coord[i,i] for i in range(n)])
						###################################
						# ASSEGNA PUNTO e ALFA
						###################################
						x = e.x
						f = e.fobs
						xi = e.xi
						for i in range(n):
							alfa_d[i] = e.alfa_c[i]
						alfa_dense = e.alfa
						alfa_tilde = e.alfa_tilde
						if ndir > len(alfa_tilde):
							alfaint = np.max(alfa_tilde)
							alfa_tilde = np.hstack((alfaint*np.ones(ndir-len(alfa_tilde)),alfa_tilde))
							self.FILTER.Lkappa[pn].alfa_tilde = np.copy(alfa_tilde)

						#trovato = True
						if iprint > 1:
							if n < 4:
								print('         x= ', x)
							print('alfa_dense= ', alfa_dense)
							print('alfa_coord= ', alfa_d)
							print('alfa_tilde= ', alfa_tilde)
							print('         f= ', f)

						z = np.copy(x)
						#########################################
						# ANALISI SINGOLO PUNTO
						#########################################
						if iprint > 1:
							print('PROVA', alfa_max)
							print('----------------------------------------------')
							print('nf=%5d   alfamax=%12.5e',nf,alfa_max)
						if iprint >= 2:
							for i in range(n):
								print(' x[',i,']=',x[i])

						''' esplorazione delle nint variabili discrete mediante
						    le ndir direzioni primitive
						'''
						flag_successo_discrete = False
						for idir in range(ndir):
							if (idir > 0) and (self.FILTER.Lkappa[pn].flag_allones):
								break
							''' estrazione della direzione '''
							d = D[:,idir]

							if iprint >= 1:
								print('discrete_search along d',f,d,alfa_tilde[idir])
								input()
							alpha, nf = self.discrete_linesearch(x,d,alfa_d,alfa_dense,alfa_tilde,idir,bl,bu,f,eps,xi,m,ncont,nf,iprint)

							if alpha > 0:
								''' discrete_linesearch SUCCESS '''
								if iprint > 1:
									print('primitive con successo ============')
								flag_successo_discrete = True
								flag_successo_discrete_GLOB = True

								#print('****************   flag_successo_discrete = ',flag_successo_discrete,self.FILTER.sizes())
								#input()
							else:
								d = -d
								if iprint >= 1:
									print('discrete_search along -d',f,d,alfa_tilde[idir])
									input()
								alpha, nf = self.discrete_linesearch(x,d,alfa_d,alfa_dense,alfa_tilde,idir,bl,bu,f,eps,xi,m,ncont,nf,iprint)

								if alpha > 0:
									''' discrete_linesearch SUCCESS '''
									if iprint > 1:
										print('primitive con successo ============')
									flag_successo_discrete = True
									flag_successo_discrete_GLOB = True
									#print('****************   flag_successo_discrete = ',flag_successo_discrete,self.FILTER.sizes())
									#input()
								else:
									#for i in range(ndir):
									if self.FILTER.Lkappa[pn].alfa_tilde[idir] > 1.0:
										self.FILTER.Lkappa[pn].alfa_tilde[idir] = max(1,np.floor(self.FILTER.Lkappa[pn].alfa_tilde[idir]/2))

						found_incache,ind,fob,ciq = self.CACHE.find(x)
						if not found_incache:
							if self.prob.name == "Problem SOLAR 8":
								fob,ciq = self.prob.run_bb(x)
							else:
								fob = self.prob.functs(x)
								ciq = self.prob.fconstriq(x)
							#fob = self.prob.functs(x)
							#ciq = self.prob.fconstriq(x)
							self.CACHE.insert(x,fob,ciq)
						f, viol = self.functs_pen(fob,ciq,eps)

						if iprint > 1:
							print('merge Lnew Ltilde')
						self.FILTER.Ltilde, flag_change = merge(self.FILTER.Lnew,self.FILTER.Ltilde,iprint)
						if iprint >= 0:
							print('%5d %5d %5d' % (self.FILTER.sizes()),end='')
							print(' %5d %5d %6s %8d' % (nf, self.CACHE.ncache_hits,'Ltilde',ndir))
							print_line += 1
							if np.mod(print_line,30) == 0:
								print(' Lkap  Ltil  Lnew    nf ncach   type      ndir')
								print('-----------------------------------------------')
						self.FILTER.init_Lnew()
						#input()

						if not flag_successo_discrete and np.max(self.FILTER.Lkappa[pn].alfa_tilde) == 1:
							'''
							REDUCE THE STEPSIZES: for the dense, coord, primitive directions
							'''
							self.FILTER.Lkappa[pn].xi /= 2
							self.FILTER.Lkappa[pn].flag_allones = True


						istop = self.stop(alfa_d,alfa_dense,alfa_tilde,alfa_max,nf,f,alfa_stop,nf_max)
						if istop >= 2:
							''' stops if nf >= nf_max '''
							break
						############################################
						# FINE ANALISI SINGOLO PUNTO
						############################################

					#endif if condizione
				#endfor pnum,e in enumerate(self.FILTER.Lkappa):

				if iprint >= 0:
					print('------------------------------------------------')
					print('  FINE esplorazione DISCRETE: flag_successo_discrete_GLOB = ',flag_successo_discrete_GLOB)
					print('------------------------------------------------')

			if iprint >= 1:
				print('merge Ltilde Lkappa')
				print('Hit return to continue')
				input()
			self.FILTER.Lkappa, flag_change = merge(self.FILTER.Ltilde,self.FILTER.Lkappa,iprint)
			if not flag_change:
				flag_successo_discrete_GLOB = False
			if (nint > 0) and (iprint >= 0):
				print('------------------------------------------------')
				print('  NUOVO FILTRO: flag_successo_discrete_GLOB = ',flag_successo_discrete_GLOB)
				print('------------------------------------------------')
			num_merge += 1
			self.FILTER.init_Lnew()
			self.FILTER.init_Ltilde()

			if iprint > 0:
				self.FILTER.print_filter(self.prob)
				#input()

			if iprint >= 0:
				print('%5d %5d %5d' % (self.FILTER.sizes()),end='')
				print(' %5d %5d %6s %8d' % (nf, self.CACHE.ncache_hits,'Lkappa',ndir))
				print_line += 1
				if np.mod(print_line,30) == 0:
					print(' Lkap  Ltil  Lnew    nf ncach   type      ndir')
					print('-----------------------------------------------')

			#print('flag_successo_discrete = ',flag_successo_discrete)
			if (nint > 0) and (not flag_successo_discrete_GLOB):

				'''
				Generate a new primitive direction and update matrix D by adding a
				column as first column
				'''
				if iprint >= 2:
					print('(BEFORE) size(D):',np.shape(D),index_primitive)
					print(D)
				iexit = 0
				while iexit == 0:
					# enrich set D
					D, iexit, index_primitive = self.generate_dirs(ncont,nint,D,eta,0,Phalton_int,index_primitive)
					if iexit == 0:
						eta = eta + 0.5
					if eta >= 2.5*(np.linalg.norm(bu - bl)/2):
						#stop execution
						iexit = 1

				ndir = np.shape(D)[1]
				if iprint >= 2:
					print('(AFTER) size(D):',np.shape(D),index_primitive)
					print(D)
					print('Hit return to continue...')
					input()

			#if cambio_eps:
				#manca tutto cambio eps

			if istop >= 2:
				#print('istop = ',istop)
				break

			##############################################
			#              STAMPA PUNTI
			##############################################
			if not trovato:
				coef_delta *= 0.95
				if iprint >= 1:
					print('not trovato',coef_delta)
				if coef_delta < 1.e-16:
					print('coef_delta troppo piccolo')
					break

			#break
			#input()
		#end While True  (ciclo principale di DFMOINT)

		'''
		l'esplorazione è terminata. Si stampano i risultati
		'''
		print('     LISTA  ')
		pn = 0
		j  = 0

		print('     id  FOBs        ',end='')
		for i in range(q-1):
			print('              ',end='')
		print('  VARs     ',end='')
		for i in range(ncont-1):
			print('           ',end='')
		for i in range(nint):
			print('    ',end='')
		print()
		      #      0  2.994938e+03  1.333333e-02  3.000e+00 100 100 100

		for i,e in enumerate(self.FILTER.Lkappa):
			j += 1
			found_incache,ind,fob,ciq = self.CACHE.find(e.x)
			if not found_incache:
				if self.prob.name == "Problem SOLAR 8":
					fob,ciq = self.prob.run_bb(e.x)
				else:
					fob = self.prob.functs(e.x)
					ciq = self.prob.fconstriq(e.x)
				#fob = self.prob.functs(zdelta)
				#ciq = self.prob.fconstriq(zdelta)
				self.CACHE.insert(e.x,fob,ciq)
				
			if m >= 1:
				#ciq = self.prob.fconstriq(e.x)
				viol = np.maximum(0.0,np.max(ciq))
			else:
				viol = 0.0
			if iprint > 1:
				print(i, e.x, viol)
			if viol <= 1.e-3:
				pn += 1
				print(' %6d' % i,end='')
				for fob in e.fobs:
					print(' %13.6e' % fob,end='')
				for i in range(ncont):
					print(' %10.3e' % e.x[i],end='')
				for i in range(nint):
					print(' %3d' % e.x[ncont+i],end='')
				print()

		print('Nondom. = ',j,' di cui amm.',pn)
		print('Numero di Generazioni = ',num_gen,' numero di merge = ',num_merge)
		print('n.f.evals(inc. cache hits) = ',self.CACHE.ncache_hits,' n.f.evals(eff) = ',nf)

		if not print_all_filters:
			num_merge = 0
		fid_Lkappa = open('Lkappa_'+str(num_merge)+'.txt','w')
		for i,e in enumerate(self.FILTER.Lkappa):
			if m >= 1:
				ciq = self.prob.fconstriq(e.x)
				viol = np.maximum(0.0,np.max(ciq))
			else:
				viol = 0.0
			if viol <= 1.e-3:
				print(' %6d' % i,end='',file=fid_Lkappa)
				for fob in e.fobs:
					print(' %13.6e' % fob,end='',file=fid_Lkappa)
				for i in range(ncont):
					print(' %10.3e' % e.x[i],end='',file=fid_Lkappa)
				for i in range(nint):
					print(' %3d' % e.x[ncont+i],end='',file=fid_Lkappa)
				print(file=fid_Lkappa)
		fid_Lkappa.close()

		#print(D)
		return

	def check_pop(self):
		F = np.zeros((0,self.prob.q))
		fid_Lkappa = open('DFMOINT_out.txt','w')
		for i,e in enumerate(self.FILTER.Lkappa):
			if self.prob.m >= 1:
				ciq = self.prob.fconstriq(e.x)
				viol = np.maximum(0.0,np.max(ciq))
			else:
				viol = 0.0
			if viol <= 1.e-3:
				print(' %6d' % i,end='',file=fid_Lkappa)
				F = np.vstack((F,e.fobs))
				for fob in e.fobs:
					print(' %13.6e' % fob,end='',file=fid_Lkappa)
				for i in range(self.prob.ncont):
					print(' %10.3e' % e.x[i],end='',file=fid_Lkappa)
				for i in range(self.prob.nint):
					print(' %3d' % e.x[self.prob.ncont+i],end='',file=fid_Lkappa)
				print(file=fid_Lkappa)
		fid_Lkappa.close()

		return F
