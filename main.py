from datetime import datetime
from algorithm_factory import dfmoint
from problem_factory import problem_factory
from filter_factory import filter_factory
from utility import pareto_efficient
import NOMAD_wrapper as nomad
import numpy as np
import pygmo as pg
import sys
import argparse
import time
import os
import scipy.io as sio

def check_pop(pop,prob):
	#FILTER = filter_factory()
	F = np.zeros((0,prob.q))

	pn = 0
	j  = 0
	m = prob.m
	for x in pop.get_x():
		j += 1
		if m >= 1:
			#ciq = prob.fconstriq(x)
			if prob.name == "Problem SOLAR 8":
				fobs,ciq = prob.run_bb(x)
			else:
				fobs= prob.functs(x)
				ciq = prob.fconstriq(x)
			viol = np.maximum(0.0,np.max(ciq))
		else:
			viol = 0.0
		if viol <= 1.e-3:
			pn += 1
			print(' %6d' % j,end='')
			#fobs = prob.functs(x)
			F = np.vstack((F,fobs))
			for fob in fobs:
				print(' %13.6e' % fob,end='')
			for i in range(prob.ncont):
				print(' %10.3e' % x[i],end='')
			for i in range(prob.nint):
				print(' %3d' % x[prob.ncont+i],end='')
			print()

	neff,_ = F.shape
	if neff > 0:
		efficient_indices = pareto_efficient(F.T)
		#print(efficient_indices)
		F = F[efficient_indices,:]

	j,_ = F.shape
	print('Amm. = ',pn,' di cui nondom.',j)
	return F

parser = argparse.ArgumentParser(prog='main',description='Algorithm DFMOINT')
parser.add_argument('-v','--version',action='store_true',help='Print version number')
parser.add_argument('--constrained',action='store_true',help='Solve constrained problems')
parser.add_argument('--print_filters',action='store_true',help='Print all filters while iterating')
parser.add_argument('--alg',nargs=1,type=str,
                             default=['DFMOINT'],
                             choices=["DFMOINT", "NSGA2", "MACO", "NOMAD"],
                             help='Name of algo to be run')
parser.add_argument('--tol',nargs=1,type=float,
                             default=[1.e-6],
                             help='Tolerance in the stopping condition')
parser.add_argument('--max_fun',nargs=1,type=int,
                             default=[5000],
                             help='Maximum number of function avaluations')
parser.add_argument('--seed',nargs=1,type=int,
                             default=[13578],
                             help='seed for random seq. initialization')
parser.add_argument('--outlev',nargs=1,type=int,
                             default=[0],
                             help='Output level')
parser.add_argument('--export',action='store_true',help='export results')
parser.add_argument('--probid',nargs=1,type=int,
                             default=[1],
                             help='id of problem to be solved')
parser.add_argument('--nvar',nargs=1,type=int,
                             default=[10],
                             choices=[10,15,20,25,30],
                             help='Number of variables')
parser.add_argument('--nint',nargs=1,type=int,
                             default=[-1],
                             help='Number of integer variables (must be < nvar)')
parser.add_argument('--ctype',nargs=1,type=str,
                             default=['a'],
                             choices=['a','b','c','d','e','f','z'],
                             help='type of constraints (z for unconstrained)')
parser.add_argument('--prbsel',nargs=1,type=str,
                             default=['bk1'],
                             #choices=['a','b','c','d','e','f','z'],
                             help='prb of DMS collection')
parser.add_argument('--resdir',nargs=1,type=str,default=[''],help='path were results for Pareto front reference construction can be found')
parser.add_argument('--runall',action='store_true',help='run code on all CEC09 problems')
parser.add_argument('--runall_dms',action='store_true',help='run code on all DMS problems')
args = parser.parse_args(sys.argv[1:])

#print(args.__dict__)
#for key in args.__dict__.keys():
#    print(key.ljust(25), args.__dict__[key])

###################################################################
# PRINT PARAMETERS
###################################################################
print()
print('------------------------------------------------------------')
print('Parameter           : value')
print('------------------------------------------------------------')
print('CONSTRAINED         : %s' % args.constrained)
print('TOL                 : %e' % args.tol[0])
print('MAX_FUN             : %d' % args.max_fun[0])
print('OUTLEV              : %d' % args.outlev[0])
print('PRINT FILETRS       : %s' % args.print_filters)
print('RANDOM SEED         : %d' % args.seed[0])
print('ALGORITHM           : %s' % args.alg[0])
print('EXPORT              : %s' % args.export)
print('RUNALL              : %s' % args.runall)
print('RUNALL_DMS          : %s' % args.runall_dms)
if not args.runall and not args.runall_dms:
	print('PROB ID             : %d' % args.probid[0])
	n       = args.nvar[0]
	nint    = args.nint[0]
	ctype   = args.ctype[0]
	prob_id = args.probid[0]
	prbsel  = args.prbsel[0]
	probl = problem_factory(prob_id,n,ctype,nint,prbsel)
	n,ncont,nint,m,q = probl.prob.getdim()
	print('PROBLEM             : %s (n=%d nint=%d q=%d m=%d)' % (probl.prob.name,n,nint,q,m))
	if args.probid[0] > 1000 and args.probid[0] < 2000:
		print('N.VARS              : %d' % args.nvar[0])
		print('N.INTS              : %d' % args.nint[0])
		print('CONSRT. TYPE        : %s' % args.ctype[0])
print('------------------------------------------------------------')
print()

if args.version:
    print('\nmain.py version 0.1\n')

#print()
#problem_factory().describe()

np.random.seed(args.seed[0])

tol = args.tol[0]
max_fun = args.max_fun[0]
outlev = args.outlev[0]
prob_id = args.probid[0]
prbsel = args.prbsel[0]
alg = args.alg[0]

### SEEDS for NSGA2
### 13578   | 2980    | 109713  | 10739   | 49386   | 18387   | 29648   | 287197  | 67565   | 22656
### NSGA2_1 | NSGA2_2 | NSGA2_3 | NSGA2_4 | NSGA2_5 | NSGA2_6 | NSGA2_7 | NSGA2_8 | NSGA2_9 | NSGA2_10
n       = args.nvar[0]
nint    = args.nint[0]
ctype   = args.ctype[0]
resdir  = args.resdir[0]

if args.constrained:
	ctypes = ['a', 'b', 'c', 'd', 'e', 'f']
else:
	ctypes  = ['z']

if args.runall:
	nvars   = [10,15,20,25,30]
	prids   = [1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010]
if args.runall_dms:
	nvars = [-1]
	prids   = [2000]

TEMPPRB = problem_factory(prob_id=2000)
LISTDMS = TEMPPRB.prob.list.list
LISTUTILI = []
for p in LISTDMS:
	TEMPPRB = problem_factory(prob_id=2000,prbsel=p)
	n1,nc1,ni1,m1,q1 = TEMPPRB.prob.getdim()
	n1,q1 = TEMPPRB.prob.prb.getdim()
	#print(p," ",TEMPPRB.prob.getdim()," -- ",TEMPPRB.prob.prb.getdim(), \
	#			TEMPPRB.prob.functs(TEMPPRB.prob.startp()),TEMPPRB.prob.startp())
	if n1 >= 4:
		LISTUTILI.append(p)
#print(LISTUTILI,len(LISTUTILI),len(LISTDMS))
#input()

if args.runall:
	probs   = [(n,c,id,'bk1') for n in nvars for c in ctypes for id in prids]
elif args.runall_dms:
	probs   = [(-1,c,2000,p) for c in ctypes for p in LISTUTILI]
else:
    if prob_id == 2000:
        probs   = [(n,ctype,prob_id,prbsel)]
    else:
        probs   = [(n,ctype,prob_id,'bk1')]

#print(probs)
#input()

count_problems = 0
for (n,ctype,prob_id,prbsel) in probs:
	count_problems += 1

#print(count_problems)
#input()

if args.export:

	date = datetime.today().strftime("%Y-%m-%d_%H%M%S")
	date = datetime.today().strftime("%Y-%m-%d")
	path = "{}/{}_{}".format(date, alg,args.seed[0])

	print("\nWARNING\n\nI risultati su {} problemi saranno posti nella cartella: {}".format(count_problems,path))
	#print("Hit return to continue")
	#input()
	if not os.path.exists(path):
		if not os.path.exists("{}/".format(date)):
			os.mkdir("{}/".format(date))
		os.mkdir(path)

pop_size = 100
versione_FORTRAN = False

if resdir == '':
	print('\nWARNING\n\nParameter <resdir> is empty. Reference Pareto fronts will not be computed')
else:
	print('\nReference Pareto fronts will be computed using data in: {}'.format(resdir))

	solvers = []
	for subdir in os.listdir(resdir):
		if os.path.isdir(resdir+'\\'+subdir):
			print(subdir)
			solvers.append(subdir)

	print(solvers)

for (n,ctype,prob_id,prbsel) in probs:

	probl = problem_factory(prob_id=prob_id,n=n,ctype=ctype,nint=nint,prbsel=prbsel)

	#print(len(probl.prob.HV))
	#input()
	print(probl.prob.getdim())
	print(probl.prob.get_bounds())
	print(probl.prob.name)

	##########################################################
	## compute of the ref_point for hypervolume computation
	#
	# compute the NADIR point and the IDEAL point
	##########################################################
	prob_front = [None] * len(solvers)
	for (s, solver_name) in enumerate(solvers):

		if (('NOMAD' in solver_name) and (probl.prob.q >= 3)):
			prob_front[s] = np.empty((probl.prob.q, 0))
			continue
		if (('NOMADcon' in solver_name) and (probl.prob.q >= 3)):
			prob_front[s] = np.empty((probl.prob.q, 0))
			continue

		# prob_front_ = sio.loadmat("2022-02-03_20000/{}/{}.mat".format(solver_name, probl.prob.name))
		try:
			prob_front_ = sio.loadmat("{}/{}/{}.mat".format(resdir, solver_name, probl.prob.name))
			npoints, nobj = prob_front_['F'].shape
			prob_front[s] = np.transpose(prob_front_['F'])
			print(prob_front[s].shape)
		except:
			prob_front[s] = np.empty((probl.prob.q, 0))
			continue

		if npoints == 0:
			prob_front[s] = np.reshape(prob_front[s], (probl.prob.q, 0))

	F_union = np.concatenate(prob_front, axis=1)

	efficient_indices = pareto_efficient(F_union)
	F_true = F_union[:, efficient_indices]

	nobj, npoints = F_union.shape
	if npoints > 0:
		ref_point = np.max(F_true, axis=1)
		idl_point = np.min(F_true, axis=1)
	else:
		ref_point = np.inf*np.ones(probl.prob.q)
		idl_point = np.inf * np.ones(probl.prob.q)
	probl.prob.ref_point = ref_point
	probl.prob.idl_point = idl_point
	#print(probl.prob)
	#print(F_union.shape)
	#print('ref_point:',probl.prob.ref_point)
	#input()
	#############################################################
	#############################################################

	#probl.prob.functs(probl.prob.startp())
	#probl.prob.fconstriq(probl.prob.startp())
	#input()

	if args.export:
		print('{}/{}.mat'.format(path, probl.prob.name))
		#if os.path.exists('{}/{}.mat'.format(path, probl.prob.name)):
		#	continue

	#if (prob_id == 1007) and (n == 25) and (ctype=='d'):
	#	F = np.zeros((0,probl.prob.q))
	#	sio.savemat('{}/{}.mat'.format(path, probl.prob.name), {'F': F})
	#	print('Salvato file mat con F=[]')
	#	input()
	#	continue

	if alg == 'NSGA2':
	    algo = pg.nsga2(gen=int(max_fun / pop_size)-1,seed=args.seed[0])
	elif alg == 'MACO':
	    algo = pg.maco(gen=int(max_fun / pop_size)+1,seed=args.seed[0])

	if alg == 'DFMOINT':
		algo = dfmoint(problem=probl.prob,n=n,ctype=ctype)
		start = time.time()
		algo.run(alfa_stop=tol,iprint=outlev,nf_max=max_fun,print_all_filters=args.print_filters,versione_FORTRAN=versione_FORTRAN)
		stop = time.time()
		F = algo.check_pop()
		if args.export:
			#sio.savemat('{}/{}.mat'.format(path, probl.prob.name), {'F': F})
			sio.savemat('{}/{}.mat'.format(path, probl.prob.name), {'F': F, 'HV': probl.prob.HV})
		print('Total running time = ',stop-start)
		print(algo.CACHE.store_pd)
	elif alg == 'NOMAD':
		if probl.prob.q < 3:
			start = time.time()
			nomad.alg_NOMAD(probl.prob, max_fun)
			stop = time.time()
			F = nomad.check_pop(probl.prob)
			if args.export:
				#sio.savemat('{}/{}.mat'.format(path, probl.prob.name), {'F': F})
				sio.savemat('{}/{}.mat'.format(path, probl.prob.name), {'F': F, 'HV': probl.prob.HV})
			print('Total running time = ',stop-start)
	else:
	    alg = pg.algorithm(algo)
	    alg.set_verbosity(1)
	    #print(alg)
	    pgprob = pg.problem(probl.prob)
	    #print(len(probl.prob.HV))
	    #input()
	    pop = pg.population(pgprob,size=pop_size,seed=args.seed[0])
	    #print(len(probl.prob.HV))
	    #input#()
	    start = time.time()
	    pop = alg.evolve(pop)
	    #for x in pop.get_x():
	    #    print(x)
	    #input()
	    stop = time.time()
	    #print(len(probl.prob.HV))
	    #input()
	    F = check_pop(pop,probl.prob)
	    #print(probl.prob.FRONT.F)
	    #print(len(probl.prob.HV))
	    #input()
	    print('Total running time = ',stop-start)
	    if args.export:
	    	sio.savemat('{}/{}.mat'.format(path, probl.prob.name), {'F': F, 'HV': probl.prob.HV})

		#print("{}/{}.txt".format(path,probl.prob.name))
	    #sio.savemat('{}/{}_{}.mat'.format(path, problem, n), {'F': f_list, 'func_eval': total_func_eval})
