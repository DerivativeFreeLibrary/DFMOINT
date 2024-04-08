# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:36:34 2020

@author: tomma
"""

#import PyNomad
from filter_factory import filter_factory
import numpy as np
import pandas as pd
import math
import os

def check_pop(prob):
	F = np.zeros((0,prob.q))
	if os.path.exists('nomad_out.0.txt'):
		data = np.array(pd.read_csv('nomad_out.0.txt',sep=",",header=None))
		m = prob.m
		n = prob.n
		q = prob.q
		F = data[:,:q]
		X = data[:,q:q+n]

		#FILTER = filter_factory()

		pn = 0
		j  = 0
		for x in X:
			j += 1
			if m >= 1:
				ciq = prob.fconstriq(x)
				viol = np.maximum(0.0,np.max(ciq))
			else:
				viol = 0.0
			if viol <= 1.e-3:
				pn += 1
				print(' %6d' % j,end='')
				for fob in prob.functs(x):
					print(' %13.6e' % fob,end='')
				for i in range(prob.ncont):
					print(' %10.3e' % x[i],end='')
				for i in range(prob.nint):
					print(' %3d' % x[prob.ncont+i],end='')
				print()

		print('Nondom. = ',j,' di cui amm.',pn)

	return F

def alg_NOMAD(prob,max_fun):


    n = prob.n
    m = prob.m
    q = prob.q
    ncont = prob.ncont
    nint = prob.nint
    funct = prob.functs
    constr = prob.fconstriq
    (lb,ub) = prob.get_bounds()
    x0 = np.copy(prob.startp())

    def fitfun(x):
        #n    = x.get_n()
        xx   = np.array([x.get_coord(i) for i in range(n)])
        floc = funct(xx)
        gloc = constr(xx)

        #print('fitfun: ',f)
        #
        for i in range(q):
            x.set_bb_output(i,floc[i])

        if m > 0:
            for i in range(m):
                x.set_bb_output(i+q,gloc[i])

        return 1

    vtype = np.array(['R' for i in range(n)])
    for i in range(nint):
        vtype[ncont+i] = 'I'

    bb_input_type = 'BB_INPUT_TYPE ('
    for i in range(n):
        bb_input_type += ' '+vtype[i]
    bb_input_type += ')'

    print(bb_input_type)

    bb_output_type = 'BB_OUTPUT_TYPE'
    for i in range(q):
        bb_output_type += ' OBJ'

    for i in range(m):
        bb_output_type += ' PEB'

    params = ['ASYNCHRONOUS FALSE',bb_output_type,                 'MAX_BB_EVAL '+str(max_fun),bb_input_type,'DISPLAY_STATS BBE OBJ CONS_H','STATS_FILE nomad_out.txt']
    params = ['ASYNCHRONOUS FALSE',bb_output_type,'DISABLE MODELS','MAX_BB_EVAL '+str(max_fun),bb_input_type,'DISPLAY_STATS BBE OBJ CONS_H','STATS_FILE nomad_out.txt']
    params = [bb_output_type,'DISABLE MODELS','MAX_BB_EVAL '+str(max_fun),bb_input_type,'DISPLAY_STATS BBE OBJ CONS_H','STATS_FILE nomad_out.txt']
    params = [bb_output_type,                 'MAX_BB_EVAL '+str(max_fun),bb_input_type,'DISPLAY_STATS BBE OBJ CONS_H','STATS_FILE nomad_out.txt']
    params = [bb_output_type,                 'MAX_BB_EVAL '+str(max_fun),bb_input_type,'DISPLAY_STATS BBE BBO CONS_H','STATS_FILE nomad_out.txt']
    params = [bb_output_type,
             'MULTI_OVERALL_BB_EVAL '+str(max_fun),
             'MAX_BB_EVAL '+str(math.ceil(max_fun/10)),
             bb_input_type,
			 'NM_SEARCH yes','DISABLE MODELS',
             'DISPLAY_STATS OBJ CONS_H - SOL',
             'STATS_FILE nomad_out.txt OBJ, , SOL, ,CONS_H']

    print(params)

    [ x_return , f_return , h_return, nb_evals , nb_iters ,  stopflag ] = PyNomad.optimize(fitfun,x0,lb,ub,params)
    print("ciaooooooo")
    print(x_return)
    print(f_return)
    return
