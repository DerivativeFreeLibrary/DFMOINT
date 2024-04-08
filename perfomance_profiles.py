import glob
import pygmo as pg
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import scipy.io as sio
import sys
import pickle

import time
import os
from matplotlib.backends.backend_pdf import PdfPages

import argparse


# performance must be a matrix n * m where n is the number of problems and m is the number of solvers
from purity import compute_purity
from spread import compute_gamma_delta
from hypervol import compute_hypervol
from utility import pareto_efficient
from problem_factory import problem_factory, transform_point, transform_front
from filter_factory import front_factory

nomi_x_legend = {'NSGA2_worst': 'NSGA-II (worst)', \
                 'NSGA2_best': 'NSGA-II (best)', \
                 'NSGA2_median': 'NSGA-II (median)', \
                 'NSGA2_worsthv': 'NSGA-II (worst)', \
                 'NSGA2_besthv': 'NSGA-II (best)', \
                 'NSGA2_medianhv': 'NSGA-II (median)', \
                 'NOMAD': 'BiMADS (3.9.1)', \
                 'NOMAD_13578': 'BiMADS (3.9.1)', \
                 'DFMOINTold': 'DFMOINT (old results)', \
                 'DFMOINT': 'DFMOINT', \
                 'DFMOINT_13578': 'DFMOINT', \
                 'DFMOINTcon': 'DFMOINT', \
                 'DFMOINTcon_2022-03-28': 'DFMOINT', \
                 'DFMOINT_2022-03-30': 'DFMOINT', \
                 'NSGA2con_worst': 'NSGA-II (worst)', \
                 'NSGA2con_best': 'NSGA-II (best)', \
                 'NSGA2con_median': 'NSGA-II (median)', \
                 'NSGA2con_1': 'NSGA-II (run 1)', \
                 'NSGA2con_2': 'NSGA-II (run 2)', \
                 'NSGA2con_3': 'NSGA-II (run 3)', \
                 'NSGA2con_4': 'NSGA-II (run 4)', \
                 'NSGA2con_5': 'NSGA-II (run 5)', \
                 'NSGA2con_6': 'NSGA-II (run 6)', \
                 'NSGA2con_7': 'NSGA-II (run 7)', \
                 'NSGA2con_8': 'NSGA-II (run 8)', \
                 'NSGA2con_9': 'NSGA-II (run 9)', \
                 'NSGA2con_10': 'NSGA-II (run 10)', \
                 'NOMAD_noNM': 'BiMADS (3.9.1 w/o N.M.)', \
                 'NOMADnomodel': 'BiMADS (3.9.1 w/o mod.)', \
                 'NOMADcon': 'BiMADS (3.9.1)'}

colori_x_plot = {'NSGA2_worst': 'green', \
                 'NSGA2_best': 'green', \
                 'NSGA2_median': 'green', \
                 'NSGA2_worsthv': 'green', \
                 'NSGA2_besthv': 'green', \
                 'NSGA2_medianhv': 'green', \
                 'NSGA2con_worst': 'green', \
                 'NSGA2con_best': 'green', \
                 'NSGA2con_median': 'green', \
                 'NOMAD': 'orange', \
                 'NOMAD_13578': 'orange', \
                 'NOMAD_noNM': 'orange', \
                 'NOMADnomodel': 'orange', \
                 'NOMADcon': 'orange', \
                 'DFMOINTold': 'blue', \
                 'DFMOINT': 'blue', \
                 'DFMOINT_13578': 'blue', \
                 'DFMOINTcon': 'blue', \
                 'DFMOINTcon_2022-03-28': 'blue', \
                 'DFMOINT_2022-03-30': 'blue'}

tratti_x_plot = {'NSGA2_worst': '-.', \
                 'NSGA2_best': '--', \
                 'NSGA2_median': '-', \
                 'NSGA2_worsthv': '-.', \
                 'NSGA2_besthv': '--', \
                 'NSGA2_medianhv': '-', \
                 'NSGA2con_worst': '-.', \
                 'NSGA2con_best': '--', \
                 'NSGA2con_median': '-', \
                 'NOMAD': '-', \
                 'NOMAD_13578': '-', \
                 'NOMAD_noNM': '-', \
                 'NOMADnomodel': '-', \
                 'NOMADcon': '-', \
                 'DFMOINTold': '-', \
                 'DFMOINT': '-', \
                 'DFMOINT_13578': '-', \
                 'DFMOINTcon': '-', \
                 'DFMOINTcon_2022-03-28': '-', \
                 'DFMOINT_2022-03-30': '-'}

def compute_profiles(performance_in):
    performance = np.copy(performance_in)
    #performance deve essere una matrice di dimensione: n_problems x n_solvers
    inf_factor = 2
    worst_value = 0
    performance = performance[~np.isinf(performance)[:, 0],:]
    n_problems, n_solvers = performance.shape
    best_values = np.nanmin(performance, axis=1)

    for (i, bv) in enumerate(best_values):
        if np.isnan(bv):
            bvv = np.inf
        else:
            bvv = bv
        performance[i, :] = np.divide(performance[i, :], bvv, dtype=float)
        current_worst = max(performance[i, :])

        if current_worst > worst_value and current_worst != np.inf:
            worst_value = max(performance[i, :])

    performance[np.isnan(performance)] = inf_factor * worst_value
    performance[np.isinf(performance)] = inf_factor * worst_value

    if True:
        # Compute cumulative
        taus = np.arange(1, min(11, np.ceil(worst_value)), 0.001)
        cumulatives = np.zeros(shape=(len(taus), n_solvers))
        for t, tau in enumerate(taus):
            for s in range(n_solvers):
                cumulatives[t, s] = len(np.where(performance[:, s] <= tau)[0]) / n_problems
    else:
        performance = np.sort(performance,axis=0)

    return cumulatives, taus

def plot_data_profiles(H,N,gate,names=None,pdfpagefile=None):
    # This subroutine produces a data profile as described in:
    #
    # Benchmarking Derivative - Free Optimization Algorithms
    # Jorge J.More' and Stefan M. Wild
    # SIAM J.Optimization, Vol. 20(1), pp.172 - 191, 2009.
    #
    # H contains a three dimensional array of function values.
    # H(f, p, s) = function value  # f for problem p and solver s.
    # N is an np - by - 1 vector of(positive) budget units.If simplex
    # gradients are desired, then N(p) would be n(p) + 1, where n(p) is
    # the number of variables for problem p.
    # gate is a positive constant reflecting the convergence tolerance.

    nf, npp, ns = H.shape # Grab the dimensions
    #prob_min = np.min(np.min(H, axis=0), axis=1) # The global minimum for each problem
    #Q = ma.masked_array(H, np.isnan(H))
    #prob_min = Q.min(axis=0).min(axis=1)
    #prob_max = H[0,:,0] # The starting value for each problem
    #prob_max = Q.max(axis=0).max(axis=1)
    #alg_max = np.zeros(ns)
    #prob_max = np.zeros(npp)
    #for p in range(npp):
    #    for s in range(ns):
    #        if (np.sum(H[np.isfinite(H[:, p, s]), p, s]) == 0):
    #            alg_max[s] = -np.inf
    #        else:
    #            alg_max[s] = np.max(H[np.isfinite(H[:, p, s]), p, s])

    #    if (np.max(alg_max) == -np.inf):
    #        prob_max[p] = np.inf
    #    else:
    #        prob_max[p] = np.max(alg_max)

    # For each problem and solver, determine the number of
    # N - function bundles(e.g. - gradients) required to reach the cutoff value
    print('0:',N)
    T = np.zeros((npp, ns))
    for p in range(npp):
        print('problem: {}'.format(p),end='')
        #cutoff = prob_min[p] + gate * (prob_max[p] - prob_min[p])
        ##cutoff = (1/(1-gate))*(1/prob_min[p])
        #cutoff = (1-gate)*(1/prob_min[p])
        #print('1:',prob_min[p],' ',prob_max[p],' ',cutoff)
        for s in range(ns):
            #nfevs = np.where(H[:,p,s] <= cutoff)
            #nfevs = np.where(1/H[:,p,s] >= cutoff)
            nfevs = np.where(H[:, p, s] >= (1-gate))
            #print('2:',H[0,p,s])
            #print('3:',nfevs[0])
            if len(nfevs[0]) == 0:
                T[p, s] = np.nan
            else:
                T[p, s] = nfevs[0][0] / N[p] + 20
            print(' {}'.format(T[p,s]),end='')
        print()
    print('Hit return to continue')
    input()

    # Replace all NaN's with twice the max_ratio and sort.
    max_data = np.max(T[~np.isnan(T)])
    T[np.isnan(T)] = 2 * max_data
    T = np.sort(T,axis=0)

    #print()
    #print(nf,npp,ns)
    #print(max_data)
    # For each solver, plot stair graphs with markers.

    plt.figure(figsize=(4.68, 4.68), dpi=300, layout='constrained')
    plt.title(r"$\tau = {}$".format(gate))
    plt.ylabel("Cumulative")
    plt.xlabel(r'groups of $n+1$ BB calls')

    #cumulative, tau = compute_profiles(T)
    #plot_performance_profiles(cumulative, tau, pdfpagefile=pdfpagefile)
    #return

    edges = [(i+1)/npp for i in range(npp)]
    edges = edges[1:]
    max_ratio = 0.0
    for s in range(ns):
        if 'NSGA2' in names[s]:
            continue
        else:
            if names is not None:
                try:
                    lgd = nomi_x_legend[names[s]]
                except:
                    lgd = names[s]
                try:
                    clr = colori_x_plot[names[s]]
                except:
                    clr = (1 / (s + 1), 1 - 1 / (s + 1), 0.5)
                try:
                    lst = tratti_x_plot[names[s]]
                except:
                    lst = '-'
            else:
                lgd = "Solver_" + str(s + 1)
                clr = (1/(s+1),1-1/(s+1),0.5)
                lst = '-'
            # aggiunta la riga di sotto per fare si che il plot non torni a 0
            T[-1,s] = 2*max_data
            ind = np.where(T[:,s] >= 1.1 * max_data)
            #print(ind[0])
            if (len(ind[0])>0):
                #print(ind[0][0])
                max_ratio = np.maximum(max_ratio,edges[ind[0][0]])
            plt.stairs(edges,T[:,s],linestyle=lst,label=lgd, color=clr)
            #print(edges[0],edges[-1])
            #print(T[0,s],T[-3:,s])
            #for p in range(npp-1):
            #    print(edges[p],T[p,s])
            #print(T[-1,s])

    #input()
    #edges.insert(0,0)
    linestyles = ['-', '--']
    is_nsga2 = False
    for s in range(ns):
        if 'NSGA2' in names[s]:
            is_nsga2 = True
            break

    if is_nsga2:
        count = 0
        edges = [(i + 1) / npp for i in range(npp)]
        edges = edges[1:]
        #compute profile for the runs of NSGA2
        y2 = -np.inf*np.ones(npp)
        y1 =  np.inf*np.ones(npp)
        ym =  np.zeros(npp)
        for s in range(ns):
            if 'NSGA2' in names[s]:
                if names is not None:
                    try:
                        lgd = nomi_x_legend[names[s]]
                    except:
                        lgd = names[s]
                else:
                    lgd = "Solver_" + str(s + 1)
                #plt.stairs(edges,T[:,s],label=lgd)
                #plt.plot(T[:,s],edges,'-')
                T[-1,s] = 2 * max_data
                y1 = np.minimum(y1,T[:,s])
                y2 = np.maximum(y2,T[:,s])
                ym += T[:,s]
                count += 1
        #ym = (y1+y2)/2.0
        #ym = np.median(T,axis=1)
        ym /= count
        #plt.stairs(edges, y1)
        #plt.stairs(edges, y2)
        ym[-1] = 2 * max_data
        plt.stairs(edges, ym, color='green', label='NSGA2 (average)')
        ind = np.where(ym >= 1.1 * max_data)
        #print(ind[0])
        if (len(ind[0]) > 0):
            #print(ind[0][0])
            max_ratio = np.maximum(max_ratio, edges[ind[0][0]])

        edges = [(i+1)/npp for i in range(npp)]
        #plt.plot(ym, edges, label='NSGA2 mean')
        plt.fill_betweenx(edges, y1, y2, facecolor='gray',alpha=0.3,step='post')

    # Axis properties are set so that failures are not shown, but with the
    # max_ratio data points shown. This highlights the "flatline" effect.
    plt.xlim((-100, 1.1*max_data))
    plt.ylim((0, np.minimum(1.0,1.2*max_ratio)))
    #if (gate == 0.1):
    #    plt.ylim((0,0.4))
    #elif (gate == 5.e-2):
    #    plt.ylim((0,0.3))
    #elif (gate == 1.e-2):
    #    plt.ylim((0,0.2))
    #else:
    #    plt.ylim((0, 1))
    if (gate == 0.5):
        plt.legend(loc='lower right',ncol=1)
    if pdfpagefile is not None:
        #plt.savefig(pdfpagefile, format="pdf", papertype='a4')
        plt.savefig(pdfpagefile, format="pdf")
    else:
        plt.show()
    plt.close()

def plot_performance_profiles_wild(H,title,names=None,pdfpagefile=None):
    npp, ns = H.shape # Grab the dimensions
    print('{}: npp = {} ns = {}\n'.format(title,npp,ns))

    # Compute ratios and divide by smallest element in each row.
    T = np.copy(H)
    T = T[~np.isinf(T)[:, 0], :]
    npp, ns = T.shape
    print('{}: npp = {} ns = {}\n'.format('  T',npp,ns))

    #T[np.where(T <= 1.e-3)] = 1.e-3
    best = np.nanmin(T,axis=1)
    worst= np.nanmax(T,axis=1)
    print('best worst')
    print(np.hstack((best.reshape(-1,1),worst.reshape(-1,1))))
    input()
    bestM = np.tile(np.reshape(best,(npp,1)),(1,ns))
    r = T/bestM
    print(bestM)
    print()
    print(r)
    input()

    max_ratio = np.max(r[~np.isnan(r)])
    r[np.isnan(r)] = 2 * max_ratio
    r = np.sort(r,axis=0)

    plt.figure(figsize=(4.68, 4.68), dpi=300)
    plt.title("{}".format(title))

    edges = [(i+1)/npp for i in range(npp)]
    edges = edges[1:]
    for s in range(ns):
        if False:
            continue
        else:
            if names is not None:
                try:
                    lgd = nomi_x_legend[names[s]]
                except:
                    lgd = names[s]
                try:
                    clr = colori_x_plot[names[s]]
                except:
                    clr = (1 / (s + 1), 1 - 1 / (s + 1), 0.5)
                try:
                    lst = tratti_x_plot[names[s]]
                except:
                    lst = '-'
            else:
                lgd = "Solver_" + str(s + 1)
                clr = (1/(s+1),1-1/(s+1),0.5)
                lst = '-'
            # aggiunta la riga di sotto per fare si che il plot non torni a 0
            r[-1, s] = 2 * max_ratio
            plt.stairs(edges,r[:,s],linestyle=lst,label=lgd,color=clr)
            #plt.plot(r[:,s],edges,label=lgd)

    # Axis properties are set so that failures are not shown, but with the
    # max_ratio data points shown. This highlights the "flatline" effect.
    plt.ylabel("Cumulative")
    plt.xlabel(r'$\tau$')
    if max_ratio < 11:
        plt.xlim((1, 1.1*max_ratio))
    else:
        plt.xlim((1, 11))
    plt.ylim((0, 1))
    plt.legend(loc='lower right',ncol=1)
    if pdfpagefile is not None:
        #plt.savefig(pdfpagefile, format="pdf", papertype='a4')
        plt.savefig(pdfpagefile, format="pdf")
    else:
        plt.show()
    plt.close()

def plot_performance_profiles(cumulative, tau, names=None, metric_name=None, pdfpagefile=None):

    n_solvers = cumulative.shape[1]
    #cumulative = np.insert(cumulative, 0, 0, axis=0)
    #tau = np.insert(tau, 0, 1)
    plt.figure(figsize=(4.68, 4.68), dpi=300)
    plt.title(metric_name if metric_name is not None else "Performance")
    linestyles = ['-', '--']

    for s in range(n_solvers):
        if names is not None:
            try:
                lgd = nomi_x_legend[names[s]]
            except:
                lgd = names[s]
            try:
                clr = colori_x_plot[names[s]]
            except:
                clr = (1/(s+1),1-1/(s+1),0.5)
            try:
                lst = tratti_x_plot[names[s]]
            except:
                lst = '-'
        else:
            lgd = "Solver_"+str(s+1)
            clr = (1/(s+1),1-1/(s+1),0.5)
            lst = '-'
        #plt.plot(tau, cumulative[:, s], '-'*(s+1),  label=names[s] if names is not None else "Solver_"+str(s+1))
        #plt.plot(tau, cumulative[:, s], '-',  label=names[s] if names is not None else "Solver_"+str(s+1))
        plt.plot(tau, cumulative[:, s], lst,  label=lgd, color=clr)

    if "Purity" in metric_name:
        plt.ylabel("Cumulative")

    plt.xlabel(r'$\tau$')
    plt.xlim((1, tau[-1]))
    plt.ylim((0, 1))

    if "Purity" in metric_name:
        plt.ylabel("Cumulative")
        plt.legend(loc='lower right',ncol=1)
        #plt.legend(loc='upper right')

    if pdfpagefile is not None:
        #plt.savefig(pdfpagefile, format="pdf", papertype='a4')
        plt.savefig(pdfpagefile, format="pdf")
    else:
        plt.show()
    plt.close()

def plot_performance_profiles_new(cumulative, tau, names=None, metric_name=None, pdfpagefile=None):

    n_solvers = cumulative.shape[1]
    npp = cumulative.shape[0]
    print('{}: npp = {} ns = {}\n'.format(metric_name,npp,n_solvers))

    #cumulative = np.insert(cumulative, 0, 0, axis=0)
    #tau = np.insert(tau, 0, 1)
    plt.figure(figsize=(4.68, 4.68), dpi=300)
    plt.title(metric_name if metric_name is not None else "Performance")
    linestyles = ['-', '--']

    for s in range(n_solvers):
        if 'NSGA2' in names[s]:
            continue
        if names is not None:
            try:
                lgd = nomi_x_legend[names[s]]
            except:
                lgd = names[s]
            try:
                clr = colori_x_plot[names[s]]
            except:
                clr = (1/(s+1),1-1/(s+1),0.5)
            try:
                lst = tratti_x_plot[names[s]]
            except:
                lst = '-'
        else:
            lgd = "Solver_"+str(s+1)
            clr = (1/(s+1),1-1/(s+1),0.5)
            lst = '-'
        #plt.plot(tau, cumulative[:, s], '-'*(s+1),  label=names[s] if names is not None else "Solver_"+str(s+1))
        #plt.plot(tau, cumulative[:, s], '-',  label=names[s] if names is not None else "Solver_"+str(s+1))
        plt.plot(tau, cumulative[:, s], lst,  label=lgd, color=clr)

    is_nsga2 = False
    for s in range(n_solvers):
        if 'NSGA2' in names[s]:
            is_nsga2 = True
            break

    if is_nsga2:
        count = 0
        y2 = -np.inf*np.ones(npp)
        y1 =  np.inf*np.ones(npp)
        ym =  np.zeros(npp)
        for s in range(n_solvers):
            if 'NSGA2' in names[s]:
                #plt.stairs(edges,T[:,s],label=lgd)
                #plt.plot(T[:,s],edges,'-')
                y1 = np.minimum(y1,cumulative[:,s])
                y2 = np.maximum(y2,cumulative[:,s])
                ym += cumulative[:,s]
                count += 1
        ym /= count
        plt.plot(tau, ym, '-',  label='NSGA2 (average)', color='green')
        plt.plot(tau, y1, '-',  label='NSGA2 (lower)', color='black')
        plt.plot(tau, y2, '-',  label='NSGA2 (upper)', color='red')
        edges = [(i+1)/npp for i in range(npp)]
        print(len(tau))
        print(len(edges))
        #plt.fill_betweenx(edges, y1, y2, facecolor='gray',alpha=0.3,step='post')
        plt.fill_betweenx(edges, y1, y2, facecolor='gray',step='post')

    #if "Purity" in metric_name:
    plt.ylabel("Cumulative")

    plt.xlabel(r'$\tau$')
    plt.xlim((1, tau[-1]))
    plt.ylim((0, 1))

    if "Purity" in metric_name:
        plt.ylabel("Cumulative")
        plt.legend(loc='lower right',ncol=1)
        #plt.legend(loc='upper right')

    if pdfpagefile is not None:
        #plt.savefig(pdfpagefile, format="pdf", papertype='a4')
        plt.savefig(pdfpagefile, format="pdf")
    else:
        plt.show()
    plt.close()

def intersect(a,b):
    av = a.view([('', a.dtype)] * a.shape[1]).ravel()
    bv = b.view([('', b.dtype)] * b.shape[1]).ravel()
    return np.intersect1d(av, bv).view(a.dtype).reshape(-1, a.shape[1])

def compute_best(probs,resdir,maxfun,solvername=None,solvers=None):
    '''
    estrae da N prove di NSGA2 il "best" e il "worst" e crea dei risultati
    come se questi due fossero nuovi solutori
    '''
    
    randruns = range(1,11)
    #randsolvers = ['NSGA2con_{}'.format(i) for i in randruns]
    randsolvers = ['{}_{}'.format(solvername,i) for i in randruns]
    randsolvers = solvers

    for (n,ctype,prob_id,prbsel) in probs:
        probl = problem_factory(prob_id=prob_id,n=n,ctype=ctype,nint=nint,prbsel=prbsel)
        print(probl.prob.name,end='')
        prob_front = [None] * len(randsolvers)
        HV = np.zeros((maxfun,len(randsolvers)))
        for (s, solver_name) in enumerate(randsolvers):
            #print(s,' ',solver_name)
            #print(' ',s,end='')
            #prob_front_ = sio.loadmat("2022-03-15/{}/{}.mat".format(solver_name, probl.prob.name))
            prob_front_ = sio.loadmat("{}/{}/{}.mat".format(resdir,solver_name, probl.prob.name))

            npoints,nobj = prob_front_['F'].shape
            prob_front[s] = np.transpose(prob_front_['F'])
            HV[:,s] = (np.transpose(prob_front_['HV'])[1:maxfun+1,0])
            HV[np.where(HV <= 0),s] = np.nan

            if npoints == 0:
                prob_front[s] = np.reshape(prob_front[s],(probl.prob.q,0))
            #print(prob_front[s].shape,solver_name,s)
            if npoints > 0:
                efficient_indices = pareto_efficient(prob_front[s])
                prob_front[s] = prob_front[s][:, efficient_indices]
            _,npt = prob_front[s].shape
            print(' {:1d} {:4d} |'.format(s,npt),end='')

        #print()
        #print(probl.prob.name,end='')

        F_union = np.concatenate(prob_front, axis=1)
        _,n_points = F_union.shape
        if n_points == 0:
            F_true = np.copy(F_union)
        else:
            efficient_indices = pareto_efficient(F_union)
            F_true = F_union[:, efficient_indices]
        _,n_points = F_true.shape
        #print(' {}'.format(n_points),end='')
        if n_points > 0:
            ratios = []
            for (s, solver_name) in enumerate(randsolvers):
                #print(s,' ',solver_name)
                F_prime = intersect(F_true.T,prob_front[s].T)
                F_prime = F_prime.T
                _,dimps = F_prime.shape
                #print(' {:1d} {:4d} {:10f}|'.format(s,dimps,dimps/n_points),end='')
                ratios.append(dimps/n_points)
            ind = np.argsort(ratios)
            indmin = ind[0]
            indmax = ind[-1]
            media  = np.sum(ratios)/len(randsolvers)
            indmed = -1
            dist = np.inf
            for i,r in enumerate(ratios):
                 if np.abs(media-r) < dist:
                     dist = np.abs(media-r)
                     indmed = i
        else:
            indmax = 0
            indmin = 0
            indmed = 0
        #print()
        #print(ratios)
        #print(ind)
        #print(indmin,' ',indmax)
        print(prob_front[indmin].shape,end='')
        print(' ',prob_front[indmed].shape,end='')
        print(' ',prob_front[indmax].shape,end='')

        pathw = '{}/{}'.format(resdir,solvername+'_worst')
        pathm = '{}/{}'.format(resdir,solvername+'_median')
        pathb = '{}/{}'.format(resdir,solvername+'_best')
        if not os.path.exists(pathw):
            os.mkdir(pathw)
        if not os.path.exists(pathm):
            os.mkdir(pathm)
        if not os.path.exists(pathb):
            os.mkdir(pathb)

        sio.savemat('{}/{}.mat'.format(pathw, probl.prob.name), {'F': prob_front[indmin].T,'HV': HV[:,indmin]})
        sio.savemat('{}/{}.mat'.format(pathm, probl.prob.name), {'F': prob_front[indmed].T,'HV': HV[:,indmed]})
        sio.savemat('{}/{}.mat'.format(pathb, probl.prob.name), {'F': prob_front[indmax].T,'HV': HV[:,indmax]})
        #input()
        print()

    print('End of compute_best. Hit RETURN to continue')
    input()


def compute_best_hv(probs, resdir, maxfun, solvername=None, solvers=None, refs=None):
    '''
    estrae da N prove di NSGA2 il "best" e il "worst" e crea dei risultati
    come se questi due fossero nuovi solutori
    '''

    if (refs == None):
        refs = solvers

    randsolvers = solvers

    for (n, ctype, prob_id, prbsel) in probs:
        probl = problem_factory(prob_id=prob_id, n=n, ctype=ctype, nint=nint, prbsel=prbsel)
        print(probl.prob.name, end='')

        # first thing: compute the Pareto reference front for current problem
        prob_front = [None] * (len(refs))
        for (s,solver_name) in enumerate(refs):
            if ('NOMAD' in solver_name) and (probl.prob.q > 2):
                prob_front[s] = np.empty((probl.prob.q,0))
            else:
                prob_front_ = sio.loadmat("{}/{}/{}.mat".format(resdir, solver_name, probl.prob.name))
                npoints, nobj = prob_front_['F'].shape
                prob_front[s] = np.transpose(prob_front_['F'])
                if npoints == 0:
                    prob_front[s] = np.reshape(prob_front[s], (probl.prob.q, 0))

        F_union = np.concatenate(prob_front, axis=1)
        _, n_points = F_union.shape
        if (True):
            if n_points == 0:
                F_true = np.copy(F_union)
            elif n_points == 1:
                F_true = np.copy(F_union)
            else:
                # print('ATT(2):',n_points,F_union.shape)
                efficient_indices = pareto_efficient(F_union)
                F_true = F_union[:, efficient_indices]
        else:
            F_true = np.copy(F_union)

        nobj, npoints = F_true.shape
        if npoints > 0:
            ref_point = np.max(F_true, axis=1)
            idl_point = np.min(F_true, axis=1)
        else:
            ref_point = np.inf * np.ones(nobj)
            idl_point = np.inf * np.ones(nobj)

        prob_front = [None] * (len(randsolvers))
        HV = np.zeros((maxfun, len(randsolvers)))

        for (s, solver_name) in enumerate(randsolvers):
            # print(s,' ',solver_name)
            # print(' ',s,end='')
            prob_front_ = sio.loadmat("{}/{}/{}.mat".format(resdir, solver_name, probl.prob.name))

            npoints, nobj = prob_front_['F'].shape
            prob_front[s] = np.transpose(prob_front_['F'])
            HV[:, s] = (np.transpose(prob_front_['HV'])[1:maxfun+1, 0])
            HV[np.where(HV[:,s] <= 0), s] = np.nan
            #print(HV[:,s])
            if npoints == 0:
                prob_front[s] = np.reshape(prob_front[s], (probl.prob.q, 0))
            # print(prob_front[s].shape,solver_name,s)

        #input()
        # print()
        # print(probl.prob.name,end='')

        ratios   = []
        ratios   = [None] * n_solvers
        hypervol = [None] * n_solvers

        for (s, solver_name) in enumerate(randsolvers):
            current_front = prob_front[s]
            _,n_points = current_front.shape
            if n_points > 0:
                #print(current_front.shape)
                _, index = np.unique(current_front, return_index=True, axis=1)
                index = sorted(index)
                current_front = current_front[:,index].T
                npoints, nobj = current_front.shape
                print(current_front.shape)
                #if npoints > 0:
                #   ref_point = np.max(current_front, axis=0)
                #   idl_point = np.min(current_front, axis=0)
                #else:
                #    ref_point = np.inf * np.ones(nobj)
                #    idl_point = np.inf * np.ones(nobj)
                Ftemp = front_factory(nobj, Fin=current_front)
                FT = Ftemp.get_clean_front(ref_point)
                FT = transform_front(FT, idl_point, ref_point, nobj)
                print(current_front.shape, FT.shape, len(FT))
                nadir = transform_point(ref_point, idl_point, ref_point, nobj)
                print(FT)
                print('nadir = ', nadir)
                try:
                    #print(ref_point.shape)
                    #hv = pg.hypervolume(current_front)
                    #ratios[s] = ( 1. / hv.compute(ref_point) )
                    if len(FT) > 0:
                        hv = pg.hypervolume(FT)
                        print(hv.compute(nadir))
                        ratios[s] = ( 1. / hv.compute(nadir) )
                    else:
                        ratios[s] = np.inf
                except ZeroDivisionError:
                    ratios[s] = np.inf
                except ValueError:
                    ratios[s] = np.inf
            else:
                ratios[s] = np.inf

        print(ratios)
        #input()
        ind = np.argsort(ratios)
        indmin = ind[0]
        indmax = ind[-1]
        media = np.sum(ratios) / len(randsolvers)
        media = 0.0
        for i,r in enumerate(ratios):
            if r < np.inf:
                media += r
        media /= len(randsolvers)
        indmed = -1
        dist = np.inf
        for i, r in enumerate(ratios):
            if np.abs(media - r) < dist:
                dist = np.abs(media - r)
                indmed = i
        #print(ratios)
        print(indmin,indmed,indmax)
        #print(' ratios: worst({}) med({}) best({})'.format(ratios[indmin],ratios[indmed],ratios[indmax]))
        print(' ratios: best({}) med({}) worst({})'.format(ratios[indmin],ratios[indmed],ratios[indmax]))
        #input()
        pathw = '{}/{}'.format(resdir, solvername + '_worsthv')
        pathm = '{}/{}'.format(resdir, solvername + '_medianhv')
        pathb = '{}/{}'.format(resdir, solvername + '_besthv')
        if not os.path.exists(pathw):
            os.mkdir(pathw)
        if not os.path.exists(pathm):
            os.mkdir(pathm)
        if not os.path.exists(pathb):
            os.mkdir(pathb)

        # se usa hypervol
        #sio.savemat('{}/{}.mat'.format(pathw, probl.prob.name), {'F': prob_front[indmin].T, 'HV': HV[:, indmin]})
        #sio.savemat('{}/{}.mat'.format(pathm, probl.prob.name), {'F': prob_front[indmed].T, 'HV': HV[:, indmed]})
        #sio.savemat('{}/{}.mat'.format(pathb, probl.prob.name), {'F': prob_front[indmax].T, 'HV': HV[:, indmax]})

        # se usa 1 / hypervol
        sio.savemat('{}/{}.mat'.format(pathw, probl.prob.name), {'F': prob_front[indmax].T, 'HV': HV[:, indmax]})
        sio.savemat('{}/{}.mat'.format(pathm, probl.prob.name), {'F': prob_front[indmed].T, 'HV': HV[:, indmed]})
        sio.savemat('{}/{}.mat'.format(pathb, probl.prob.name), {'F': prob_front[indmin].T, 'HV': HV[:, indmin]})
        # input()
        print()

    print('End of compute_best (HV). Hit RETURN to continue')
    input()


# -----------------------------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------------------------

# Argument parser
parser = argparse.ArgumentParser(description="Compute performance profiles given two or more algs")
# Positional arguments
parser.add_argument('--recompute',action='store_true',help='recompute using data in dir')
parser.add_argument('--compute_feas_stat',action='store_true',help='compute feasibility statistics ONLY')
parser.add_argument('--dir', type=str, help='Dir in which the results are stored', default=".")

parser.add_argument('--names', type=str, nargs="+", help='Names of the methods', default=False)

args = parser.parse_args(sys.argv[1:])

export = True
res_dir = args.dir
recompute = args.recompute
feas_stat = args.compute_feas_stat

if feas_stat:
    recompute = False
    FEASSTAT = {}

if not recompute:

    if False:
        solvers = ['NOMAD', 'NSGA2_best']
        solvers = ['DFMOINT', 'NSGA2_best', 'NSGA2_worst']
        solvers = ['DFMOINT', 'NSGA2_worst', 'NSGA2_median', 'NSGA2_best']
        solvers = ['DFMOINT', 'NOMAD']
        solvers = ['DFMOINT_versione_FORTRAN(True)', 'DFMOINT_versione_FORTRAN(False)']
        solvers = ['DFMOINT_versione_FORTRAN(False)', 'NOMAD']
        solvers = ['DFMOINT_versione_FORTRAN(True)', 'NSGA2_best', 'NSGA2_median', 'NSGA2_worst']
        solvers = ['NOMAD', 'NSGA2_median']
        solvers = ['DFMOINT_versione_FORTRAN(False)', 'NSGA2_median']
        solvers = ['NOMADcon', 'NSGA2con_best']
        solvers = ['NOMADcon', 'NSGA2con_median']
        solvers = ['DFMOINT', 'NSGA2con_median']
        solvers = ['DFMOINT', 'NSGA2con_best']
        solvers = ['DFMOINT', 'DFMOINT_2022-03-28']
        solvers = ['DFMOINT_versione_FORTRAN(False)', 'DFMOINT_2022-03-30']

        solvers = ['DFMOINT_versione_FORTRAN(False)', 'NSGA2_best', 'NSGA2_median', 'NSGA2_worst']
        solvers = ['DFMOINT_versione_FORTRAN(False)', 'NOMAD']
        solvers = ['DFMOINT', 'NSGA2con_best', 'NSGA2con_median', 'NSGA2con_worst']
        solvers = ['DFMOINT', 'NOMADcon']
        solvers = ['DFMOINT_2022-03-30', 'NSGA2_best', 'NSGA2_median', 'NSGA2_worst']
        solvers = ['DFMOINT_2022-03-30', 'NOMAD_noNM']
        solvers = ['DFMOINT_2022-03-30', 'NOMADnomodel']
        solvers = ['DFMOINT_2022-03-30', 'NOMAD']
        solvers = ['DFMOINTcon_2022-03-28', 'NOMADcon']
        solvers = ['DFMOINTcon_2022-03-28', 'NSGA2con_best', 'NSGA2con_median', 'NSGA2con_worst']
        solvers = ['DFMOINT', 'NSGA2_worst', 'NSGA2_median', 'NSGA2_best']
        solvers = ['DFMOINT', 'NOMAD']

        solvers = ['DFMOINTcon', 'NOMADcon', 'NSGA2con_best', 'NSGA2con_median', 'NSGA2con_worst']
        solvers = ['DFMOINTcon', 'NSGA2con_1', 'NSGA2con_2', 'NSGA2con_3', 'NSGA2con_4', 'NSGA2con_5', 'NSGA2con_7', \
                   'NSGA2con_8', 'NSGA2con_9', 'NSGA2con_10']
        solvers = ['DFMOINT', 'NOMAD', 'NSGA2_1', 'NSGA2_2', 'NSGA2_3', 'NSGA2_4', 'NSGA2_5', 'NSGA2_7', \
                   'NSGA2_8', 'NSGA2_9', 'NSGA2_10']
        solvers = ['DFMOINT', 'NSGA2_best', 'NSGA2_median', 'NSGA2_worst']
        solvers = ['DFMOINTcon', 'NOMADcon']
        solvers = ['DFMOINTcon', 'NSGA2con_best', 'NSGA2con_median', 'NSGA2con_worst']
        solvers = ['NSGA2con_1', 'NSGA2con_2', 'NSGA2con_3', 'NSGA2con_4', 'NSGA2con_5', 'NSGA2con_6', 'NSGA2con_7', \
                   'NSGA2con_8', 'NSGA2con_9', 'NSGA2con_10']
        #solvers = ['NSGA2con_1', 'NSGA2con_2', 'NSGA2con_3', 'NSGA2con_4']
        solvers = ['DFMOINT', 'NOMAD']

        solvers = ['NSGA2_10739', 'NSGA2_13578', 'NSGA2_22656', 'NSGA2_29648', 'NSGA2_49386', \
                   'NSGA2_109713', 'NSGA2_18387', 'NSGA2_287197', 'NSGA2_2980', 'NSGA2_67565']

####################################################
#  number of problems:
#  Unconstrained   q<=2 |    Constrained  q<=2 |  Feasible   q<=2
#-----------------------|----------------------|-------------------
#  total    96     58   |    576          348  |  401        213
#  n <= 12  44     19   |    264          114  |  197         84
#  n >= 13  52     39   |    312          234  |  204        129
    comp_best = False
    constrained = True
    onlyfeasible = True
    dimensions = 'big'  #must be one of {'all', 'big', 'small'}

    maxfun = 20000
    gates = [1.e-2, 5.e-2, 1.e-1, 5.e-1]
    prob_dims = [-np.inf, np.inf]
    if dimensions == 'big':
        prob_dims = [13, np.inf]
    elif dimensions == 'small':
        prob_dims = [1, 12]

    if comp_best:
        solvers = ['NSGA2_10739', 'NSGA2_13578', 'NSGA2_22656', 'NSGA2_29648', 'NSGA2_49386', \
                   'NSGA2_109713', 'NSGA2_18387', 'NSGA2_287197', 'NSGA2_2980', 'NSGA2_67565']
    else:
        solvers = ['DFMOINT_13578', 'NSGA2_best', 'NSGA2_median', 'NSGA2_worst']
        solvers = ['DFMOINT_13578', 'NOMAD_13578']
        solvers = ['DFMOINT_13578', 'NOMAD_13578']
        solvers = ['DFMOINT_13578', 'NOMAD_13578', 'NSGA2_10739',   'NSGA2_13578',  'NSGA2_22656', \
                   'NSGA2_29648',  'NSGA2_49386', 'NSGA2_109713',  'NSGA2_18387',  'NSGA2_287197', \
                   'NSGA2_2980',   'NSGA2_67565']
        solvers = ['DFMOINT_13578', 'NOMAD_13578', 'NSGA2_10739',   'NSGA2_13578',  'NSGA2_22656', \
                   'NSGA2_29648',  'NSGA2_49386', 'NSGA2_109713',  'NSGA2_18387',  'NSGA2_287197', \
                   'NSGA2_2980',   'NSGA2_67565']
        solvers = ['DFMOINT_13578', 'NSGA2_10739',   'NSGA2_13578',  'NSGA2_22656', \
                   'NSGA2_29648',  'NSGA2_49386', 'NSGA2_109713',  'NSGA2_18387',  'NSGA2_287197', \
                   'NSGA2_2980',   'NSGA2_67565']
        solvers = ['DFMOINT_13578', 'NSGA2_besthv', 'NSGA2_medianhv', 'NSGA2_worsthv']
        solvers = ['DFMOINT_13578', 'NOMAD_13578', 'NSGA2_besthv', 'NSGA2_medianhv', 'NSGA2_worsthv']

    solvers_for_ref = ['DFMOINT_13578', 'NOMAD_13578']
    solvers_for_ref = ['DFMOINT_13578', 'NSGA2_10739',   'NSGA2_13578',  'NSGA2_22656', \
                       'NSGA2_29648',  'NSGA2_49386', 'NSGA2_109713',  'NSGA2_18387',  'NSGA2_287197', \
                       'NSGA2_2980',   'NSGA2_67565']
    solvers_for_ref = ['DFMOINT_13578', 'NOMAD_13578']
    solvers_for_ref = ['DFMOINT_13578', 'NSGA2_10739',   'NSGA2_13578',  'NSGA2_22656', \
                       'NSGA2_29648',  'NSGA2_49386', 'NSGA2_109713',  'NSGA2_18387',  'NSGA2_287197', \
                       'NSGA2_2980',   'NSGA2_67565']
    solvers_for_ref = ['DFMOINT_13578', 'NOMAD_13578', 'NSGA2_10739',   'NSGA2_13578',  'NSGA2_22656', \
                       'NSGA2_29648',  'NSGA2_49386', 'NSGA2_109713',  'NSGA2_18387',  'NSGA2_287197', \
                       'NSGA2_2980',   'NSGA2_67565']

    #resdir  = '2022-03-23'
    #resdir  = '2022-11-10'
    #resdir  = '2022-12-22'
    if constrained:
        resdir  = '2022-12-24con'
        resdir  = '2024-01-28con'
        resdir  = '2024-02-14_DFMOINT_NOMAD_con'
        resdir  = '2024-02-14_DFMOINT_NSGA2_con'

        refdir  = '2024-02-27con'
        resdir  = '2024-03CONRES20000'

    else:
        resdir  = '2022-12-15'
        resdir  = '2024-02-12_NOMAD_DFMOINT'
        resdir  = '2024-02-14_DFMOINT_NOMAD_unc'
        resdir  = '2024-02-14_DFMOINT_NSGA2_unc'

        refdir  = '2024-02-27unc'
        resdir  = '2024-03UNCRES20000'

    if feas_stat:
        for s in solvers:
            FEASSTAT[s] = 0

    metrics = ["Purity", "Spread Gamma", "Spread Delta", "Hypervolume"]
    var_range = range(5, 51, 5)
    n_solvers = len(solvers)

    nint    = -1
    nvars   = [10,15,20,25,30]
    prids   = [1001, 1002, 1003, 1004, 1005, 1006, 1007]
    prids   = [1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010]
    if constrained:
        ctypes  = ['a','b','c','d','e','f']
    else:
        ctypes  = ['z']

    probs_cec   = [(n,c,id,'bk1') for n in nvars for c in ctypes for id in prids]

    TEMPPRB = problem_factory(prob_id=2000)
    LISTDMS = TEMPPRB.prob.list.list
    #LISTDMS = ['dpam1', 'fes1', 'l2zdt3', 'l2zdt4', 'l3zdt1', 'l3zdt2', 'l3zdt3', \
    #           'l3zdt4', 'mop2', 'qv1', 'sk2', 'zdt2', 'zdt3', 'zdt4', 'zdt6']

    LISTUTILI = []
    ndms = []
    qdms = []
    ndms_min = np.inf
    ndms_max = 0
    qdms_min = np.inf
    qdms_max = 0
    for p in LISTDMS:
        TEMPPRB = problem_factory(prob_id=2000,prbsel=p)
        n1,nc1,ni1,m1,q1 = TEMPPRB.prob.getdim()
        n1,q1 = TEMPPRB.prob.prb.getdim()
        if not (n1 in ndms):
            ndms.append(n1)
        if not (q1 in qdms):
            qdms.append(q1)
        if n1 < ndms_min:
            ndms_min = n1
        if n1 > ndms_max:
            ndms_max = n1
        if q1 < qdms_min:
            qdms_min = q1
        if q1 > qdms_max:
            qdms_max = q1
            
        #print(p," ",TEMPPRB.prob.getdim()," -- ",TEMPPRB.prob.prb.getdim())
    	#print(p," ",TEMPPRB.prob.getdim()," -- ",TEMPPRB.prob.prb.getdim(), \
    	#			TEMPPRB.prob.functs(TEMPPRB.prob.startp()),TEMPPRB.prob.startp())
        if n1 >= 4:
            LISTUTILI.append(p)

    probs_dms   = [(-1,c,2000,p) for c in ctypes for p in LISTUTILI]

    probs = probs_cec
    probs = probs_dms
    probs = probs_cec + probs_dms

    mdms = []
    for (n,ctype,prob_id,prbsel) in probs:
        TEMPPRB = problem_factory(prob_id=prob_id,n=n,ctype=ctype,prbsel=prbsel)
        n1,nc1,ni1,m1,q1 = TEMPPRB.prob.getdim()
        if not (m1 in mdms):
            mdms.append(m1)
        if not feas_stat:
            if prob_id == 2000:
                print(prbsel," ",ctype," ",TEMPPRB.prob.getdim()," -- ",TEMPPRB.prob.prb.getdim())
            else:
                print(prbsel," ",ctype," ",TEMPPRB.prob.getdim())
    
    print("problemi DMS: n \in [{},{}], q\in [{},{}]\n".format(ndms_min,ndms_max,qdms_min,qdms_max))
    print(ndms,qdms,mdms)
    print(len(probs))
    #print(probs)
    #input()

    count_problems = 0
    for (n,ctype,prob_id,prbsel) in probs:
        probl = problem_factory(prob_id=prob_id,n=n,ctype=ctype,nint=nint,prbsel=prbsel)
        nn,_,_,_,q = probl.prob.getdim()
        #print(nn)
        if (('NOMAD' in solvers) or \
            ('NOMADcon' in solvers) or \
            ('NOMAD_noNM' in solvers) or \
            ('NOMAD_13578' in solvers) or \
            ('NOMADnomodel' in solvers)) and (q > 2):
            continue
        if not (nn >= prob_dims[0] and nn <= prob_dims[1]):
            #print(prob_dims[0],n,prob_dims[1])
            continue
        #if n <= 12:
        #    continue
        print(prob_dims[0], nn, prob_dims[1],' q = ',probl.prob.q)

        # calcola il fronte approssimato del problema corrente F_true
        prob_front = [None] * len(solvers_for_ref)
        for (s, solver_name) in enumerate(solvers_for_ref):
            if ('NOMAD' in solver_name) and (probl.prob.q > 2):
                prob_front[s] = np.empty((probl.prob.q,0))
            else:
                prob_front_ = sio.loadmat("{}/{}/{}.mat".format(resdir,solver_name, probl.prob.name))
                npoints,nobj = prob_front_['F'].shape
                prob_front[s] = np.transpose(prob_front_['F'])
                if npoints == 0:
                    prob_front[s] = np.reshape(prob_front[s],(probl.prob.q,0))

        lenpf = [None] * len(solvers_for_ref)
        for i,pf in enumerate(prob_front):
            nobj, lenpf[i] = pf.shape

        if np.max(lenpf) > 0 or not onlyfeasible:
            count_problems += 1

    path_name = "_".join(solvers)
    if len(path_name) > 50:
        path_name = path_name[:50]
        path_name = path_name + "cont"
    print(path_name)
    input()

    #path_name = "2OBJ_" + path_name

    date = datetime.today().strftime("%Y-%m-%d_%H%M%S")
    datesuffix = "{}_".format(date)
    path_name_date = datesuffix + path_name

    print()
    print('-----------------------------------------------------')
    print('computing cumulative results')
    print('          maxfun: ',maxfun)
    print('     for solvers: ',solvers)
    print('with solvers_ref: ',solvers_for_ref)
    print('   using metrics: ',metrics)
    print('using results in: ',resdir)
    print('   only feasible: ==> ',onlyfeasible,' <==')
    print('for problems defined by:')
    print('\t         nvars = ',nvars)
    print('\t        ctypes = ',ctypes)
    print('\t         prids = ',prids)
    print('\t problem count = ',count_problems,end='')
    if onlyfeasible:
        print(' feasible problems')
    else:
        print()
    print('-----------------------------------------------------')
    print()
    print('plots will be put in :',path_name_date)
    print()
    print('-----------------------------------------------------')
    print('If this is ok, hit return to continue',end='')
    input()
    ####################################################
    # computebest per NSGA2 va chiamato una sola volta
    # perchè scrive i risultati su file
    ####################################################
    if comp_best:
        #compute_best(probs,resdir,maxfun,solvers=solvers,solvername='NSGA2')
        compute_best_hv(probs,resdir,maxfun,solvers=solvers,solvername='NSGA2',refs=solvers_for_ref)
        print('Done compute_best')
        input()

    A = 0 # n.prob. where no solvers find feasible points
    B = 0 # n.prob. where at least a solver finds feasible points
    C = 0 # n.prob. where all solvers find feasible points

    k = 0
    p = -1
    ns = len(solvers)
    debug = False
    N = []
    H = np.zeros((maxfun,count_problems,ns))
    matrix = {metric_name: np.Inf + np.ones(shape=(count_problems, ns)) for metric_name in metrics}
    gate = 1.e-1
    #count_problems = 0
    for (n,ctype,prob_id,prbsel) in probs:
        probl = problem_factory(prob_id=prob_id,n=n,ctype=ctype,nint=nint,prbsel=prbsel)
        nn,_,_,_,q = probl.prob.getdim()

        if (('NOMAD' in solvers) or \
            ('NOMADcon' in solvers) or \
            ('NOMAD_noNM' in solvers) or \
            ('NOMAD_13578' in solvers) or \
            ('NOMADnomodel' in solvers)) and (q > 2):
            continue
        if not (nn >= prob_dims[0] and nn <= prob_dims[1]):
            #print(prob_dims[0],n,prob_dims[1])
            continue
        #if n <= 12:
        #    continue

        print('\nProblem name: {}  nobj: {}'.format(probl.prob.name, q))

        p += 1
        #if not feas_stat:
        #    print(probl.prob.name,end='')

        # calcola il fronte approssimato del problema corrente F_true
        print('\tcomputing approximate True Pareto front')
        print('\t---------------------------------------')
        prob_front = [None] * len(solvers_for_ref)
        for (s, solver_name) in enumerate(solvers_for_ref):
            if ('NOMAD' in solver_name) and (probl.prob.q > 2):
                prob_front[s] = np.empty((probl.prob.q,0))
            else:
                prob_front_ = sio.loadmat("{}/{}/{}.mat".format(resdir,solver_name, probl.prob.name))
                npoints,nobj = prob_front_['F'].shape
                prob_front[s] = np.transpose(prob_front_['F'])
                if npoints == 0:
                    prob_front[s] = np.reshape(prob_front[s],(probl.prob.q,0))

        lenpf = [None] * len(solvers_for_ref)
        for i,pf in enumerate(prob_front):
            nobj, lenpf[i] = pf.shape

        if np.max(lenpf) > 0 or not onlyfeasible:
            N.append(probl.prob.n + 1)
            #count_problems += 1

            # calcola il fronte approssimato F_true
            F_union = np.concatenate(prob_front, axis=1)
            _, n_points = F_union.shape
            if n_points == 0:
                F_true = np.copy(F_union)
            elif n_points == 1:
                F_true = np.copy(F_union)
            else:
                #print('ATT(2):',n_points,F_union.shape)
                #print('ATT(3):',F_union)
                efficient_indices = pareto_efficient(F_union)
                F_true = F_union[:, efficient_indices]

            # calcola l'ipervolume normalizzato di F_true
            _, npoints = F_true.shape
            if npoints > 0:
                nadir = np.max(F_true, axis=1)
                ideal = np.min(F_true, axis=1)
            else:
                nadir = np.inf * np.ones(probl.prob.q)
                ideal = np.inf * np.ones(probl.prob.q)

            print('\tFtrue(shape): {}  F1(shape): {}'.format(F_true.shape,prob_front[0].shape),end='')
            Ftemp = front_factory(probl.prob.q, Fin=F_true.T)
            FT = Ftemp.get_clean_front(nadir)
            FT = transform_front(FT, ideal, nadir, probl.prob.q)
            print('  FT(shape): {}'.format(FT.shape),end='')
            if (npoints == 1 and len(FT) == 0):
                print('  clean front is EMPTY')
            else:
                print()
            #print('---------------------')
            nadir1 = transform_point(nadir, ideal, nadir, probl.prob.q)
            ideal1 = transform_point(ideal, ideal, nadir, probl.prob.q)
            if len(FT) > 0:
                #print(FT)
                #print(nadir1)
                try:
                    hv = pg.hypervolume(FT)
                    HVp = hv.compute(nadir1)
                except:
                    HVp = 0.0
            else:
                HVp = 0.0

            prob_front = [None] * len(solvers)
            for (s, solver_name) in enumerate(solvers):
                print('\textract results for {}'.format(solver_name),end='')
                #prob_front_ = sio.loadmat("2022-02-03_20000/{}/{}.mat".format(solver_name, probl.prob.name))
                prob_front_ = sio.loadmat("{}/{}/{}.mat".format(resdir,solver_name, probl.prob.name))
                #print(solver_name, probl.prob.name)

                npoints,nobj = prob_front_['F'].shape
                if debug:
                    print(' debug: ',np.transpose(prob_front_['HV'])[:,0].shape,p,count_problems)
                hv = (np.transpose(prob_front_['HV'])[1:maxfun+1,0]).astype('float64')
                hv[np.where(hv <= 0)] = np.nan
                print(' F({}).shape: {}  HV({}).shape: {}'.format(solver_name,prob_front_['F'].shape,solver_name,hv.shape))

                if (len(hv) < maxfun):
                    Q = ma.masked_array(hv, np.isnan(hv))
                    min_hv = Q.min(axis=0)
                    for i in range(maxfun-len(hv)):
                        hv = np.append(hv,min_hv)

                #H[:,p,s] = 1. / hv
                H[:,k,s] = hv / HVp
                #prob_front[s] = np.reshape(prob_front_['F'],(npoints,probl.prob.q))
                #prob_front[s] = np.transpose(np.reshape(prob_front_['F'],(npoints,probl.prob.q)))
                prob_front[s] = np.transpose(prob_front_['F'])
                if npoints == 0:
                    prob_front[s] = np.reshape(prob_front[s],(probl.prob.q,0))
                #if not feas_stat:
                #    print(prob_front[s], prob_front[s].shape)

                if ~debug:
                    if feas_stat:
                        if npoints > 0:
                            FEASSTAT[solver_name] += 1
                        #else:
                        #    print(probl.prob.name, end='')
                        #    print(' ',solver_name)

                    else:
                        if npoints > 1:
                            #print('ATT(1):',npoints,prob_front[s].shape)
                            efficient_indices = pareto_efficient(prob_front[s])
                            prob_front[s] = prob_front[s][:, efficient_indices]

            #lenpf = [None] * len(solvers)
            #for i,pf in enumerate(prob_front):
            #    #if not feas_stat:
            #    #    print(pf.shape,end='')
            #    nobj, lenpf[i] = pf.shape

            if feas_stat:
                if np.max(lenpf) <= 0:
                    A += 1
                    print(probl.prob.name,q,' A **********')
                else: #if np.max(lenpf) > 0:
                    B += 1
                    print(probl.prob.name,q,' B')
                if np.min(lenpf) > 0:
                    C += 1
                    print(probl.prob.name,q,' C')
            #print()

            if ~debug:
                if not feas_stat:
                    #if np.max(lenpf) > 0 or not onlyfeasible:
                    #if np.max(lenpf) > 0 or not onlyfeasible:
                    #print(prob_front)
                    print('\n\tcompute purity')
                    purity = compute_purity(*prob_front,Ftrue=F_true)
                    #print(purity)

                    matrix["Purity"][k, :] = purity

                    print('\tcompute spreads')
                    gamma, delta = compute_gamma_delta(*prob_front,Ftrue=F_true)

                    matrix["Spread Gamma"][k, :] = gamma
                    matrix["Spread Delta"][k, :] = delta

                    print('\tcompute hypervol')
                    hypervol = compute_hypervol(*prob_front,nadir1=nadir,ideal1=ideal,Ftrue=F_true)
                    matrix["Hypervolume"][k, :] = hypervol

            k += 1
        else:
            print('infeasibles')

    print(k,count_problems)
    input()

    if feas_stat:
        print('\n',FEASSTAT)
        print('no one:',A,' at least one:',B,' all:',C)

    if export and not feas_stat:
        pdfpagefile = None

        if not os.path.exists(path_name_date):
            os.mkdir(path_name_date)

        ##################### ########## SAVE DATA to pickle file
        DATA = {'solvers': solvers, \
                'probs' : probs, \
                'metrics': metrics, \
                'matrix':  matrix,  \
                }
        filenamepickle = "./{}/DATA_{}.npy".format(path_name_date,path_name)
        filepickle = open(filenamepickle,'wb')
        pickle.dump(DATA,filepickle)
        filepickle.close()
        ##################################################

        for m in metrics:

            if export:
                pdfpagefile = PdfPages("./{}/profiles_{}.pdf".format(path_name_date, m))
            if False:
                cumulative, tau = compute_profiles(matrix[m])
                print('matrix[m].shape = ',matrix[m].shape)
                plot_performance_profiles(cumulative, tau, names=solvers, metric_name=m, pdfpagefile=pdfpagefile)
                #plot_performance_profiles_new(cumulative, tau, names=solvers, metric_name=m, pdfpagefile=pdfpagefile)
            else:
                plot_performance_profiles_wild(matrix[m],m,names=solvers,pdfpagefile=pdfpagefile)
            if pdfpagefile is not None:
                pdfpagefile.close()

        if export:
            k=0
            for gate in gates:
                pdfpagefile = PdfPages("./{}/data_profile{}.pdf".format(path_name_date,k))
                plot_data_profiles(H,N,gate,names=solvers,pdfpagefile=pdfpagefile)
                k += 1
                pdfpagefile.close()

            #if pdfpagefile is not None:
            #    pdfpagefile.close()

else:
    print('will use data in {}'.format(res_dir))

    filenamepickle = ""
    for file in glob.glob('./{}/*.npy'.format(res_dir)):
        #print(file)
        filenamepickle = file

    print('will retrive data from pickle file: {}'.format(filenamepickle))
    print('If this is ok, hit return to continue',end='')
    input()

    #filenamepickle = "./{}/DATA_{}.npy".format(res_dir,path_name)
    filepickle = open(filenamepickle,'rb')
    DATA = pickle.load(filepickle)
    filepickle.close()

    solvers = DATA['solvers']
    probs = DATA['probs']
    metrics = DATA['metrics']
    matrix = DATA['matrix']

    onlyfeasible = True

    #count_problems = 0
    #for (n,ctype,prob_id) in probs:
    #	count_problems += 1

    print()
    print('-----------------------------------------------------')
    print('computing cumulative results')
    print('     for solvers: ',solvers)
    print('   using metrics: ',metrics)
    print('using results in: ',res_dir)
    print('   only feasible: ==> ',onlyfeasible,' <==')
    print('for problems defined by:')
    #print('\t         nvars = ',nvars)
    #print('\t        ctypes = ',ctypes)
    #print('\t         prids = ',prids)
    #print('\t problem count = ',count_problems)
    print('-----------------------------------------------------')
    print()
    print('plots will be put in :',res_dir)
    print()
    print('-----------------------------------------------------')
    print('If this is ok, hit return to continue',end='')
    input()

    pdfpagefile = None

    for m in metrics:

        pdfpagefile = PdfPages("./{}/profiles_{}.pdf".format(res_dir, m))

        cumulative, tau = compute_profiles(matrix[m])
        plot_performance_profiles(cumulative, tau, names=solvers, metric_name=m, pdfpagefile=pdfpagefile)

        if pdfpagefile is not None:
            pdfpagefile.close()
