import numpy as np
import pygmo as pg
from utility import fast_non_dominated_filter, pareto_efficient
from problem_factory import transform_front, transform_point
from filter_factory import front_factory

def compute_hypervol(*args,nadir1=None,ideal1=None,Ftrue=None):
    assert len(args) >= 2
    n_solvers = len(args)

    # Merge fronts
    hypervol = [None] * n_solvers

    if (Ftrue is None):
        F_union = np.concatenate(args, axis=1)
        _,n_points = F_union.shape
        if(True):
            if n_points == 0:
                F_true = np.copy(F_union)
            elif n_points == 1:
                F_true = np.copy(F_union)
            else:
                #print('ATT(2):',n_points,F_union.shape)
                efficient_indices = pareto_efficient(F_union)
                F_true = F_union[:, efficient_indices]
        else:
            F_true = np.copy(F_union)
    else:
        F_true = Ftrue

    nobj, npoints = F_true.shape
    if npoints > 0:
        ref_point = np.max(F_true, axis=1)
        idl_point = np.min(F_true, axis=1)
    else:
        ref_point = np.inf * np.ones(nobj)
        idl_point = np.inf * np.ones(nobj)

    #print(F_true.shape)
    #print(ref_point)
    #print(idl_point)
    #print(ref_point.shape)
    #input()
    for s in range(n_solvers):
        current_front = args[s]
        _,n_points = current_front.shape
        if n_points > 0:
            #print(current_front.shape)
            _, index = np.unique(current_front, return_index=True, axis=1)
            index = sorted(index)
            current_front = current_front[:,index].T
            npoints, nobj = current_front.shape
            #print(current_front.shape)
            #print(current_front.shape)
            #print(ref_point.shape)
            #hv = pg.hypervolume(current_front)
            #hypervol[s] = 1. / hv.compute(ref_point)
            ##FT = transform_front(current_front, idl_point, ref_point, nobj)
            ##nadir = transform_point(ref_point, idl_point, ref_point, nobj)
            Ftemp = front_factory(nobj,Fin=current_front)
            FT = Ftemp.get_clean_front(nadir1)
            FT = transform_front(FT, ideal1, nadir1, nobj)
            nadir = transform_point(nadir1, ideal1, nadir1, nobj)
            try:
                if len(FT) > 0:
                    hv = pg.hypervolume(FT)
                #    #print(hv.compute(nadir))
                    hypervol[s] = 1. / hv.compute(nadir)
                else:
                    hypervol[s] = np.nan
            except ZeroDivisionError:
                hypervol[s] = np.nan
            except ValueError:
                hypervol[s] = np.nan
        else:
            hypervol[s] = np.nan
        #print(hypervol[s])
    #print(hypervol)
    #input()
    return hypervol
