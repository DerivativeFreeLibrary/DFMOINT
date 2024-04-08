import numpy as np
from utility import fast_non_dominated_filter, pareto_efficient
import matplotlib.pyplot as plt


def compute_purity(*args,Ftrue=None):
    assert len(args) >= 2
    n_solvers = len(args)

    # Merge fronts
    purity = [None] * n_solvers

    if (Ftrue is None):
        F_union = np.concatenate(args, axis=1)
        _,n_points = F_union.shape
        if n_points == 0:
            F_true = np.copy(F_union)
        elif n_points == 1:
            F_true = np.copy(F_union)
        else:
            #print('ATT(2):',n_points,F_union.shape)
            efficient_indices = pareto_efficient(F_union)
            F_true = F_union[:, efficient_indices]
    else:
        F_true = Ftrue

    #print('F_true.shape = ',F_true.shape)
    q,n_points = F_true.shape
    if q == 6:
        #ax = plt.axes(projection='3d')
        #ax.scatter3D(F_true[0,:],F_true[1,:],F_true[2,:],'.')
        plt.plot(F_true[0, :], F_true[1, :], '.')
        plt.show()
        input()

    for s in range(n_solvers):
        current_front = args[s]
        _,n_points = current_front.shape
        #print(current_front.shape)
        if n_points == 0:
            current_efficient = np.copy(current_front)
        else:
            _, index = np.unique(current_front, return_index=True, axis=1)
            index = sorted(index)
            current_front = current_front[:,index]
            n_points = current_front.shape[1]
            current_efficient = fast_non_dominated_filter(current_front, F_true, equal_is_good=True)
        try:
            purity[s] = n_points / len(np.where(current_efficient)[0])
        except ZeroDivisionError:
            purity[s] = np.nan
        #print(purity[s])
    #input()
    return purity
