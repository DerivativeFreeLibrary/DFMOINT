import numpy as np
from utility import pareto_efficient


def compute_gamma_delta(*args,Ftrue=None):
    assert len(args) >= 2
    n_solvers = len(args)
    n_obj = args[0].shape[0]
    #
    # Merge Pareto front
    #
    spread_gamma, spread_delta = [None] * n_solvers, [None] * n_solvers

    if Ftrue is None:
        F_temp = np.concatenate(args, axis=1)
        _, n_points = F_temp.shape
        if n_points > 0:
            efficient_indices = pareto_efficient(F_temp)
            F_union = F_temp[:, efficient_indices]
    else:
        F_union = Ftrue

    _, n_points = F_union.shape
    if n_points > 0:
        # Compute extreme points of true front
        extreme_points = np.zeros((n_obj, 2))
        for j in range(n_obj):
            F_j = F_union[j, :]
            np.matrix.sort(F_j)
            #print(F_j.shape)
            extreme_points[j, 0] = F_j[0]
            extreme_points[j, 1] = F_j[-1]

        # Compute gamma, delta
        for s in range(n_solvers):
            delta_max = -1
            tsum_max = -1
            current_front = np.sort(args[s])
            for j in range(n_obj):
                F_j = current_front[j, :]
                n_points = len(F_j)
                if n_points >= 2:
                    delta = np.array([0] * (n_points + 1), dtype=float)
                    delta[0] = abs(F_j[0] - extreme_points[j, 0])

                    delta[1:n_points] = abs(F_j[1:] - F_j[:n_points - 1])
                    delta[n_points] = abs(extreme_points[j, 1] - F_j[-1])

                    if max(delta) > delta_max:
                        delta_max = max(delta)

                    delta_avg = np.mean(delta[1:n_points])
                    tsum = delta[0] + delta[-1] + sum(abs(delta[1:n_points] - delta_avg))
                    tsum = tsum / (delta[0] + delta[-1] + (n_points - 1) * delta_avg)
                else:
                    tsum = np.nan
                    delta_max = np.nan

                if tsum > tsum_max:
                    tsum_max = tsum

            if delta_max <= 0:
                delta_max = np.nan

            if tsum_max <= 0:
                tsum_max = np.nan

            spread_gamma[s] = delta_max
            spread_delta[s] = tsum_max
    else:
        for s in range(n_solvers):
            spread_gamma[s] = np.nan
            spread_delta[s] = np.nan

    return spread_gamma, spread_delta
