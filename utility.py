import numpy as np
from filter_type import filter_elem

def pareto_efficient(costs):
    """
    :param costs: the list to be analyzed, it has to be a np.array of shape (n_obj, n_points)
    :return: boolean vector of efficient indices (True efficient, False not efficient)
    """

    n_obj, n_points = costs.shape
    efficient = [False] * n_points

    # Delete duplicates
    #print(costs)
    #input()
    _, index = np.unique(costs, return_index=True, axis=1)
    index = sorted(index)
    duplicates = [el for el in np.arange(n_points) if el not in index]
    indices = np.arange(n_points)

    for i in range(n_points):
        partial_ix = duplicates + [i]
        # Comparison with only unique points and different from the point itself
        partial_matrix = costs[:, np.delete(indices, partial_ix)]
        dominance_matrix = partial_matrix - np.reshape(costs[:, i], newshape=(n_obj, 1))
        is_dominated = (sum(dominance_matrix <= 0) == n_obj).any()
        if not is_dominated:
            efficient[i] = True

    return efficient

def fast_non_dominated_filter(curr_list: np.ndarray, new_list: np.ndarray, equal_is_good=False):
    """
    :param curr_list: list of non_dominated points
    :param new_list: new list of non_dominated solutions to be added
    :param equal_is_good: if True same pareto_points found in new_list and curr_list are mantained in curr_list else they are discarded
    :return: boolean vector of efficient indices (True efficient, False not efficient)
    """
    m, n_new_points = new_list.shape
    efficient = np.array([True] * curr_list.shape[1])
    neq_ix = None

    for i in range(n_new_points):
        dominance_matrix = curr_list - np.reshape(new_list[:, i], newshape=(m, 1))
        dominated_idx = (sum(dominance_matrix >= 0) == m)

        if equal_is_good:
            neq_ix = (sum(dominance_matrix == 0) != m)
            assert len(neq_ix.shape) == 1

        assert len(dominated_idx.shape) == 1
        dom_indices = np.where((dominated_idx * neq_ix if equal_is_good else dominated_idx) == True)[0]

        if len(dom_indices) > 0:
            efficient[dom_indices] = False

    return efficient

def minore_uguale(f1,f2):
	for i in range(len(f1)):
		if f1[i] > f2[i]:
			return False
	return True

def domina(f1,f2):
	for i in range(len(f1)):
		if f1[i] > f2[i]:
			return False

	for i in range(len(f1)):
		if f1[i] < f2[i]:
			return True

	return False

def maggiore_stretto(f1,f2):
	for i in range(len(f1)):
		if f1[i] <= f2[i]:
			return False
	return True

def insert_row(row,Lout):
	'''
	possono succedere TRE cose:
	1) Lkappa[i] < della tuple (x,f)
	2) (x,f) < Lkappa
	3) (x,f) non domina e non è dominato da Lkappa[i]
	'''
	flag = 0
	while flag < 3:
		if len(Lout) == 0:
			flag = 3
			continue
		for i, e in enumerate(Lout):
			#print('insert: ',i,e['fobs'],row['fobs'])
			if minore_uguale(e.fobs,row.fobs) and minore_uguale(row.fobs,e.fobs):
				flag = 1
				return False, Lout
			if(domina(e.fobs,row.fobs)):
				flag = 1
				return False, Lout
			if(domina(row.fobs,e.fobs)):
				flag = 2
				#print('pop',i)
				Lout = np.delete(Lout,i)
				break
			else:
				flag = 3

	Lout = np.append(Lout,filter_elem(row.x,row.fobs,row.alfa,row.alfa_c,row.alfa_tilde,row.xi,row.flag_allones))
	#for i,e in enumerate(Lout):
	#	print('insert(Lout):',i,e['fobs'])

	return True, Lout

def merge(Lin,Lout,iprint):
	flag_change = False
	if iprint > 1:
		print('---------- MERGE START ----------')
		for i,e in enumerate(Lin):
			print('Lin:',i,e.fobs)
		for i,e in enumerate(Lout):
			print('Lou:',i,e.fobs)

	for i,e in enumerate(Lin):
		flag, Lout = insert_row(e,Lout)
		flag_change = flag_change or flag
	if iprint > 1:
		for i,e in enumerate(Lout):
			print('Lou:',i,e.fobs)
		print('---------- MERGE END ----------')
	return Lout,flag_change
