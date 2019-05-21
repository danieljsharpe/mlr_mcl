'''
Calculate the matrix exponential at a given lag time. For use with MLR-MCL algorithm C++ code.

Daniel J. Sharpe, April 2019
'''

import numpy as np
from scipy.linalg import expm

''' calculate the exponential of the matrix passed in sparse format, i.e. lists of:
    non-zero elems, column indices for elems, cumulative row lengths '''
def calc_expm(A,A_ci,A_rl,tau,eps):
    print "lag time is: %E" % tau
    tau = np.float128(tau)
    A = np.array(A,dtype=np.float128)
    A_ci = np.array(A_ci,dtype=int)
    A_rl = np.array(A_rl,dtype=int)
    K = np.zeros((len(A_rl),len(A_rl)),dtype=np.float128) # transition rate matrix
    row_no = 0
    print "lengths of arrays:\n", np.shape(A)[0], np.shape(A_ci)[0], np.shape(A_rl)[0]
    for i in range(np.shape(A)[0]):
        K[A_ci[i],row_no] = tau*np.exp(A[i])
        while True:
            if i==A_rl[row_no]-1 and i<np.shape(A)[0]-1: row_no += 1
            else: break
    print "calculating diagonal elements"
    for i in range(np.shape(A_rl)[0]):
        assert K[i,i] == 0., "Error: row %i has a non-zero diagonal element (self-loop)" % i
        k_row = np.sort(K[i,:].copy())
        K[i,i] = -np.sum(k_row) # set diagonal elements of transition rate matrix
        assert abs(np.sum(np.sort(K[i,:].copy()))) < 1.0E-08, "Error: row %i of K does not sum to 0" % i
    print "transition rate matrix:"
    print K
    T = expm(K)
    print "transition matrix:"
    print T
    for i in range(np.shape(A_rl)[0]):
        assert abs(np.sum(np.sort(T[i,:].copy()))-1.) < 1.0E-04, "Error: row %i of T does not sum to 1" % i
    T = np.transpose(T) # transition matrix in MCL is column-stochastic
    # pass matrix back in CSR sparse format
    Tsp = []
    Tci = []
    Trl = np.zeros(np.shape(A_rl)[0],dtype=int)
    ci = 0 # col index
    rl = 0 # row length
    rn = 0 # row no.
    for t in T.flatten(order="C"):
        if ci==np.shape(A_rl)[0]: # new row
            Trl[rn] = rl
            ci = 0
            rl = 0
            rn += 1
        if t > eps:
            Tsp.append(t)
            Tci.append(ci)
            rl += 1
        ci += 1
    Trl[rn] = rl
    for i in range(1,np.shape(A_rl)[0]):
        Trl[i] += Trl[i-1]
    assert Trl[-1]==np.shape(Tsp)[0], "Error: row lengths not consistent with size of transition matrix"
    Trl = list(Trl)
    print "Number of non-zero elements of the transition matrix: %i" % np.shape(Tsp)[0]
    return Tsp, Tci, Trl
