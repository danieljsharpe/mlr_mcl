'''
Calculate the matrix exponential at a given lag time. For use with MLR-MCL algorithm C++ code.

Daniel J. Sharpe, April 2019
'''

import numpy as np
from scipy.linalg import expm

''' calculate the exponential of the matrix passed in sparse format, i.e. lists of:
    non-zero elems, column indices for elems, cumulative row lengths '''
def calc_expm(A,A_ci,A_rl,tau):
    print "lag time is: %.6f" % tau
    A = np.array(A)
    A_ci = np.array(A_ci)
    A_rl = np.array(A_rl)
    K = np.zeros((len(A_rl),len(A_rl)),dtype=float) # transition rate matrix
    row_no = 0
    print "dimensions of K:", np.shape(K)
    print "lengths of arrays:"
    print len(A), len(A_ci), len(A_rl)
    for i in range(np.shape(A)[0]):
#        print "i:", i, "col_idx:", A_ci[i], "row_no:", A_rl[row_no], "    rl:", row_no
        K[A_ci[i],row_no] = A[i]
        while True:
            if i==A_rl[row_no]-1 and i<np.shape(A)[0]-1: row_no += 1
            else: break
    print "calculating diagonal elements"
    for i in range(np.shape(A_rl)[0]):
        K[i,i] = -np.sum(K[i,:]) # set diagonal elements of transition rate matrix
        print "row:", i, "sum of row:", np.sum(K[i,:])
    print "transition rate matrix:"
    print K
    print "taking matrix exponential"
    T = expm(tau*K);
    print "transition matrix:"
    print T
    for i in range(np.shape(A_rl)[0]):
        print "row:", i, "sum of row:", np.sum(T[i,:])
#        assert abs(np.sum(A[:,i])-1) < 1.0E-06, "Error, row %i of transition matrix does not sum to 1" % i

    print "transition matrix:"
    print T
    '''
    print "first elems of arrays:"
    print "%2.6f %5i %5i" % (A[0], A_ci[0], A_rl[0])
    '''
    return 5000
