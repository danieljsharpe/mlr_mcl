'''
Calculate the matrix exponential at a given lag time. For use with MLR-MCL algorithm C++ code.

Daniel J. Sharpe, April 2019
'''

import numpy as np
from scipy.linalg import expm

# calculate the exponential of the matrix passed in sparse format, i.e. lists of:
# non-zero elems, column indices for elems, cumulative row lengths
def calc_expm(A,A_ci,A_rl,tau):
    print "lag time is: %.6f" % tau
    K = np.zeros((len(A_rl),len(A_rl)),dtype=float) # transition rate matrix
    row_no = 0
    i = 0
    '''
    for a in A:
        K[A_ci[i],row_no] = np.exp(A[i])
        i += 1
        if i <= A_rl[row_no]: row_no += 1
    for i in range(np.shape(M)[0]):
        K[i,i] = -np.sum(M[i,:]) # set diagonal elements of transition rate matrix
    T = expm(tau*K);

    print "transition matrix:"
    print T
    '''
    print "lengths of arrays:"
    print len(A), len(A_ci), len(A_rl)
    print "first elems of arrays:"
    print A[0], A_ci[0], A_rl[0]
    print "goodbye"
    return 5000
