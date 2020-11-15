#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:55:53 2020

@author: jackfreestone
"""
import numpy as np
import cvxopt as cvx

def optimal_entropy_cvxopt_batch(I, Sig, nsko):
    if type(I) == int:
        lenI = 1
        I = [I]
    else:
        lenI = len(I)
    
    p = len(Sig)
    I_PP = np.diag(np.tile(1, p))
    I_PI = I_PP[:, I] 
    
    def F(x = None, z = None):
        if x is None: return (0, cvx.matrix(min(np.linalg.eigvals(Sig)), (lenI, 1)))
        mat = ((nsko + 1)/nsko)*Sig - I_PI@np.diag(list(x))@I_PI.transpose()
        f = -cvx.log(np.linalg.det(mat)) - nsko*cvx.sum(cvx.log(x)) 
        mat_inv = cvx.matrix(np.linalg.inv(mat))
        Df = mat_inv[cvx.matrix(I)*p + cvx.matrix(I)] - cvx.div(nsko, x) 
        if z is None: return f, Df.T
        Matrix_list = []
        for i in range(lenI):
            I0 = cvx.matrix(0, (p, p))
            I0[I[i], I[i]] = 1
            A = -mat_inv*I0*mat_inv
            Matrix_list.append(cvx.matrix(A[cvx.matrix(I)*p + cvx.matrix(I)], (1, lenI)))
        DDf = -cvx.matrix(Matrix_list, (lenI, lenI))
        H = cvx.sparse(DDf + cvx.spdiag(cvx.div(nsko, x**2)))
        return f, Df.T, z[0]*H
    
    G = []
    for i in range(lenI):
        I0 = [0]*(p**2)
        I0[I[i]*p + I[i]] = 1.0
        I1 = [0]*(lenI)
        I1[i] = 1.0
        I2 = [0]*(lenI)
        I2[i] = -1.0
        G.append(I1 + I2 + I0)
    G = cvx.matrix(G)
    
    
    Sig1 = np.reshape(Sig, (1, p**2))
    I2 = np.array([0]*(lenI))
    I3 = np.array([1]*(lenI))
    h = cvx.matrix(np.concatenate((I3, I2,((nsko + 1)/nsko)*Sig1), axis = None))
    
    dims = {'l': 2*lenI, 'q': [], 's': [p]}
    cvx.solvers.options['show_progress'] = False
    sol = cvx.solvers.cp(F, G, h, dims)
    return np.array(sol['x'])
