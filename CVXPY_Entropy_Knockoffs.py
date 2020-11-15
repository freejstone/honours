#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:53:39 2020

@author: jackfreestone
"""

import numpy as np
import cvxpy as cp


def optimal_entropy_cvxpy_batch(I, Sig, nsko, warm_start = None):
    if type(I) == int:
        lenI = 1
        I = [I]
    else:
        lenI = len(I)
    p = len(Sig)
    I_PP = np.diag(np.tile(1, p)) #identity matrix
    I_PI = I_PP[:, I]
    I_IP = I_PI.transpose()
    s0 = cp.Variable(lenI)
    if warm_start is None:
        Objective = cp.Minimize(-cp.log_det( ((nsko + 1)/nsko)*Sig - I_PI@cp.diag(s0)@I_IP) - nsko*cp.sum(cp.log(s0)))
        prob = cp.Problem(Objective, [s0 >= 0.00001, ((nsko + 1)/nsko)*Sig >> I_PI@cp.diag(s0)@I_IP])
        prob.solve()
        return s0.value
    else:
        if lenI == 1:
            s0.value = np.array([warm_start])
        else:
            s0.value = warm_start
        Objective = cp.Minimize(-cp.log_det( ((nsko + 1)/nsko)*Sig - I_PI@cp.diag(s0)@I_IP) - nsko*cp.sum(cp.log(s0)))
        prob = cp.Problem(Objective, [s0 >= 0, ((nsko + 1)/nsko)*Sig >> I_PI@cp.diag(s0)@I_IP])
        prob.solve(solver = cp.SCS)
        return s0.value