# -*- coding: utf-8 -*-
"""
Created on Wed May 17 12:48:27 2017

@author: mihai
"""

import numpy as np
import scipy as sc

R1 = 1.25
R2 = 1.25
c0 = 0.6
c1 = sc.sqrt(R1**(-1) * (c0**2 + 5/6))
c2 = sc.sqrt(R2**(-1) * (c0**2 + 5/6))
cT = sc.sqrt(c0**2 / (1 + c0**2))
K = 0.5
M_A = 0
mode = 'kink'

def m1(W):
    return sc.sqrt(1 - W**2 / c1**2)

def m2(W):
    return sc.sqrt(1 - W**2 / c2**2)

def m0(W, M_A):
    return sc.sqrt((1 - (W - M_A)**2) * (c0**2 - (W - M_A)**2) / ((1 + c0**2) * (cT**2 - (W - M_A)**2)))

def disp_rel(W, K, M_A):
    return m0(W, M_A)**2 * W**4 + 1/R1 * 1/R2 * m1(W) * m2(W) * (1 - (W - M_A)**2)**2 - \
        0.5 * W**2 * m0(W, M_A) * (1 - (W - M_A)**2) * (1/R1 * m1(W) + 1/R2 * m2(W)) * \
        (sc.tanh(m0(W, M_A) * K) + sc.tanh(m0(W, M_A) * K)**(-1))

def lambda0(W, M_A):
    return 1j * (1 - (W - M_A)**2)/(m0(W, M_A) * (W - M_A))

def lambda1(W):
    return 1j * R1 * W / m1(W)

def A(mode, W, M_A):
    return (sc.cosh(m1(W) * K) - sc.sinh(m1(W) * K))**(-1) * \
        (sc.cosh(m0(W, M_A) * K) * B(mode, W, M_A) - sc.sinh(m0(W, M_A) * K) * C(mode, W, M_A)) * W / (W - M_A)

def B(mode, W, M_A):
    if mode == 'sausage':
        return C(mode, W, M_A) * (lambda0(W, M_A) * sc.cosh(m0(W, M_A) * K) * (W - M_A) - lambda1(W) * sc.sinh(m0(W, M_A) * K) * W) / \
            (lambda0(W, M_A) * sc.sinh(m0(W, M_A) * K) * (W - M_A) - lambda1(W) * sc.cosh(m0(W, M_A) * K) * W)
    else:
        return 0.1

def C(mode, W, M_A):
    if mode == 'sausage':
        return 0.1
    else:
        return B(mode, W, M_A) * (lambda0(W, M_A) * sc.sinh(m0(W, M_A) * K) * (W - M_A) - lambda1(W) * sc.cosh(m0(W, M_A) * K) * W) / \
            (lambda0(W, M_A) * sc.cosh(m0(W, M_A) * K) * (W - M_A) - lambda1(W) * sc.sinh(m0(W, M_A) * K) * W)

def D(mode, W, M_A):
    return (sc.cosh(m2(W) * K) - sc.sinh(m2(W) * K) )**(-1) * \
        (sc.cosh(m0(W, M_A) * K) * B(mode, W, M_A) + sc.sinh(m0(W, M_A) * K) * C(mode, W, M_A)) * W / (W - M_A)

def v(mode, x, W, K, M_A):
    if np.real(x) < -K:
        return A(mode, W, M_A) * (sc.cosh(m1(W) * x) + sc.sinh(m1(W) * x))
    elif np.abs(np.real(x)) <= K:
        return B(mode, W, M_A) * sc.cosh(m0(W, M_A) * x) + C(mode, W, M_A) * sc.sinh(m0(W, M_A) * x)
    elif np.real(x) > K:
        return D(mode, W, M_A) * (sc.cosh(m2(W) * x) - sc.sinh(m2(W) * x))

def displacement(mode, x, W, K, M_A):
    if np.real(x) < -K:
        return v(mode, x, W, K, M_A) / W
    elif np.abs(np.real(x)) <= K:
        return v(mode, x, W, K, M_A) / (W - M_A)
    elif np.real(x) > K:
        return v(mode, x, W, K, M_A) / W

