# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 13:17:02 2016

@author: app13jfm
"""

from __future__ import division
#import matplotlib
import numpy as np
#import matplotlib.cm as cm
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import math as mt
import mpmath as mp
import scipy.optimize as sp
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.axes_grid1 import host_subplot
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d.axes3d import Axes3D
#from matplotlib import rc
#import pdb
#import pickle
gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
mp.dps=10
mp.pretty= True

import mpmath as mp


d = 2
theta = 0.1

g=5/6

G=10

#define gammas
def gamma_0_a0(W):
    return (-G * W**2 + np.sqrt( G**2 * W**4 - 4*(W**2 - b) * (W**2 - 1) * ((1+b) * W**2 - b) +0j)) / (2 * ((1+b) * W**2 - b) + G * W**2)
def gamma_1_a0(W):
    g0 = gamma_0_a0(W)
    return 1j * (-2 * b * g0**3 + G * g0**2 + 2 * b * g0 + G) / (2 * ((1+b) * W**2 - b) * g0 + G * W**2)
def gamma_a0(W):
    return gamma_0_a0(W) + theta * gamma_1_a0(W)

def gamma_0_b0(W):
    return (-G * W**2 - np.sqrt(G**2 * W**4 - 4 * (W**2 - b) * (W**2 - 1) * ((1+b) * W**2 - b)+0j)) / (2 * ((1+b) * W**2 - b) + G * W**2)
def gamma_1_b0(W):
    g0 = gamma_0_b0(W)
    return 1j * (-2 * b * g0**3 + G * g0**2 + 2 * b * g0 + G) / (2 * ((1+b) * W**2 - b) * g0 + G * W**2)
def gamma_b0(W):
    return gamma_0_b0(W) + theta*gamma_1_b0(W)

def gamma_0_ae(W):
    return (-G * (W**2)/d + np.sqrt(G**2 * (W**2 / d)**2 - 4 * ((W**2 / d)-b)*(W**2-1) * ((1+b) * (W**2 / d) - b) +0j)) / (2 * ((1+b) * (W**2 / d) - b) + G * (W**2 / d))
def gamma_1_ae(W):
    g0 = gamma_0_ae(W)
    return 1j * (-2 * b * g0**3 + G * g0**2 + 2 * b * g0 + G) / (2 * ((1+b) * (W**2 / d) - b) * g0 + G * (W**2 / d))
def gamma_ae(W):
    return gamma_0_ae(W) + theta*gamma_1_ae(W)

def gamma_0_be(W):
    return (-G * (W**2) / d - np.sqrt(G**2 * (W**4 / d**2) - 4 * ((W**2 / d) - b) * (W**2 - 1) * ((1+b) * (W**2 / d) - b) +0j)) / (2 * ((1+b) * (W**2 / d) - b) + G * (W**2 / d))
def gamma_1_be(W):
    g0 = gamma_0_be(W)
    return 1j * (-2 * b * g0**3 + G * g0**2 + 2 * b * g0 + G) / (2 * ((1+b) * (W**2 / d) - b) * g0 + G * (W**2 / d))
def gamma_be(W):
    return gamma_0_be(W) + theta*gamma_1_be(W)

#define ratios between vx and vz values of exponentials    
def a0(W):
    ga = gamma_a0(W)
    return (theta * (1-ga**2) + 1j*g*b*ga) / (W**2 - b)
    
def b0(W):
    gb = gamma_b0(W)
    return (theta * (1-gb**2) + 1j*g*b*gb) / (W**2 - b)
    
def ae(W):
    ga = gamma_ae(W)
    return (theta * (1-ga**2) + 1j*g*b*ga) / ((W**2/d) - b)

def be(W):
    gb = gamma_be(W)
    return (theta * (1-gb**2) + 1j*g*b*gb) / ((W**2/d) - b)    

#critical speeds
def vA(b):
    return 1
def vAe(b):
    return np.sqrt(d)
def cT(b):
    return np.sqrt(g*b/ (1+g*b))
def cTe(b):
    return np.sqrt(d * g*b / (1+g*b))
def c0(b):
    return np.sqrt(g*b)
def ce(b):
    return np.sqrt(g*b*d)

#define all of the dispersionrelations!!!

# dispersion relations with 2 unknowns
#aeb0
def disp_rel_0_2_aeb0(W):
        return 1j*b*(gamma_0_ae(W) / (W**2 +b) + d * gamma_0_b0(W)) / (W**2 -d * b)

def disp_rel_1_2_aeb0(W):
        return (1j * b * gamma_1_b0(W) - gamma_0_b0(W)**2 + 1) / (W**2 -b) +\
                d * (1j * b * gamma_1_ae(W) - gamma_0_ae(W)**2 + 1) / (W**2 - d*b) 

def deriv_dr_0_2_aeb0(W):
    return 1j * W * (d ** 2 * (-2 * ((1 + b) ** 3 * W ** 12 - b * (b + 2) * (1 + b) ** 2 * (d + 1) * W ** 10 + b ** 2 * (d ** 2 + 1) * (1 + b) * (b ** 2 + 5 * b + 5) * W ** 8 / 2 - b ** 3 * (d + 1) * (d - 1) ** 2 * (b + 2) * (1 + b) * W ** 6 - (b ** 2 * d ** 2 + (-d ** 4 / 2 + 6 * d ** 2 - 1 / 2) * b - d ** 4 / 2 + 6 * d ** 2 - 1 / 2) * b ** 4 * W ** 4 + b ** 5 * d ** 2 * (d + 1) * (b + 2) * W ** 2 - b ** 6 * d ** 2 * (d ** 2 + 1) / 2) * G * np.sqrt((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) + (W ** 2 - b * d) ** 2 * ((1 + b) * W ** 2 - b * d) ** 2 * (-2 * (1 + b) ** 2 * W ** 8 + (1 + b) * (G ** 2 + 2 * b ** 2 + 8 * b + 4) * W ** 6 + (-6 * b ** 3 - 18 * b ** 2 - 12 * b) * W ** 4 - b ** 2 * (G ** 2 - 8 * b - 12) * W ** 2 - 4 * b ** 3)) * np.sqrt((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) - (-2 * d * (1 + b) ** 2 * W ** 8 + (1 + b) * (2 * b ** 2 * d ** 2 + 8 * b * d ** 2 + G ** 2 + 4 * d ** 2) * W ** 6 - 6 * b * d ** 3 * (b + 2) * (1 + b) * W ** 4 - b ** 2 * d ** 2 * (-8 * b * d ** 2 + G ** 2 - 12 * d ** 2) * W ** 2 - 4 * b ** 3 * d ** 5) * np.sqrt((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) * (W ** 2 - b) ** 2 * ((1 + b) * W ** 2 - b) ** 2) * ((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) ** (-1 / 2) * b * ((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) ** (-1 / 2) / d ** 2 / (W ** 2 - b * d) ** 2 / (W ** 2 - b) ** 2 / ((1 + b) * W ** 2 - b) ** 2 / ((1 + b) * W ** 2 - b * d) ** 2
    
#aea0        
def disp_rel_0_2_aea0(W):
#    if np.real(gamma_a0(W)) > 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) < 0 and np.real(gamma_be(W)) > 0:
    if W < cTe(b) and cT(b) < W < min(c0(b), vA(b)):
        return 1j*b*(gamma_0_ae(W) / (W**2 +b) + d * gamma_0_a0(W)) / (W**2 -d * b)
    else:
        raise Exception('incorrect range')
def disp_rel_1_2_aea0(W):
#    if np.real(gamma_a0(W)) > 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) < 0 and np.real(gamma_be(W)) > 0:
    if W < cTe(b) and cT(b) < W < min(c0(b), vA(b)):
        return (1j * b * gamma_1_a0(W) - gamma_0_a0(W)**2 + 1) / (W**2 -b) +\
        d * (1j * b * gamma_1_ae(W) - gamma_0_ae(W)**2 + 1) / (W**2 - d*b) 
    else:
        raise Exception('incorrect range')
def deriv_dr_0_2_aea0(W):
 #   if np.real(gamma_a0(W)) > 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) < 0 and np.real(gamma_be(W)) > 0:
    if W < cTe(b) and cT(b) < W < min(c0(b), vA(b)):
        return 1j* W * (d ** 2 * (2 * ((1 + b) ** 3 * W ** 12 - b * (b + 2) * (1 + b) ** 2 * (d + 1) * W ** 10 + b ** 2 * (d ** 2 + 1) * (1 + b) * (b ** 2 + 5 * b + 5) * W ** 8 / 2 - b ** 3 * (d + 1) * (d - 1) ** 2 * (b + 2) * (1 + b) * W ** 6 - (b ** 2 * d ** 2 + (-d ** 4 / 2 + 6 * d ** 2 - 1 / 2) * b - d ** 4 / 2 + 6 * d ** 2 - 1 / 2) * b ** 4 * W ** 4 + b ** 5 * d ** 2 * (d + 1) * (b + 2) * W ** 2 - b ** 6 * d ** 2 * (d ** 2 + 1) / 2) * G * np.sqrt((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) + (W ** 2 - b * d) ** 2 * ((1 + b) * W ** 2 - b * d) ** 2 * (-2 * (1 + b) ** 2 * W ** 8 + (1 + b) * (G ** 2 + 2 * b ** 2 + 8 * b + 4) * W ** 6 + (-6 * b ** 3 - 18 * b ** 2 - 12 * b) * W ** 4 - b ** 2 * (G ** 2 - 8 * b - 12) * W ** 2 - 4 * b ** 3)) * np.sqrt((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) + (-2 * d * (1 + b) ** 2 * W ** 8 + (1 + b) * (2 * b ** 2 * d ** 2 + 8 * b * d ** 2 + G ** 2 + 4 * d ** 2) * W ** 6 - 6 * b * d ** 3 * (b + 2) * (1 + b) * W ** 4 - b ** 2 * d ** 2 * (-8 * b * d ** 2 + G ** 2 - 12 * d ** 2) * W ** 2 - 4 * b ** 3 * d ** 5) * np.sqrt((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) * (W ** 2 - b) ** 2 * ((1 + b) * W ** 2 - b) ** 2) * ((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) ** (-1 / 2) * b * ((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) ** (-1 / 2) / d ** 2 / (W ** 2 - b * d) ** 2 / (W ** 2 - b) ** 2 / ((1 + b) * W ** 2 - b) ** 2 / ((1 + b) * W ** 2 - b * d) ** 2
    else:
        raise Exception('incorrect range')
#beb0
def disp_rel_0_2_beb0(W):
#    if np.real(gamma_a0(W))  < 0 and np.real(gamma_b0(W)) > 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) < 0:
    if W < cT(b) and cTe(b) < W < min(ce(b), vAe(b)):
        return 1j*b*(gamma_0_be(W) / (W**2 +b) + d * gamma_0_b0(W)) / (W**2 -d * b)
    else:
        raise Exception('incorrect range')
def disp_rel_1_2_beb0(W):
#    if np.real(gamma_a0(W))  < 0 and np.real(gamma_b0(W)) > 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) < 0:
    if W < cT(b) and cTe(b) < W < min(ce(b), vAe(b)):
        return  (1j * b * gamma_1_b0(W) - gamma_0_b0(W)**2 + 1) / (W**2 -b) +\
        d * (1j * b * gamma_1_be(W) - gamma_0_be(W)**2 + 1) / (W**2 - d*b)
    else:
        raise Exception('incorrect range')
def deriv_dr_0_2_beb0(W):
#    if np.real(gamma_a0(W))  < 0 and np.real(gamma_b0(W)) > 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) < 0:
    if W < cT(b) and cTe(b) < W < min(ce(b), vAe(b)):
        return -1j * b * (d ** 2 * (-2 * ((1 + b) ** 3 * W ** 12 - b * (b + 2) * (1 + b) ** 2 * (d + 1) * W ** 10 + b ** 2 * (d ** 2 + 1) * (1 + b) * (b ** 2 + 5 * b + 5) * W ** 8 / 2 - b ** 3 * (d + 1) * (d - 1) ** 2 * (b + 2) * (1 + b) * W ** 6 - b ** 4 * (b ** 2 * d ** 2 + (-d ** 4 / 2 + 6 * d ** 2 - 1 / 2) * b - d ** 4 / 2 + 6 * d ** 2 - 1 / 2) * W ** 4 + b ** 5 * d ** 2 * (d + 1) * (b + 2) * W ** 2 - b ** 6 * d ** 2 * (d ** 2 + 1) / 2) * G * np.sqrt((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) + ((1 + b) * W ** 2 - b * d) ** 2 * (W ** 2 - b * d) ** 2 * (-2 * (1 + b) ** 2 * W ** 8 + (1 + b) * (G ** 2 + 2 * b ** 2 + 8 * b + 4) * W ** 6 + (-6 * b ** 3 - 18 * b ** 2 - 12 * b) * W ** 4 - b ** 2 * (G ** 2 - 8 * b - 12) * W ** 2 - 4 * b ** 3)) * np.sqrt((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) + np.sqrt((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) * ((1 + b) * W ** 2 - b) ** 2 * (W ** 2 - b) ** 2 * (-2 * d * (1 + b) ** 2 * W ** 8 + (1 + b) * (2 * b ** 2 * d ** 2 + 8 * b * d ** 2 + G ** 2 + 4 * d ** 2) * W ** 6 - 6 * b * d ** 3 * (b + 2) * (1 + b) * W ** 4 - b ** 2 * d ** 2 * (-8 * b * d ** 2 + G ** 2 - 12 * d ** 2) * W ** 2 - 4 * b ** 3 * d ** 5)) * ((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) ** (-1 / 2) * ((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) ** (-1 / 2) * W / d ** 2 / ((1 + b) * W ** 2 - b * d) ** 2 / (W ** 2 - b * d) ** 2 / ((1 + b) * W ** 2 - b) ** 2 / (W ** 2 - b) ** 2
    else:
        raise Exception('incorrect range')
# bea0
def disp_rel_0_2_bea0(W):
#    if np.real(gamma_a0(W))  > 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) < 0:
    if cT(b) < W < min(c0(b), vA(b)) and cTe(b) < W < min(ce(b), vAe(b)):
        return 1j*b*(gamma_0_be(W) / (W**2 +b) + d * gamma_0_a0(W)) / (W**2 -d * b)
    else:
        raise Exception('incorrect range')
def disp_rel_1_2_bea0(W):
#    if np.real(gamma_a0(W))  > 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) < 0:
    if cT(b) < W < min(c0(b), vA(b)) and cTe(b) < W < min(ce(b), vAe(b)):
        return (1j * b * gamma_1_a0(W) - gamma_0_a0(W)**2 + 1) / (W**2 -b) +\
        d * (1j * b * gamma_1_be(W) - gamma_0_be(W)**2 + 1) / (W**2 - d*b)
    else:
        raise Exception('incorrect range')
def deriv_dr_0_2_bea0(W):
#    if np.real(gamma_a0(W))  > 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) < 0:
    if cT(b) < W < min(c0(b), vA(b)) and cTe(b) < W < min(ce(b), vAe(b)):
        return 1j * (d ** 2 * (2 * ((1 + b) ** 3 * W ** 12 - b * (b + 2) * (1 + b) ** 2 * (d + 1) * W ** 10 + b ** 2 * (d ** 2 + 1) * (1 + b) * (b ** 2 + 5 * b + 5) * W ** 8 / 2 - b ** 3 * (d + 1) * (d - 1) ** 2 * (b + 2) * (1 + b) * W ** 6 - b ** 4 * (b ** 2 * d ** 2 + (-d ** 4 / 2 + 6 * d ** 2 - 1 / 2) * b - d ** 4 / 2 + 6 * d ** 2 - 1 / 2) * W ** 4 + b ** 5 * d ** 2 * (d + 1) * (b + 2) * W ** 2 - b ** 6 * d ** 2 * (d ** 2 + 1) / 2) * G * np.sqrt((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) + ((1 + b) * W ** 2 - b * d) ** 2 * (W ** 2 - b * d) ** 2 * (-2 * (1 + b) ** 2 * W ** 8 + (1 + b) * (G ** 2 + 2 * b ** 2 + 8 * b + 4) * W ** 6 + (-6 * b ** 3 - 18 * b ** 2 - 12 * b) * W ** 4 - b ** 2 * (G ** 2 - 8 * b - 12) * W ** 2 - 4 * b ** 3)) * np.sqrt((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) - np.sqrt((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) * ((1 + b) * W ** 2 - b) ** 2 * (W ** 2 - b) ** 2 * (-2 * d * (1 + b) ** 2 * W ** 8 + (1 + b) * (2 * b ** 2 * d ** 2 + 8 * b * d ** 2 + G ** 2 + 4 * d ** 2) * W ** 6 - 6 * b * d ** 3 * (b + 2) * (1 + b) * W ** 4 - b ** 2 * d ** 2 * (-8 * b * d ** 2 + G ** 2 - 12 * d ** 2) * W ** 2 - 4 * b ** 3 * d ** 5)) * b * ((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) ** (-1 / 2) * ((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) ** (-1 / 2) * W / d ** 2 / ((1 + b) * W ** 2 - b * d) ** 2 / (W ** 2 - b * d) ** 2 / ((1 + b) * W ** 2 - b) ** 2 / (W ** 2 - b) ** 2
    else:
        raise Exception('incorrect range')

#a0b0        
def disp_rel_0_2_a0b0(W):
#    if np.real(gamma_a0(W))  > 0 and np.real(gamma_b0(W)) > 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) > 0:
    if min(c0(b), vA(b)) < W < max(c0(b), vA(b)) and min(ce(b), vAe(b)) < W < max(ce(b), vAe(b)):
        return 1j*b*(gamma_0_a0(W)-gamma_0_b0(W)) / (W**2 -b)
    else:
        raise Exception('incorrect range') 
def disp_rel_1_2_a0b0(W):
#    if np.real(gamma_a0(W))  > 0 and np.real(gamma_b0(W)) > 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) > 0:
    if min(c0(b), vA(b)) < W < max(c0(b), vA(b)) and min(ce(b), vAe(b)) < W < max(ce(b), vAe(b)):
        return 1j * b * (gamma_0_a0(W) - gamma_0_b0(W)) / (W**2 -b)
    else:
        raise Exception('incorrect range')
def deriv_dr_0_2_a0b0(W):
#    if np.real(gamma_a0(W))  > 0 and np.real(gamma_b0(W)) > 0 and np.real(gamma_ae(W)) > 0 and np.real(gamma_be(W)) > 0:
    if min(c0(b), vA(b)) < W < max(c0(b), vA(b)) and min(ce(b), vAe(b)) < W < max(ce(b), vAe(b)):
        return 2j * ((-4 * b - 4) * W ** 6 + (G ** 2 + 4 * b ** 2 + 12 * b + 4) * W ** 4 + (-8 * b ** 2 - 8 * b) * W ** 2 + 4 * b ** 2) ** (-1 / 2) * b * (-2 * (1 + b) ** 2 * W ** 8 + (1 + b) * (G ** 2 + 2 * b ** 2 + 8 * b + 4) * W ** 6 + (-6 * b ** 3 - 18 * b ** 2 - 12 * b) * W ** 4 - b ** 2 * (G ** 2 - 8 * b - 12) * W ** 2 - 4 * b ** 3) * W / ((1 + b) * W ** 2 - b) ** 2 / (W ** 2 - b) ** 2
    else:
        raise Exception('incorrect range')
        
#aebe        
def disp_rel_0_2_aebe(W):
#    if np.real(gamma_a0(W))  < 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) < 0 and np.real(gamma_be(W)) < 0:
    if W > max(c0(b), vA(b)) and W > max(ce(b), vAe(b)):
        return 1j*b*(gamma_0_ae(W)-gamma_0_be(W)) / (W**2 -b)
    else:
        raise Exception('incorrect range')
def disp_rel_1_2_aebe(W):
#    if np.real(gamma_a0(W))  < 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) < 0 and np.real(gamma_be(W)) < 0:
    if W > max(c0(b), vA(b)) and W > max(ce(b), vAe(b)):
        return 1j * b * d *(gamma_0_ae(W) - gamma_0_be(W)) / (W**2 - d * b)
    else:
        raise Exception('incorrect range')  
def deriv_dr_0_2_aebe(W):
#    if np.real(gamma_a0(W))  < 0 and np.real(gamma_b0(W)) < 0 and np.real(gamma_ae(W)) < 0 and np.real(gamma_be(W)) < 0:
    if W > max(c0(b), vA(b)) and W > max(ce(b), vAe(b)):
        return 2j * b * ((4 * b ** 2 * d ** 4 - 8 * W ** 2 * b * (1 + b) * d ** 3 + 4 * W ** 4 * (b ** 2 + 3 * b + 1) * d ** 2 - 4 * W ** 6 * (1 + b) * d + G ** 2 * W ** 4) / d ** 4) ** (-1 / 2) * W * (-2 * d * (1 + b) ** 2 * W ** 8 + (1 + b) * (2 * b ** 2 * d ** 2 + 8 * b * d ** 2 + G ** 2 + 4 * d ** 2) * W ** 6 - 6 * b * d ** 3 * (b + 2) * (1 + b) * W ** 4 - b ** 2 * d ** 2 * (-8 * b * d ** 2 + G ** 2 - 12 * d ** 2) * W ** 2 - 4 * b ** 3 * d ** 5) / d ** 2 / ((1 + b) * W ** 2 - b * d) ** 2 / (W ** 2 - b * d) ** 2

    

# dispersion relation with 3 unknowns
def disp_rel_0_3(ga0, gb0, gc0, W):
    return (ga0 - gb0) * W**2 *gc0 / ((W**2 - b) * (b*d - W**2))

#        
def disp_rel_4(ga0, gb0, gae, gbe, W):
        return b * (gae - gbe) * (ga0 - gb0) * ( 4 *d**2 * (gb0 * ga0 + gbe * gae)* b**3\
        - d * b**2 * (d*(W**2 * (gb0*ga0 + 5*gae*gbe) + 2*ga0*gb0 + 2*gae*gbe +4) \
        + W**2 * (5*gb0*ga0 + gae*gbe)) + (d + 1) * b * W**2 * (W**2 * ga0 * gb0 \
        + d * (W**2 * gae * gbe + 2 * gb0 * ga0 + 2 * gae * gbe + 4) + W**2 * ga0 * gb0)\
        - 2 * d * W**4 * (ga0 * gb0 + gae * gbe +2)) / ((W**2 -b)**2 * (b*d - W**2)**2)


#dispersion relation for all cases
def big_disp_0(W):
    gae = gamma_0_ae(W)
    gbe = gamma_0_be(W)
    ga0 = gamma_0_a0(W)
    gb0 = gamma_0_b0(W)
    if np.real(W**2) < b / (b+1):
        if np.real(W**2) < d * b / (b+1):
            return disp_rel_0_2_aeb0(W)
        elif  d * b / (b+1) < np.real(W**2) < min(d, d*b): 
            return disp_rel_0_2_beb0(W)
        elif  max(d, d*b) < np.real(W**2):
            return disp_rel_0_3(gbe, gae, gb0, W)
        else:
            return np.Nan
    elif  b / (b+1) < np.real(W**2) < min(1, b):
        if np.real(W**2) < d * b / (b+1):
            return disp_rel_0_2_aea0(W)
        elif  d * b / (b+1) < np.real(W**2) < min(d, d*b): 
            return disp_rel_0_2_bea0(W)
        elif  max(d, d*b) < np.real(W**2):
            return disp_rel_0_3(gbe, gae, ga0, W)
        else:
            return np.Nan
    elif min(1, b) < np.real(W**2) < max(1,b):
        if np.real(W**2) < d * b / (b+1):
            return disp_rel_0_3(ga0, gb0, gae, W)
        elif  d * b / (b+1) <np.real(W**2) < min(d, d*b): 
            return disp_rel_0_3(ga0, gb0, gbe, W)
        elif min(d, d*b)< np.real(W**2) <max(d, d*b):
            return disp_rel_0_2_a0b0(W)
        elif  max(d, d*b) < np.real(W**2):
            return disp_rel_4(ga0, gb0, gae, gbe, W)
        else:
            return np.NaN
    elif  np.real(W**2) > max(1,b) and np.real(W**2) > max(d, d*b):
        return disp_rel_0_2_aebe(W)
    else:
        return np.NaN
        
def big_disp_1(W):
    if np.real(W**2) < b / (b+1):
        if np.real(W**2) < d * b / (b+1):
            return disp_rel_1_2_aeb0(W)
        elif  d * b / (b+1) < np.real(W**2) < min(d, d*b): 
            return disp_rel_1_2_beb0(W)
        elif  max(d, d*b) < np.real(W**2):
            return np.NaN #disp_rel_1_3(gbe, gae, gb0, W)
        else:
            return np.NaN
    elif  b / (b+1) < np.real(W**2) < min(1, b):
        if np.real(W**2) < d * b / (b+1):
            return disp_rel_1_2_aea0(W)
        elif  d * b / (b+1) < np.real(W**2) < min(d, d*b): 
            return disp_rel_1_2_bea0(W)
        elif  max(d, d*b) < np.real(W**2):
            return np.NAN #disp_rel_1_3(gbe, gae, ga0, W)
        else:
            return np.NaN
    elif min(1, b) < np.real(W**2) < max(1,b):
        if np.real(W**2) < d * b / (b+1):
            return np.NAN #disp_rel_1_3(ga0, gb0, gae, W)
        elif  d * b / (b+1) < np.real(W**2) < min(d, d*b): 
            return np.NAN #disp_rel_1_3(ga0, gb0, gbe, W)
        elif min(d, d*b)< np.real(W**2) <max(d, d*b):
            return disp_rel_1_2_a0b0(W)
        elif  max(d, d*b) < np.real(W**2):
            return np.NAN #disp_rel_1_4(ga0, gb0, gae, gbe, W)
        else:
            return np.NaN
    elif  np.real(W**2) > max(1,b) and np.real(W**2)> max(d, d*b):
        return disp_rel_1_2_aebe(W)
    else:
        return np.NaN
        
def big_dr_0_deriv(W):
    if np.real(W**2) < b / (b+1):
        if np.real(W**2) < d * b / (b+1):
            return deriv_dr_0_2_aeb0(W)
        elif  d * b / (b+1) < np.real(W**2) < min(d, d*b): 
            return deriv_dr_0_2_beb0(W)
        elif  max(d, d*b) < np.real(W**2):
            return np.NaN #disp_rel_1_3(gbe, gae, gb0, W)
        else:
            return np.NaN
    elif  b / (b+1) < np.real(W**2) < min(1, b):
        if np.real(W**2) < d * b / (b+1):
            return deriv_dr_0_2_aea0(W)
        elif  d * b / (b+1) < np.real(W**2) < min(d, d*b): 
            return deriv_dr_0_2_bea0(W)
        elif  max(d, d*b) < np.real(W**2):
            return np.NAN #disp_rel_1_3(gbe, gae, ga0, W)
        else:
            return np.NaN
    elif min(1, b) < np.real(W**2) < max(1,b):
        if np.real(W**2) < d * b / (b+1):
            return np.NAN #disp_rel_1_3(ga0, gb0, gae, W)
        elif  d * b / (b+1) < np.real(W**2) < min(d, d*b): 
            return np.NAN #disp_rel_1_3(ga0, gb0, gbe, W)
        elif min(d, d*b)< np.real(W**2) <max(d, d*b):
            return deriv_dr_0_2_a0b0(W)
        elif  max(d, d*b) < np.real(W**2):
            return np.NAN #disp_rel_1_4(ga0, gb0, gae, gbe, W)
        else:
            return np.NaN
    elif  np.real(W**2) > max(1,b) and np.real(W**2) > max(d, d*b):
        return deriv_dr_0_2_aebe(W)
    else:
        return np.NaN

def perturbed(W):
    W1 = -big_disp_1(W) / big_dr_0_deriv(W)
    return W + theta * W1

b_num=100
b_start = -2
b_end = 2
b_range=np.logspace(b_start,b_end,b_num)
b=np.logspace(b_start,b_end,b_num)
#d=np.linspace(b_start, b_end,b_num)

snum=200
s_start =0
s_end=4
s_range = np.linspace(s_start,s_end,snum)

n_0roots=50
megaroots=np.nan*np.ones(shape=(n_0roots,b_num),dtype=np.complex64)


#Check it's not finding the same root multiple times
def check(root_temp, root):
    n = True
    if len(root_temp)==0:
        root_temp.append(root)
    else:
        for x in root_temp:
            if np.abs(np.real(x) - np.real(root)) < 10**(-6) and \
               np.abs(np.imag(x) - np.imag(root)) < 10**(-6):
                   n = False
            else:
                pass
        if n == True:
            root_temp.append(root)

# find list of roots
root_list = []
for b in b_range:
    root_temp = []
    for s in s_range:
        try:
            root = sp.newton(big_disp_0, s + 0.01*1j, tol = 10**(-10))
            print root
            check(root_temp, root)
        except Exception as error:
            print error
            pass
    root_list.append(root_temp)


np.save('RTIroots0',root_list)
root_list =np.load('RTIroots0.npy')

    
#plot all sols
plt.subplots()
plt.tight_layout
ax = plt.subplot(211)
#plt.title('Solutions, beta=0.2') 
ax2 = plt.subplot(212)
ax.set_ylabel(r'W0',size=20)
#ax.set_xlabel(r'beta',size=17)
ax2.set_ylabel(r'W',size=20)
ax2.set_xlabel(r'beta',size=17)
ax.set_ylim(0, 2)
ax2.set_ylim(0, 2)



#plot perturbed sols
for i,root in enumerate(root_list):
    for j in root:
        if np.real(j)>0:        
            W = perturbed(j)
            ax2.plot(b_range[i],  np.real(W), color='b', marker='.')
            ax2.plot(b_range[i],  np.imag(W), color='r', marker='.')
ax2.set_xscale('log')

#plot zero_order sols
for i,root in enumerate(root_list):
    for j in root:
        if np.real(j)>0:        
            W = perturbed(j)
            if np.isnan(W) == False:
                ax.plot(b_range[i],  np.real(j), color='b', marker='.')
                ax.plot(b_range[i],  np.imag(j), color='r', marker='.')
            else:
                ax.plot(b_range[i],  np.real(j), color='g', marker='.')
                ax.plot(b_range[i],  np.imag(j), color='m', marker='.')
ax.set_xscale('log')

# plot critical speeds

ax.plot(b_range, vA(b_range)*np.ones(len(b_range)) , 'g')
ax.plot(b_range, vAe(b_range)*np.ones(len(b_range)), 'r')
ax.plot(b_range, c0(b_range), 'g')
ax.plot(b_range, ce(b_range), 'r')
ax.plot(b_range, cT(b_range), 'g--')
ax.plot(b_range, cTe(b_range), 'r--')

ax2.plot(b_range, vA(b_range)*np.ones(len(b_range)) , 'g')
ax2.plot(b_range, vAe(b_range)*np.ones(len(b_range)), 'r')
ax2.plot(b_range, c0(b_range), 'g')
ax2.plot(b_range, ce(b_range), 'r')
ax2.plot(b_range, cT(b_range), 'g--')
ax2.plot(b_range, cTe(b_range), 'r--')


plt.show