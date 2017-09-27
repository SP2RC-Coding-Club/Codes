#import toolbox as tool
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy as sc
#import scipy.optimize as sp
#import mpmath as mp
from scipy.interpolate import interp1d
import pickle
from matplotlib import gridspec

#field.text = data.decode("utf8")

def m1(W):
    return sc.sqrt(1 - W**2 / c1**2)
    
def m2(W):
    return sc.sqrt(1 - W**2 / c2**2)

def m0(W, M_A):
    return sc.sqrt((1 - (W - M_A)**2) * (c0**2 - (W - M_A)**2) /
                      ((1 + c0**2) * (cT**2 - (W - M_A)**2)))

def disp_rel(W, K, M_A):
    return m0(W, M_A)**2 * W**4 + R1 * R2 * m1(W) * m2(W) * (1 - (W - M_A)**2)**2 - \
           0.5 * W**2 * m0(W, M_A) * (1 - (W - M_A)**2) * (R1 * m1(W) + R2 * m2(W)) * \
           (sc.tanh(m0(W, M_A) * K) + sc.tanh(m0(W, M_A) * K)**(-1))

c0 = 0.6
R1 = 1.5
R2 = 0.5
c1 = sc.sqrt(R1 * (c0**2 + 5/6))
c2 = sc.sqrt(R2 * (c0**2 + 5/6))
cT = sc.sqrt(c0**2 / (1 + c0**2))

K_start = 0.5   # starting value of K
K_left  = 0.02  # limit on the left
K_right = 3     # limit on the right
K_no_right = 50 # number of points to find right of K_start
K_no_left = 20  # number of points to find left of K_start
K_range_right = np.linspace(K_start, K_right, K_no_right)
K_range_left = np.linspace(K_left, K_start, K_no_left, endpoint=False)[::-1]

threshold = []
roots = []

x_start = []
y_start = []

"""
First loop:
The code increases (goes to the right of) K_start and attempts to find the
value of M_A for which the dispersion relation has complex solutions.

Remember to change the values of x_range and y_range at the start of the first
loop.

Second loop:
Now K is decreased (again, starting from K_start), and the same thing happens.
"""

#for K in K_range_right:
#    if len(threshold) == 0:
#        x_range = np.linspace(0.89, 0.93, 101)
#        y_range = np.linspace(0.4, 0.45, 101)
#    elif len(threshold) == 1:
#        x_range = np.linspace(x_start[-1], x_start[-1] + 0.1, 51)
#        y_range = np.linspace(y_start[-1], y_start[-1] + 0.1, 51)
#    elif len(threshold) > 1:
#        x_range = np.linspace(x_start[-1], x_start[-1] + (threshold[1] - threshold[0] + (threshold[1] - threshold[0])/4), 51)
#        y_range = np.linspace(y_start[-1], y_start[-1] + (roots[1] - roots[0] + (roots[1] - roots[0])/4), 51)
#    for x in x_range:
#        for y in y_range:
#            try:
#                root = sp.newton(disp_rel, y, args = (K, x), tol = 1e-20)
#                if np.imag(root) > 1e-6:
#                    threshold.append(x)
#                    roots.append(root)
#                    break
#            except RuntimeError:
#                pass
#        else: #this is executed if the for loop exited normally
#            continue
#        break #this is executed if the for loop did not exit normally
#    x_start.append(threshold[-1])
#    y_start.append(roots[-1])
#
#print('Done.')
#    
#for K in K_range_left:
#    x_range = np.linspace(x_start[0], x_start[0] - (threshold[1] - threshold[0])*2, 51)[::-1]#(threshold[1] - threshold[0] + (threshold[1] - threshold[0])/2), 51)[::-1]
#    y_range = np.linspace(y_start[0], y_start[0] - (roots[1] - roots[0])*2, 51)[::-1]
#    for x in x_range:
#        for y in y_range:
#            try:
#                root = sp.newton(disp_rel, y, args = (K, x), tol = 1e-20)
#                if np.imag(root) > 1e-6:
#                    threshold.insert(0, x)
#                    roots.insert(0, root)
#                    break
#            except RuntimeError:
#                pass
#        else: #this is executed if the for loop exited normally
#            continue
#        break #this is executed if the for loop did not exit normally
#    x_start.insert(0, threshold[0])
#    y_start.insert(0, roots[0])
#
#print('Done.')
#    
#K_range = list(K_range_left[::-1]) + list(K_range_right)
#
#pickleness = []
#pickleness.append(K_range)
#pickleness.append(threshold)
#
#pickle.dump(pickleness, open("threshold_R1={}_R2={}.p".format(R1,R2), "wb"))

plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 3)

R_list = [(0.5, 0.5), (1.0, 1.0), (1.5, 1.5)]
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

font = {'size': 15}
matplotlib.rc('font', **font)

ax1 = plt.subplot(gs[:, 0])

ax1.set_ylabel(r'$M_A$', fontsize=30, rotation=0, labelpad=20)

for i, (R1, R2) in enumerate(R_list):
    K_range = pickle.load(open('pickles/threshold_R1={}_R2={}.p'.format(R1,R2), "rb"), encoding='latin1')[0]
    threshold = pickle.load(open("pickles/threshold_R1={}_R2={}.p".format(R1,R2), "rb"), encoding='latin1')[1]
    f = interp1d(K_range, threshold, kind='cubic')
    ax1.plot(K_range, threshold, '.', color = colors[i])
    if R2 == 1.5 and R1 == R2:
        ax1.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 2/3', color = colors[i])
    elif R2 == 1.0 and R1 == R2:
        ax1.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 1', color = colors[i])
    elif R2 == 0.5 and R1 == R2:
        ax1.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 2', color = colors[i])
    elif R2 == 1.0:
        ax1.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0$ = 2/3, $\rho_2/\rho_0$ = 1'.format(R1, R2), color = colors[i])
    elif R2 == 0.5:
        ax1.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0$ = 2/3, $\rho_2/\rho_0$ = 2'.format(R1, R2), color = colors[i])

plt.legend(loc = 'lower right')
ax1.set_xlim(0.02, 2)
ax1.set_ylim(0.5, 1.4)

R_list = [(0.5, 0.5), (1.0, 1.0), (1.5, 1.0)]
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

ax2 = plt.subplot(gs[:, 1])

for i, (R1, R2) in enumerate(R_list):
    K_range = pickle.load(open('pickles/threshold_R1={}_R2={}.p'.format(R1,R2), "rb"), encoding='latin1')[0]
    threshold = pickle.load(open("pickles/threshold_R1={}_R2={}.p".format(R1,R2), "rb"), encoding='latin1')[1]
    f = interp1d(K_range, threshold, kind='cubic')
    ax2.plot(K_range, threshold, '.', color = colors[i])
    if R2 == 1.5 and R1 == R2:
        ax2.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 2/3', color = colors[i])
    elif R2 == 1.0 and R1 == R2:
        ax2.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 1', color = colors[i])
    elif R2 == 0.5 and R1 == R2:
        ax2.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 2', color = colors[i])
    elif R2 == 1.0:
        ax2.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0$ = 2/3, $\rho_2/\rho_0$ = 1'.format(R1, R2), color = colors[i])
    elif R2 == 0.5:
        ax2.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0$ = 2/3, $\rho_2/\rho_0$ = 2'.format(R1, R2), color = colors[i])

ax2.set_xlabel(r'$k x_0$', fontsize=30)

plt.legend(loc = 'lower right')
ax2.set_xlim(0.02, 2)
ax2.set_ylim(0.5, 1.4)

R_list = [(0.5, 0.5), (1.0, 1.0), (1.5, 0.5)]
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

ax3 = plt.subplot(gs[:, 2])

for i, (R1, R2) in enumerate(R_list):
    K_range = pickle.load(open('pickles/threshold_R1={}_R2={}.p'.format(R1,R2), "rb"), encoding='latin1')[0]
    threshold = pickle.load(open("pickles/threshold_R1={}_R2={}.p".format(R1,R2), "rb"), encoding='latin1')[1]
    f = interp1d(K_range, threshold, kind='cubic')
    ax3.plot(K_range, threshold, '.', color = colors[i])
    if R2 == 1.5 and R1 == R2:
        ax3.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 2/3', color = colors[i])
    elif R2 == 1.0 and R1 == R2:
        ax3.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 1', color = colors[i])
    elif R2 == 0.5 and R1 == R2:
        ax3.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0 = \rho_2/\rho_0$ = 2', color = colors[i])
    elif R2 == 1.0:
        ax3.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0$ = 2/3, $\rho_2/\rho_0$ = 1'.format(R1, R2), color = colors[i])
    elif R2 == 0.5:
        ax3.plot(K_range, f(K_range), label=r'$\rho_1/\rho_0$ = 2/3, $\rho_2/\rho_0$ = 2'.format(R1, R2), color = colors[i])

plt.legend(loc = 'lower right')
ax3.set_xlim(0.02, 2)
ax3.set_ylim(0.5, 1.4)

plt.savefig('images_traced/threshold.png', bbox_inches='tight')