import solvers as sol
from AS1_class import Asym_slab

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle

def save():
    with open('pickles/flow_c0={}_R1={}_R2={}_K={}_M_A={}.p'.format(
    				slab.c0, slab.R1, slab.R2, slab.K, slab.M_A), 'wb') as f:
    	    pickle.dump(root_array, f)

slab = Asym_slab(c0=0.6, R1=0.1, R2=0.1, K=10.0, M_A=None)

x_range = np.linspace(1, 3, 201)
y_range = np.linspace(0, slab.c1, 201)

load = False

if load:
    root_array = pickle.load(open('pickles/flow_c0={}_R1={}_R2={}_K={}_M_A={}.p'.format(
							slab.c0, slab.R1, slab.R2, slab.K, slab.M_A), 'rb'))
else:
    root_array = sol.point_finder_sp(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))

font = {'size': 15}
matplotlib.rc('font', **font)

plt.figure(num=None, figsize=(15, 11), dpi=80, facecolor='w', edgecolor='k')
ax = plt.subplot()

ax.plot(root_array[:,0], np.real(root_array[:,1]), '.', color = 'b')
ax.plot(root_array[:,0], np.imag(root_array[:,1]), '.', color = 'r')

ax.set_xlim(x_range[0], x_range[-1])
ax.set_ylim(y_range[0], y_range[-1])

ax.set_ylabel(r'$\bar c_{ph}$', fontsize = 30, rotation = 0, labelpad=20)
ax.set_xlabel(r'$M_A$', fontsize = 30)

ax.plot([0, x_range[-1]],[slab.cT, slab.cT + x_range[-1]], color='0.5')
ax.plot([0, x_range[-1]],[-slab.cT, -slab.cT + x_range[-1]], color='0.5')
ax.plot([0, x_range[-1]],[slab.c0, slab.c0 + x_range[-1]], color='0.5')
ax.plot([0, x_range[-1]],[-slab.c0, -slab.c0 + x_range[-1]], color='0.5')
ax.plot([0, x_range[-1]],[slab.c2, slab.c2], color='0.5')
ax.plot([0, x_range[-1]],[-slab.c2, -slab.c2], color='0.5')
ax.plot([0, x_range[-1]],[0, 0], color='0.5', linestyle='--')

ax.annotate(r'$\bar c_1, \bar c_2$', xy=(x_range[-1]+0.01, slab.c2-0.02), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$-\bar c_1, -\bar c_2$', xy=(x_range[-1]+0.01, -slab.c2-0.02), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$\bar c_T + M_A$', xy=(-0.13, slab.cT - 0.04), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$-\bar c_T + M_A$', xy=(-0.16, -slab.cT - 0.02), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$\bar c_0 + M_A$', xy=(-0.13, slab.c0 - 0.02), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$-\bar c_0 + M_A$', xy=(-0.16, -slab.c0 - 0.04), xycoords='data', annotation_clip=False, fontsize=20)