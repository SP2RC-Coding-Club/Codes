import solvers as sol
from AS1_class import Asym_slab

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle

def save():
    with open('pickles/wavenumber_c0={}_R1={}_R2={}_K={}_M_A={}.p'.format(
    				slab.c0, slab.R1, slab.R2, slab.K, slab.M_A), 'wb') as f:
    	    pickle.dump(root_array, f)

slab = Asym_slab(c0=0.6, R1=5/3, R2=2.5, K=None, M_A=0.4)

x_range = np.linspace(0, 4, 101)
y_range = np.linspace(-max(slab.c1, slab.c2, 1), max(slab.c1, slab.c2, 1), 101)

load = False

if load:
    root_array = pickle.load(open('pickles/wavenumber_c0={}_R1={}_R2={}_K={}_M_A={}.p'.format(
							slab.c0, slab.R1, slab.R2, slab.K, slab.M_A), 'rb'))
else:
    root_array = sol.point_finder_sp(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))

font = {'size': 15}
matplotlib.rc('font', **font)

plt.figure(num=None, figsize=(8, 13), dpi=80, facecolor='w', edgecolor='k')
ax = plt.subplot()

#ax.plot(root_array[:,0], np.real(root_array[:,1]), '.', color = 'b')
#ax.plot(root_array[:,0], np.imag(root_array[:,1]), '.', color = 'r')

ax.set_xlim(x_range[0], x_range[-1])
ax.set_ylim(y_range[0], y_range[-1])

ax.set_ylabel(r'$\bar c_{ph}$', fontsize = 30, rotation=0)
ax.set_xlabel(r'$k x_0$', fontsize = 30)

line1 = sol.line_trace_sp(slab.disp_rel, 1, 0.635, 0.01, 0.37, 4, args=(slab.K, slab.M_A))
line2 = sol.line_trace_sp(slab.disp_rel, 1, 0.495, 0.01, 0, 4, args=(slab.K, slab.M_A))
line3 = sol.line_trace_sp(slab.disp_rel, 1, -0.383, 0.01, 0, 4, args=(slab.K, slab.M_A))
line4 = sol.line_trace_sp(slab.disp_rel, 1, -0.64, 0.01, 0.55, 4, args=(slab.K, slab.M_A))
line5 = sol.line_trace_sp(slab.disp_rel, 2.32, -0.159, 0.005, 0.06, 4, args=(slab.K, slab.M_A))
line6 = sol.line_trace_sp(slab.disp_rel, 2.24, -0.12, 0.005, 0.17, 4, args=(slab.K, slab.M_A))

ax.plot(line1[:,0], line1[:,1], color = 'b', linewidth = 2, linestyle='-', label='Quasi-Sausage')
ax.plot(line2[:,0], line2[:,1], color = 'g', linewidth = 2, linestyle='-', label='Quasi-Kink')
ax.plot(line2[:,0], np.imag(line2[:,1]), color = 'r', linewidth = 2, linestyle='-', label='Imaginary')
ax.plot(line3[:,0], line3[:,1], color = 'g', linewidth = 2, linestyle='-')
ax.plot(line4[:,0], line4[:,1], color = 'b', linewidth = 2, linestyle='-')
ax.plot(line5[:,0], line5[:,1], color = 'g', linewidth = 2, linestyle='-')
ax.plot(line6[:,0], line6[:,1], color = 'b', linewidth = 2, linestyle='-')

ax.plot([x_range[0], x_range[-1]], [slab.c0 + slab.M_A, slab.c0 + slab.M_A],
        [x_range[0], x_range[-1]], [-slab.c0 + slab.M_A, -slab.c0 + slab.M_A],
        [x_range[0], x_range[-1]], [slab.cT + slab.M_A, slab.cT + slab.M_A],
        [x_range[0], x_range[-1]], [-slab.cT + slab.M_A, -slab.cT + slab.M_A],
        [x_range[0], x_range[-1]], [1 + slab.M_A, 1 + slab.M_A],
        [x_range[0], x_range[-1]], [-1 + slab.M_A, -1 + slab.M_A],
        color = '0.5', linestyle=':', linewidth=2)

if slab.cT + slab.M_A < y_range[-1]:
    ax.annotate(r'$\bar c_T + M_A$', xy=(x_range[-1] + 0.03, slab.cT + slab.M_A - 0.01),
                xycoords='data', annotation_clip=False, fontsize=20)
if slab.c0 + slab.M_A < y_range[-1]:
    ax.annotate(r'$\bar c_0 + M_A$', xy=(x_range[-1] + 0.03, slab.c0 + slab.M_A - 0.01),
                xycoords='data', annotation_clip=False, fontsize=20)
if 1 + slab.M_A < y_range[-1]:
    ax.annotate(r'$1+M_A$', xy=(x_range[-1] + 0.03, 1+slab.M_A - 0.01),
			xycoords='data', annotation_clip=False, fontsize=20)

ax.annotate(r'$-\bar c_T + M_A$', xy=(x_range[-1] + 0.03, -slab.cT + slab.M_A - 0.01),
			xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$-\bar c_0 + M_A$', xy=(x_range[-1] + 0.03, -slab.c0 + slab.M_A - 0.01),
			xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$-1+M_A$', xy=(x_range[-1] + 0.03, -1+slab.M_A - 0.01),
			xycoords='data', annotation_clip=False, fontsize=20)

if slab.c1 == slab.c2:
    ax.plot([x_range[0], x_range[-1]], [slab.c1, slab.c1],
             [x_range[0], x_range[-1]], [-slab.c1, -slab.c1],
             color = '0.5', linestyle=':', linewidth=2)
    ax.annotate(r'$\bar c_1, \bar c_2$', xy=(x_range[-1] + 0.03, slab.c1 - 0.01),
                 xycoords='data', annotation_clip=False, fontsize=20)
    ax.annotate(r'$-\bar c_1, -\bar c_2$', xy=(x_range[-1] + 0.03, -slab.c1 - 0.01),
                 xycoords='data', annotation_clip=False, fontsize=20)
else:
    ax.plot([x_range[0], x_range[-1]], [slab.c1, slab.c1],
            [x_range[0], x_range[-1]], [-slab.c1, -slab.c1],
            [x_range[0], x_range[-1]], [slab.c2, slab.c2],
            [x_range[0], x_range[-1]], [-slab.c2, -slab.c2],
            color = '0.5', linestyle=':', linewidth=2)
    ax.annotate(r'$\bar c_1$', xy=(x_range[-1] + 0.03, slab.c1 - 0.01),
                xycoords='data', annotation_clip=False, fontsize=20)
    ax.annotate(r'$-\bar c_1$', xy=(x_range[-1] + 0.03, -slab.c1 - 0.01),
                xycoords='data', annotation_clip=False, fontsize=20)
    ax.annotate(r'$\bar c_2$', xy=(x_range[-1] + 0.03, slab.c2 - 0.01),
                xycoords='data', annotation_clip=False, fontsize=20)
    ax.annotate(r'$-\bar c_2$', xy=(x_range[-1] + 0.03, -slab.c2 - 0.01),
                xycoords='data', annotation_clip=False, fontsize=20)

ax.fill_between(x_range, -slab.c0 + slab.M_A, -slab.cT + slab.M_A, color='0.9')
ax.fill_between(x_range, np.maximum(-slab.c1, -slab.c2), -1 + slab.M_A, color='0.9')
ax.fill_between(x_range, np.maximum(-slab.c1, -slab.c2), y_range[0], color='1', hatch='/', edgecolor='0')
ax.fill_between(x_range, np.minimum(slab.c1, slab.c2), y_range[-1], color='1', hatch='/', edgecolor='0')

plt.legend(loc=1)

plt.savefig('images_traced/dense2.png', bbox_inches='tight')