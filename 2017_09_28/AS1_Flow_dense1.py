import solvers as sol
from AS_funcs import Asym_slab

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle

def save():
    with open('pickles/flow_c0={}_R1={}_R2={}_K={}_M_A={}.p'.format(
    				slab.c0, slab.R1, slab.R2, slab.K, slab.M_A), 'wb') as f:
    	    pickle.dump(root_array, f)

slab = Asym_slab(c0=0.6, R1=1.4, R2=1.4, K=0.5, M_A=None)

x_range = np.linspace(0, 1.5, 201)
y_range = np.linspace(-1, 1, 201)

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

#ax.plot(root_array[:,0], np.real(root_array[:,1]), '.', color = 'b')
#ax.plot(root_array[:,0], np.imag(root_array[:,1]), '.', color = 'r')

ax.set_xlim(x_range[0], x_range[-1])
ax.set_ylim(y_range[0], y_range[-1])

ax.set_ylabel(r'$\bar c_{ph}$', fontsize = 30, rotation = 0, labelpad=20)
ax.set_xlabel(r'$M_A$', fontsize = 30)

ax.plot([0, 0.4088],[slab.cT, slab.cT + 0.4088], color='0.5')
ax.plot([0, 1.4377],[-slab.cT, -slab.cT + 1.4377], color='0.5')
ax.plot([0, 0.3232],[slab.c0, slab.c0 + 0.3232], color='0.5')
ax.plot([0, x_range[-1]],[-slab.c0, -slab.c0 + x_range[-1]], color='0.5')
ax.plot([0, x_range[-1]],[slab.c2, slab.c2], color='0.5')
ax.plot([0, x_range[-1]],[-slab.c2, -slab.c2], color='0.5')
ax.plot([0, x_range[-1]],[0, 0], color='0.5', linestyle='--')

ax.fill_between(x_range, slab.c2, 1, color='1', hatch="/", edgecolor='0.5')
ax.fill_between(x_range, -slab.c2, -1, color='1', hatch="/", edgecolor='0.5')

x=np.linspace(0, 0.4088, 101)
ax.fill_between(x, slab.cT + x, np.minimum(slab.c0 + x, slab.c2), color='0.9')
x=np.linspace(0, x_range[-1], 101)
ax.fill_between(x, -slab.c0 + x, np.minimum(-slab.cT + x, slab.c2), color='0.9')

ax.set_yticks([-1.0, -0.75, -0.25, 0, 0.25, 0.75, 1.0])
ax.set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])

ax.annotate(r'$\bar c_1, \bar c_2$', xy=(x_range[-1]+0.01, slab.c2-0.02), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$-\bar c_1, -\bar c_2$', xy=(x_range[-1]+0.01, -slab.c2-0.02), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$\bar c_T + M_A$', xy=(-0.13, slab.cT - 0.04), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$-\bar c_T + M_A$', xy=(-0.16, -slab.cT - 0.02), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$\bar c_0 + M_A$', xy=(-0.13, slab.c0 - 0.02), xycoords='data', annotation_clip=False, fontsize=20)
ax.annotate(r'$-\bar c_0 + M_A$', xy=(-0.16, -slab.c0 - 0.04), xycoords='data', annotation_clip=False, fontsize=20)

line1 = sol.line_trace_sp(slab.disp_rel, 0.0975, 0.8839, 0.005, 0, 0.3232, args=(slab.K, slab.M_A))
line2 = sol.line_trace_sp(slab.disp_rel, 0.6, 0.8197, 0.01, -0.01, 1.5, args=(slab.K, slab.M_A))
line3 = sol.line_trace_sp(slab.disp_rel, 0.6, 0.5205, 0.001, 0, 1.5, args=(slab.K, slab.M_A))
line4 = sol.line_trace_sp(slab.disp_rel, 0.6, -0.8537, 0.01, -0.01, 1.51, args=(slab.K, slab.M_A))
line5 = sol.line_trace_sp(slab.disp_rel, 0.6, -0.2779, 0.001, 0, 1.156, args=(slab.K, slab.M_A))
line6 = sol.line_trace_sp(slab.disp_rel, 0.06, 0.5841, 0.001, 0, 0.403, args=(slab.K, slab.M_A))
line7 = sol.line_trace_mp(slab.disp_rel_mp, 0.06, -0.435, 0.01, -0.01, 1, args=(slab.K, slab.M_A))
line8 = sol.line_trace_sp(slab.disp_rel, 1.1025, 0.614, 0.001, 0.99, 1.187, args=(slab.K, slab.M_A))
line9 = sol.line_trace_sp(slab.disp_rel, 0.0075, -0.4957, 0.0001, 0, 0.2, args=(slab.K, slab.M_A))
line10 = sol.line_trace_sp(slab.disp_rel, 0.9225, 0.4216, 0.0001, 0.8, 0.933, args=(slab.K, slab.M_A))
line11 = sol.line_trace_sp(slab.disp_rel, 1.11, 0.586, 0.0001, 1.027, 1.442, args=(slab.K, slab.M_A))

ax.plot(line1[:,0], line1[:,1], color = 'b', linewidth = 2, linestyle='-', label='Stable')
ax.plot(line2[:,0][(np.imag(line2[:,1]) < 2e-2)], line2[:,1][(np.imag(line2[:,1]) < 2e-2)],
        color = 'b', linewidth = 2, linestyle='-')
ax.plot(line2[:,0][(np.imag(line2[:,1]) > 1e-2)], line2[:,1][(np.imag(line2[:,1]) > 1e-2)],
        color = 'b', linewidth = 2, linestyle='--', label='Unstable')
ax.plot(line2[118:,0], np.imag(line2[118:,1]), color = 'r', linewidth = 2, linestyle='-', label='Imaginary')
ax.plot(line2[118:,0], -np.imag(line2[118:,1]), color = 'r', linewidth = 2, linestyle='-')
ax.plot(line3[:932,0], np.real(line3[:932,1]), color = 'b', linewidth = 2, linestyle='-')
ax.plot(line3[931:1028,0], np.real(line3[931:1028,1]), color = 'b', linewidth = 2, linestyle='--')
ax.plot(line3[931:1028,0], np.imag(line3[931:1028,1]), color = 'r', linewidth = 2, linestyle='-')
ax.plot(line3[931:1028,0], -np.imag(line3[931:1028,1]), color = 'r', linewidth = 2, linestyle='-')
ax.plot(line3[1027:1156,0], np.real(line3[1027:1156,1]), color = 'b', linewidth = 2, linestyle='-')
ax.plot(line3[1154:,0], np.real(line3[1154:,1]), color = 'b', linewidth = 2, linestyle='--')
ax.plot(line3[1154:,0], np.imag(line3[1154:,1]), color = 'r', linewidth = 2, linestyle='-')
ax.plot(line3[1154:,0], -np.imag(line3[1154:,1]), color = 'r', linewidth = 2, linestyle='-')
ax.plot(line4[:,0], line4[:,1], color = 'b', linewidth = 2, linestyle='-')
ax.plot(line5[:,0], line5[:,1], color = 'b', linewidth = 2, linestyle='-')
ax.plot(line6[:,0], line6[:,1], color = 'b', linewidth = 2, linestyle='-')
ax.plot(line7[:,0], line7[:,1], color = 'b', linewidth = 2, linestyle='-')
ax.plot(line8[:,0], line8[:,1], color = 'b', linewidth = 2, linestyle='-')
ax.plot(line9[:,0], line9[:,1], color = 'b', linewidth = 2, linestyle='-')
ax.plot(line10[:,0], line10[:,1], color = 'b', linewidth = 2, linestyle='-')
ax.plot(line11[:,0], line11[:,1], color = 'b', linewidth = 2, linestyle='-')

plt.legend(loc=4, fancybox=True, framealpha=0.9)

plt.subplots_adjust(bottom=0.1, top=0.95)

plt.savefig('images_traced/flow1.png')