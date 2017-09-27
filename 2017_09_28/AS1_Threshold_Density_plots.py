import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d
from matplotlib import gridspec

def load_sym(K):
    return pickle.load(open('pickles/threshold_density_sym_c0={}_K={}_M_A={}.p'.format(0.6, K, 'None'), 'rb'))

def load_asym(K):
    return pickle.load(open('pickles/threshold_density_c0={}_R2={}_K={}_M_A={}.p'.format(0.6, 2, K, 'None'), 'rb'))

gs = gridspec.GridSpec(2, 1)

plt.figure(num=None, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')

ax1 = plt.subplot(gs[0,:])

threshold1 = load_sym(0.1)
threshold2 = load_sym(1.0)
threshold3 = load_sym(10.0)

interp1 = interp1d(threshold1[:,0], threshold1[:,1], kind='cubic')
interp2 = interp1d(threshold2[:,0], threshold2[:,1], kind='cubic')
interp3 = interp1d(threshold3[:,0], threshold3[:,1], kind='cubic')

ax1.plot(threshold1[:,0], threshold1[:,1], '.', color='r')
ax1.plot(threshold2[:,0], threshold2[:,1], '.', color='g')
ax1.plot(threshold3[:,0], threshold3[:,1], '.', color='b')

ax1.plot(threshold1[:,0], interp1(threshold1[:,0]), color='r', label=r'$k x_0 = 0.1$')
ax1.plot(threshold2[:,0], interp2(threshold2[:,0]), color='g', label=r'$k x_0 = 1.0$')
ax1.plot(threshold3[:,0], interp3(threshold3[:,0]), color='b', label=r'$k x_0 = 10.0$')

ax1.legend(loc=1, fontsize=20)

ax1.set_xlim([0, 10])
ax1.set_ylim([0.55, 0.97])

ax1.set_xlabel(r'$\rho_e / \rho_0$', fontsize=20)
ax1.set_ylabel(r'$M_A$', fontsize=30, rotation=0, labelpad=20, y=-0.15)

#ax1.set_xticks([0, 2, 4, 6, 8, 10])
#ax1.set_xticklabels([0, 2, 4, 6, 8, 10], fontsize=15)
#ax1.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
#ax1.set_yticklabels([0.6, 0.8, 1.0, 1.2, 1.4], fontsize=15)

ax2 = plt.subplot(gs[1,:])

ax2.plot([2,2],[0,1], '--', color='0.8')

threshold4 = load_asym(0.1)
threshold5 = load_asym(1.0)
threshold6 = load_asym(10.0)

interp4 = interp1d(threshold4[:,0], threshold4[:,1], kind='cubic')
interp5 = interp1d(threshold5[:,0], threshold5[:,1], kind='cubic')
interp6 = interp1d(threshold6[:,0], threshold6[:,1], kind='cubic')

ax2.plot(threshold4[:,0], threshold4[:,1], '.', color='r')
ax2.plot(threshold5[:,0], threshold5[:,1], '.', color='g')
ax2.plot(threshold6[:,0], threshold6[:,1], '.', color='b')

ax2.plot(threshold4[:,0], interp4(threshold4[:,0]), color='r', label=r'$k x_0 = 0.1$')
ax2.plot(threshold5[:,0], interp5(threshold5[:,0]), color='g', label=r'$k x_0 = 1.0$')
ax2.plot(threshold6[:,0], interp6(threshold6[:,0]), color='b', label=r'$k x_0 = 10.0$')

ax2.legend(loc=1, fontsize=20)

ax2.set_xlim([0, 10])
ax2.set_ylim([0.6, 0.97])

ax2.set_xlabel(r'$\rho_1 / \rho_0$', fontsize=20)

ax2.set_xticks([0, 2, 4, 6, 8, 10])
ax2.set_xticklabels([0, 2, 4, 6, 8, 10], fontsize=15)
ax2.set_yticks([0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95])
ax2.set_yticklabels([0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95], fontsize=15)