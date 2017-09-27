import solvers as sol
from AS1_class import Asym_slab

import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d

R_range = np.append([0.05], np.linspace(0.25, 10, 79))
R2 = 2
K = 10.0
c0 = 0.6

slab = Asym_slab(c0=c0, R1=R_range[0], R2=R2, K=K, M_A=None)

"""First attempt at finding the first complex value for R_range[0].
Performed over a larger than usual linspace to make sure that it's found."""

x_range = np.linspace(0, 2, 201)
y_range = np.linspace(0, slab.c1, 201)

loc = sol.find_first_imag(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))
loc = np.round(np.real(loc), 6)

"""Second attempt at finding the first complex value for R_range[0].
Performed over a smaller linspace in order to find an accurate value."""

x_range = np.linspace(loc[0]-0.1, loc[0]+0.1, 201)
y_range = np.linspace(0, slab.c1, 201)

loc = sol.find_first_imag(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))
loc = np.round(np.real(loc), 6)

threshold = np.swapaxes(np.vstack((np.array([R_range[0]]), [loc[0]])),1,0)

space_size = 0.1

"""This loops over all other values in R_range, finds the first complex value for each element
and stacks it under threshold. Column 1 contains the values where the first complex entry was
found, and column 0 contains the value of R1 that it was found for."""

space_size = 0.03

for R1 in R_range[16:]:

    x_range = np.linspace(loc[0]-space_size, loc[0]+space_size, 201)
    y_range = np.linspace(0, slab.c1, 151)

    slab = Asym_slab(c0=c0, R1=R1, R2=R2, K=K, M_A=None)

    loc = sol.find_first_imag(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))
    loc = np.round(np.real(loc), 6)

    threshold = np.vstack((threshold, [R1, loc[0]]))

    """We adjust the size of the linspace to get better accuracy."""

    space_size = 0.02 + (threshold[-2,1] - threshold[-1,1]) * 1.1

threshold = np.real(threshold)
f = interp1d(threshold[:,0], threshold[:,1], kind='cubic')

plt.figure(num=None, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
plt.plot(threshold[:,0], threshold[:,1], '.', color='b')
plt.plot(threshold[:,0], f(threshold[:,0]), color='b')

def save():
    with open('pickles/threshold_density_c0={}_R2={}_K={}_M_A={}.p'.format(
    				slab.c0, slab.R2, slab.K, slab.M_A), 'wb') as f:
    	    pickle.dump(threshold, f)

def load(K):
    return pickle.load(open('pickles/threshold_density_c0={}_R2={}_K={}_M_A={}.p'.format(
            slab.c0, slab.R2, K, slab.M_A), 'rb'))