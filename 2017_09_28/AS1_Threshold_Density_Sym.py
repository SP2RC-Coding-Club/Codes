import solvers as sol
from AS1_class import Asym_slab

import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d

def save():
    with open('pickles/threshold_density_sym_c0={}_K={}_M_A={}.p'.format(
    				slab.c0, slab.K, slab.M_A), 'wb') as f:
    	    pickle.dump(threshold, f)

def load(K):
    return pickle.load(open('pickles/threshold_density_sym_c0={}_K={}_M_A={}.p'
                            .format(c0, K, slab.M_A), 'rb'))

#R_range = np.append([0.05], np.linspace(0.25, 10, 79))
R_range = np.linspace(7.125, 10, 24)
K = 10.0
c0 = 0.6

slab = Asym_slab(c0=c0, R1=R_range[0], R2=R_range[0], K=K, M_A=None)

"""First attempt at finding the first complex value for R_range[0].
Performed over a larger than usual linspace to make sure that it's found."""

x_range = np.linspace(0.7, 0.75, 200)
y_range = np.linspace(0, slab.c1, 200)

loc = sol.find_first_imag(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))
loc = np.round(loc, 6)

"""Second attempt at finding the first complex value for R_range[0].
Performed over a smaller linspace in order to find an accurate value."""

x_range = np.linspace(loc[0]-0.02, loc[0]+0.02, 200)
y_range = np.linspace(max([0, loc[1]-0.02]), loc[1]+0.02, 200)

loc = sol.find_first_imag(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))
loc = np.round(loc, 6)

print('Found the threshold for R={} at x={}, for y={}. \n'.format(round(R_range[0],2), loc[0], loc[1]))

threshold = np.swapaxes(np.vstack((np.array([R_range[0]]), [loc[0]])),1,0)
roots = np.swapaxes(np.vstack((np.array([R_range[0]]), [loc[1]])), 1, 0)

space_size = 0.01

debug_loop_count = 0

"""This loops over all other values in R_range, finds the first complex value for each element
and stacks it under threshold. Column 1 contains the values where the first complex entry was
found, and column 0 contains the value of R that it was found for."""

for R in R_range[1:]:

    slab = Asym_slab(c0=c0, R1=R, R2=R, K=K, M_A=None)

    loc = None

    numpoints = 500

    x_range = np.real(np.linspace(threshold[-1,1]-space_size,
                                  threshold[-1,1], numpoints))
    y_range = np.real(np.linspace(max([0, roots[-1, 1]-space_size]),
                                  roots[-1, 1], numpoints/2))

    loc = sol.find_first_imag(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))

    while loc is None or loc[0] > threshold[-1, 1]:

        debug_loop_count+=1

        space_size *= 1.5

        print('!!! Could not find threshold value for R1={}. Increasing space_size to {}'
              .format(round(R,2), space_size))

        x_range = np.linspace(threshold[-1,1]-space_size,
                              threshold[-1,1]+space_size*0.3, numpoints)
        y_range = np.linspace(max([0, roots[-1, 1]-space_size]),
                              roots[-1, 1]+space_size*0.3, numpoints)

        loc = sol.find_first_imag(slab.disp_rel, x_range, y_range, args=(slab.K, slab.M_A))

    loc = np.round(loc, 6)

    threshold = np.vstack((threshold, [R, loc[0]]))
    roots = np.vstack((roots, [R, loc[1]]))

    print('Found the threshold for R1={} at x={}, for y={}.\
          \nDifference between last two entries of {} on the x-axis, and {} on the y-axis. \n'
          .format(round(R,2), loc[0], loc[1],
          np.round(threshold[-2,1]-threshold[-1,1], 6), np.round(roots[-2,1]-roots[-1,1], 6)))

    np.save('pickles/threshold_density_sym_temp_c0={}_K={}_M_A={}.npy'.format(
    				slab.c0, slab.K, slab.M_A), threshold)

    """We adjust the size of the linspace to get better accuracy."""

    space_size = 0.01 + np.maximum(np.abs(threshold[-2,1]-threshold[-1,1]),
                                   np.abs(roots[-2,1]-roots[-1,1]))
threshold = np.real(threshold)
f = interp1d(threshold[:,0], threshold[:,1], kind='cubic')

plt.figure(num=None, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
plt.plot(threshold[:,0], threshold[:,1], '.', color='b')
plt.plot(threshold[:,0], f(threshold[:,0]), color='b')