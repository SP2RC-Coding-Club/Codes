# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 16:48:53 2017

@author: Matt

Move seed points by displacement field
"""
from scipy.optimize import fsolve
import slab_functions as sf
import numpy as np

def move_seeds(seeds, disp_x, disp_y, disp_z, mins, maxes):
    """
    Move seed points by displacement field
    
    Parameters
    ----------
    seeds: list of tuples
        original seed points, e.g. [(x,y,z),(x,y,z)]
        
    disp_x, disp_y, disp_z: 3d array
        x, y, z displacement
    
    mins, maxes: list
        mins and maxes of x, y, z coords. [xmin, ymin, zmin]
        
    """
    
    moved_seeds = []
    for seed in seeds:
        seed = list(seed)
        seed[0] = seed[0] + disp_x[seed[0],seed[1],seed[2]] * (disp_x.shape[0] / (maxes[0]-mins[0]))
        seed[1] = seed[1] + disp_z[seed[0],seed[1],seed[2]] * (disp_z.shape[2] / (maxes[2]-mins[2]))
        seed[2] = seed[2] + disp_y[seed[0],seed[1],seed[2]] * (disp_y.shape[1] / (maxes[1]-mins[1]))
        seed = tuple(seed)
        moved_seeds.append(seed)
        
    return moved_seeds
    
def move_seeds_step(seeds, mins, maxes, tmin, tmax, n, nt,
                      mode, x, z, t, W, K, R1):
    # IGNORE THIS
    """
    Move seed points by displacement field - previous displacement field
    
    Parameters
    ----------
    seeds: list of tuples
        original seed points, e.g. [(x,y,z),(x,y,z)]
    
    mins, maxes: list
        mins and maxes of x, y, z coords. [xmin, ymin, zmin]
    """
    
    def rescale_n(unscaled, axis):
        return unscaled * ((maxes[axis]-mins[axis]) / n[axis]) + mins[axis]
    def rescale_maxmin(unscaled, axis, minus_min=False):
        if minus_min == True:
            return (unscaled - mins[axis]) * n[axis] / (maxes[axis]-mins[axis])
        else:
            return unscaled * n[axis] / (maxes[axis]-mins[axis])    
            

    
    moved_seeds = []
    for seed in seeds:
        disp_1 = [np.real(sf.xix(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=2), t, W, K, R1)),
                  np.real(sf.xiz(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=2), t, W, K, R1)),
                  np.real(sf.xiy(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=2), t, W, K))]
        disp_2 = [np.real(sf.xix(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=2), t-(tmax-tmin)*nt, W, K, R1)),
                  np.real(sf.xiz(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=2), t-(tmax-tmin)*nt, W, K, R1)),
                  np.real(sf.xiy(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=2), t-(tmax-tmin)*nt, W, K))]        
        
        new_seed = list(seed)
#        print('old seed = ' + str(new_seed))  
#        print(np.real(sf.xix(mode, x[seed[0]], z[seed[1]], t, W, K, R1)))
#        print(np.real(sf.xix(mode, x[seed[0]], z[seed[1]], t, W, K, R1)) * (n[0] / (maxes[0]-mins[0])))
        new_seed[0] = seed[0] + rescale_maxmin(disp_1[0] - disp_2[0], axis=0)
        new_seed[1] = seed[1] + rescale_maxmin(disp_1[1] - disp_2[1], axis=2)
        new_seed[2] = seed[2] + rescale_maxmin(disp_1[2] - disp_2[2], axis=1)
        
#        print('seed[1] = ' + str(seed[1]))     
#        print('disp[0] = ' + str(rescale_maxmin(disp_1[0] - disp_2[0], axis=0)))  
#        print('disp[1] = ' + str(rescale_maxmin(disp_1[1] - disp_2[1], axis=2)))  
#        print('disp[2] = ' + str(rescale_maxmin(disp_1[2] - disp_2[2], axis=1)))  
        
        moved_seeds.append(tuple(new_seed))
#        print('moved seed = ' + str(moved_seeds))
        print('\n\n')
    return moved_seeds
    
H = 1.
def move_seeds_non_int(seeds, mins, maxes, n, 
                      mode, x, z, t, W, K, R1):
    # CANT GET THIS TO WORK
    """
    Move seed points by displacement field
    
    Parameters
    ----------
    seeds: list of tuples
        original seed points, e.g. [(x,y,z),(x,y,z)]
        
    disp_x, disp_y, disp_z: functions
        x, y, z displacement
    
    mins, maxes: list
        mins and maxes of x, y, z coords. [xmin, ymin, zmin]
    """
    
    def rescale_n(unscaled, axis):
        return unscaled * (maxes[axis]-mins[axis]) / n[axis] + mins[axis]
    def rescale_maxmin(unscaled, axis, minus_min=False):
        if minus_min == True:
            return (unscaled - mins[axis]) * n[axis] / (maxes[axis]-mins[axis])
        else:
            return unscaled * n[axis] / (maxes[axis]-mins[axis])    
    
    
    moved_seeds = []
    for seed in seeds:
        new_seed = list(seed)
#        print(np.real(sf.xix(mode, x[seed[0]], z[seed[1]], t, W, K, R1)))
#        print(np.real(sf.xix(mode, x[seed[0]], z[seed[1]], t, W, K, R1)) * (n[0] / (maxes[0]-mins[0])))
        new_seed[0] = seed[0] + H * rescale_maxmin(np.real(sf.xix(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=1), t, W, K, R1)), axis=0)
        new_seed[1] = seed[1]# + H * rescale_maxmin(np.real(sf.xiz(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=1), t, W, K, R1)), axis=2)
        new_seed[2] = seed[2] + H * rescale_maxmin(np.real(sf.xiy(mode, rescale_n(seed[0], axis=0), rescale_n(seed[1], axis=1), t, W, K)), axis=1)      
        moved_seeds.append(tuple(new_seed))
        
    return moved_seeds
    
def original_seeds_non_int(moved_seeds, mins, maxes, n, 
                           mode, x, z, t, W, K, R1):
    # CANT GET THIS TO WORK
    """
    Find original seed points
    
    Parameters
    ----------
    moved_seeds: list of tuples
        moved seed points, e.g. [(x,y,z),(x,y,z)]
        
    disp_x, disp_y, disp_z: 3d array
        x, y, z displacement
    
    mins, maxes: list
        mins and maxes of x, y, z coords. e.g. [xmin, ymin, zmin]
        
    """
    

    
    def rescale_n(unscaled, axis):
        return unscaled * (maxes[axis]-mins[axis]) / n[axis] + mins[axis]
    def rescale_maxmin(unscaled, axis, minus_min=False):
        if minus_min == True:
            return (unscaled - mins[axis]) * n[axis] / (maxes[axis]-mins[axis])
        else:
            return unscaled * n[axis] / (maxes[axis]-mins[axis])
    
    
    seeds = []
    for seed in moved_seeds:
        seed = list(seed)
        def function(orig_seed):
            return [orig_seed[0] - seed[0] + H * rescale_maxmin(np.real(sf.xix(mode, rescale_n(orig_seed[0], axis=0), rescale_n(orig_seed[1], axis=1), t, W, K, R1)), axis=0),
                    orig_seed[1] - seed[1], # + H * rescale_maxmin(np.real(sf.xiz(mode, rescale_n(orig_seed[0], axis=0), rescale_n(orig_seed[1], axis=1), t, W, K, R1)), axis=2),
                    orig_seed[2] - seed[2] + H * rescale_maxmin(np.real(sf.xiy(mode, rescale_n(orig_seed[0], axis=0), rescale_n(orig_seed[1], axis=1), t, W, K)), axis=1)]
        original_seed = list(np.real(fsolve(function, seed, xtol=1e-03)))
        seeds.append(tuple(original_seed))
    return seeds
    