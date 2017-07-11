# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:52:53 2017

@author: Matt
"""

import numpy as np
import slab_functions as sf
import toolbox as tool
import matplotlib.pyplot as plt
import matplotlib


def density_diagram(disp_rel, K, W, R1, R1min, R1max,
                    just_dots=False):
    # Plot omega/k against rho_1/rho_0 for eigenmodes
                    
    Wmin = sf.cT
    Wmax = sf.vA
    
    R1_range = np.linspace(R1min, R1max, 51)
    W_range = np.linspace(Wmin, Wmax, 51)
    
    # Global font change
    font = {'size': 15}
    matplotlib.rc('font', **font)
    
    plt.figure(num=None, figsize=(10, 11), facecolor='w', edgecolor='k')
    ax = plt.subplot()
        
    if just_dots == True:
        #Plot the dots
        W_array = tool.point_find(disp_rel, R1_range, W_range,
                                          args=(None))
        ax.plot(R1_range, W_array, '.', color = 'b')
        ax.set_xlim(0, R1_range[-1])
        ax.set_ylim(Wmin, Wmax)


def dispersion_diagram(mode_options, chosen_mode, disp_rel, K, W, R1,
                       just_dots=False):
    
    # Dispersion diagram for eigenmodes given by mode_options.
    # Can only do R1=1.5, 1.8, 2.0.
    
    # Create sublists of modes
    kink_mode_options = []
    saus_mode_options = []
    slow_surf_mode_options = []
    fast_surf_mode_options = []
    fast_kink_mode_options = []
    fast_saus_mode_options = []
    slow_body_1_mode_options = []
    slow_body_2_mode_options = []
    slow_body_3_mode_options = []
    slow_body_mode_options = []
    fast_body_1_mode_options = []
    fast_body_2_mode_options = []
    fast_body_3_mode_options = []
    
    for mode in mode_options:
        if 'kink' in mode:
            kink_mode_options.append(mode)
        if 'saus' in mode:
            saus_mode_options.append(mode)
        if 'slow' in mode and 'surf' in mode:
            slow_surf_mode_options.append(mode)
        if 'fast' in mode and 'surf' in mode:
            fast_surf_mode_options.append(mode)
        if 'fast' in mode and 'kink' in mode:
            fast_kink_mode_options.append(mode)
        if 'fast' in mode and 'saus' in mode:
            fast_saus_mode_options.append(mode)
        if 'fast' in mode and 'body-1' in mode:
            fast_body_1_mode_options.append(mode)
        if 'fast' in mode and 'body-2' in mode:
            fast_body_2_mode_options.append(mode)
        if 'fast' in mode and 'body-3' in mode:
            fast_body_3_mode_options.append(mode)
        if 'slow' in mode and 'body-1' in mode:
            slow_body_1_mode_options.append(mode)
        if 'slow' in mode and 'body-2' in mode:
            slow_body_2_mode_options.append(mode)
        if 'slow' in mode and 'body-3' in mode:
            slow_body_3_mode_options.append(mode)
        if 'slow' in mode and 'body' in mode:
            slow_body_mode_options.append(mode)  
    
    # Set plotting range
    Kmin = 0.
    Kmax = 23.#10.
    
    Wmin = 0.
    Wmax = sf.c2
    
    K_range = np.linspace(Kmin, Kmax, 51)
    W_range = np.linspace(Wmin, Wmax, 51)
    
    # Global font change
    font = {'size': 15}
    matplotlib.rc('font', **font)
    
    # Initialise plot
    plt.figure(num=None, figsize=(10, 11), facecolor='w', edgecolor='k')
    ax = plt.subplot()
    
    # line colours - g=chosen, b=body, r=surface
    def colour(mode):
        colours = []
        if mode in fast_surf_mode_options:
            for i in range(len(mode_options)):
                if i <= 7:
                    if mode_options[i] == mode:
                        colours.append('g')
                    elif mode_options[i] in slow_surf_mode_options + fast_surf_mode_options:
                        colours.append('r')
                    else:
                        colours.append('b')
                elif i >= 14:
                    if mode_options[i] == mode:
                        colours.append('g')
                    elif mode_options[i] in slow_surf_mode_options + fast_surf_mode_options:
                        colours.append('r')
                    else:
                        colours.append('b')
        else:
            for i in range(len(mode_options)):
                if mode_options[i] == mode:
                    colours.append('g')
                elif mode_options[i] in slow_surf_mode_options + fast_surf_mode_options:
                    colours.append('r')
                else:
                    colours.append('b')
        return colours
     
    # line styles - '--'=kink, '-'=sausage
    def ls(mode):
        linestyles = []
        if mode in fast_surf_mode_options:
            for i in range(8):
                if mode_options[i] in kink_mode_options:
                    linestyles.append('--')
                else:
                    linestyles.append('-')
            for i in range(8, len(mode_options)):
                if mode_options[i-1] in kink_mode_options:
                    linestyles.append('--')
                else:
                    linestyles.append('-')
        else:
            for i in range(len(mode_options)):
                if mode_options[i] in kink_mode_options:
                    linestyles.append('--')
                else:
                    linestyles.append('-')
        return linestyles
    
    # Basic, but much fast plot to get an idea of the solutions.
    if just_dots == True:
        #Plot the dots
        W_array = tool.point_find(disp_rel, K_range, W_range,
                                          args=(None))
        ax.plot(K_range, W_array, '.', color = 'b')
        
    else:
    # Do the full plot
        if R1 == 1.5:
            if chosen_mode in fast_surf_mode_options:
                K_guess = [1., 1., 10.12, 8.74, 5.98, 5.98, 3.68, 1.84, 6., 6.]
                W_guess = [0.64, 0.72, 0.863, 0.867, 0.855, 0.885, 0.885,
                           0.920, 1.01567, 1.05754]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01]
                K_start = [0.0001, 0.0001, 2., 1.6, 1.2, 0.9, 0.6, 0.2, 4.22, 
                           0.67]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, 4.2, Kmax, 
                         Kmax]
                K_fast_body_circles = [-10, 0.67]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = [4.2]
                W_transition_circles = [sf.c0]
            else:
                K_guess = [1., 1., 8.74, 8.74, 5.98, 5.98, 3.68, 1.84, 1.84,
                           6.44, 9.2, 15.18, 19.78, 22.01]
                W_guess = [0.473, 0.547, 0.719463, 0.733764, 0.722183, 
                           0.746659, 0.741942, 0.749326, 1.085, 1.131, 1.181,
                           1.154, 1.156, 1.181]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
                K_start = [0.0001, 0.0001, 2.3, 1.8, 1.3, 1., 0.8, 0.3, 0.47, 
                           4.51, 8.55, 12.55, 16.5, 20.53]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax,
                         Kmax, Kmax, Kmax, Kmax, Kmax]
                K_fast_body_circles = [0.47, 4.51, 8.55, 12.55, 16.5, 20.53]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = []
                W_transition_circles = []
     
                
        elif R1 == 1.8:
            if chosen_mode in fast_surf_mode_options:
                K_guess = [1., 1., 10.12, 8.74, 5.98, 5.98, 3.68, 1.84, 6., 6.]
                W_guess = [0.64, 0.72, 0.863, 0.867, 0.855, 0.885, 0.885,
                           0.920, 1.01567, 1.05754]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01]
                K_start = [0.0001, 0.0001, 2., 1.6, 1.2, 0.9, 0.6, 0.2, 4.752, 
                           0.365]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, 4.752, Kmax, 
                         Kmax]
                K_fast_body_circles = [-10, 0.365]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = [4.752]
                W_transition_circles = [sf.c0]
            else:
                K_guess = [1., 1., 8.74, 8.74, 5.98, 5.98, 3.68, 1.84, 1.84,
                           6.44, 9.2, 15.18, 19.78, 22.01]
                W_guess = [0.473, 0.547, 0.719463, 0.733764, 0.722183, 
                           0.746659, 0.741942, 0.749326, 1.085, 1.131, 1.181,
                           1.154, 1.156, 1.181]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
                K_start = [0.0001, 0.0001, 2.3, 1.8, 1.3, 1., 0.8, 0.3, 0.336,
                           4.306, 8.346, 12.318, 16.309, 20.342]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax,
                         Kmax, Kmax, Kmax, Kmax, Kmax]
                K_fast_body_circles = [0.336, 4.306, 8.346, 12.318, 16.309,
                                       20.342]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = []
                W_transition_circles = []
    
        elif R1 == 2.:
            if chosen_mode in fast_surf_mode_options:
                K_guess = [1., 1., 11.96, 6., 4., 4., 2., 2., 6., 6.]
                W_guess = [0.64, 0.72, 0.877, 0.836, 0.826, 0.847, 0.83, 
                           0.915, 1.01567, 1.05754]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01]
                K_start = [0.0001, 0.0001, 2., 1.6, 1.2, 0.9, 0.6, 0.2, 5.24, 
                           0.001]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, 5.2, Kmax, 
                         Kmax]
                K_fast_body_circles = []
                W_fast_body_circles = []
                K_transition_circles = [5.2]
                W_transition_circles = [sf.c0]
            else:
                K_guess = [1., 1., 8.74, 8.74, 5.98, 5.98, 3.68, 1.84, 2., 6., 
                           9., 14.73, 19.8, 23.4]
                W_guess = [0.473, 0.547, 0.719463, 0.733764, 0.722183, 
                           0.746659, 0.741942, 0.749326, 1.16675, 1.156, 
                           1.17, 1.1547, 1.16, 1.16]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
                K_start = [0.0001, 0.0001, 2.3, 1.8, 1.3, 1., 0.8, 0.3, 0.0001, 
                           4.1, 8.2, 12.2, 16.2, 20.2]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax,
                         Kmax, Kmax, Kmax, Kmax, Kmax]
                K_fast_body_circles = [-10, 4.1, 8.2, 12.2, 16.2, 20.2]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = []
                W_transition_circles = []
        else:
            raise ValueError('R1 can only take the values 1.5, 1.8, 2.0')
        
        # line width
        lw = 1.5
        disp_modes = range(len(K_guess))
        for i in disp_modes:
            K_values, root_array = tool.line_trace(disp_rel, K_guess[i], 
                                                         W_guess[i], step[i], K_start[i],
                                                         K_end[i], (None))
            ax.plot(K_values, root_array, color=colour(chosen_mode)[i], linewidth=lw, 
                    linestyle=ls(chosen_mode)[i])
        for i in range(len(K_fast_body_circles)):
            ax.plot(K_fast_body_circles[i], W_fast_body_circles[i], marker='o', 
                    markerfacecolor='None', markeredgecolor=colour(chosen_mode)[i+8], 
                    markeredgewidth=lw, markersize=8)
        for i in range(len(K_transition_circles)):
            ax.plot(K_transition_circles[i], W_transition_circles[i], marker='o',
                    markerfacecolor='None', markeredgecolor=colour(chosen_mode)[i+8],
                    markeredgewidth=lw, markersize=8)
        
        # Labels, annotations etc...
        ax.set_ylabel(r'$v_{ph}/c_0$', fontsize = 20)
        ax.set_xlabel(r'$k x_0$', fontsize = 20)
        
        ax.set_xlim(K_range[0], K_range[-1])
        ax.set_ylim(0., 1.41)
        
        if chosen_mode in fast_surf_mode_options:
            ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c2, sf.c2), (1.42, 1.42),
                edgecolor='gray', linestyle='-.', color='None', hatch='/',
                linewidth=2)
    
            #ax.plot([K_range[0], K_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='--', linewidth=2)
            ax.annotate(r'$v_A$', xy=(K_range[-1] + 0.03, sf.vA - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='-.', linewidth=2)
            ax.plot([K_range[0], K_range[-1]], [sf.cT, sf.cT], color = '0.5', linestyle='-.', linewidth=2)
            ax.annotate(r'$c_T$', xy=(K_range[-1] + 0.03, sf.cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.c0, sf.c0], color = '0.5', linestyle='-.', linewidth=2)
            ax.annotate(r'$c_0$', xy=(K_range[-1] + 0.03, sf.c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            #ax.plot([K_range[0], K_range[-1]], [sf.c2, sf.c2], color = '0.5', linestyle='--', linewidth=2)
    #            ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='-.', linewidth=2)
            if R1 == 2.:
                ax.annotate(r'$c_1=c_2$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            else:
                ax.annotate(r'$c_1$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                
        else:
            K_range_for_fill = np.append(K_range, K_range[-1]+0.1)
            ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c0,sf.c0), (sf.vA,sf.vA),
                            edgecolor='gray', linestyle='-.', color='None', hatch='/',
                            linewidth=2)
            ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c2, sf.c2), (1.42, 1.42),
                            edgecolor='gray', linestyle='-.', color='None', hatch='/',
                            linewidth=2)
            
            #ax.plot([K_range[0], K_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='--', linewidth=2)
            ax.annotate(r'$v_A$', xy=(K_range[-1] + 0.03, sf.vA - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.cT, sf.cT], color = '0.5', linestyle='-.', linewidth=2)
            ax.annotate(r'$c_T$', xy=(K_range[-1] + 0.03, sf.cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            #ax.plot([K_range[0], K_range[-1]], [sf.c0, sf.c0], color = '0.5', linestyle='--', linewidth=2)
            ax.annotate(r'$c_0$', xy=(K_range[-1] + 0.03, sf.c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            #ax.plot([K_range[0], K_range[-1]], [sf.c2, sf.c2], color = '0.5', linestyle='--', linewidth=2)
    #            ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='-.', linewidth=2)
            if R1 == 2.:
                ax.annotate(r'$c_1=c_2$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            else:
                ax.annotate(r'$c_1$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
    
        # Plot a green circle where the chosen mode is
        ax.plot(K, W, 'go', markersize=10)
    #        plt.tight_layout() # seems to make it chop the sides off with this