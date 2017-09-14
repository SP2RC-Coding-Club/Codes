#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 16:51:30 2017

@author: rach
"""
import numpy as np
import os
import sys
import random 
#import random
cwd = os.getcwd()
sys.path.insert(0, cwd + "/../../../Programming/Read_data/")
sys.path.insert(0, cwd + "/../../../Programming/Geometry_Functions/")
sys.path.insert(0, cwd + "/../../../Programming/Extra_functions/")
from functions_for_winds import alphas_from_ratios, get_min_beta_vec, return_nearest_values, alphas_from_betas, get_min_beta, alphas_from_betas_with_min, mag_between_winds

def sech(x):
    return np.cosh(x)**(-1)

def check(phi1, ratio1, beta, phi2, ratio2):
    ratio = ratio2
    if ~np.isreal(ratio2) :
        check = False
        winds_2 = False
    else:
        alphas1 = alphas_from_betas(ratio1, beta)
        winds_1_p = phi1+alphas1
        winds_1_n = phi1-alphas1
        beta_2_min = get_min_beta(ratio2)

        if beta_2_min > beta:
            check = False
            winds_2 = False
        else:
            alphas2 = alphas_from_betas_with_min(ratio2, beta, beta_2_min)
            winds_2_p = phi2+alphas2
            winds_2_n = phi2-alphas2
            if (np.round(winds_1_p,2)==np.round(winds_2_p,2)) | (np.round(winds_1_p,2)==np.round(winds_2_n,2)): 
                check = True
                winds_2 = winds_1_p
            elif (np.round(winds_1_n,2)==np.round(winds_2_p,2)) | (np.round(winds_1_n,2)==np.round(winds_2_n,2)):
                check = True
                winds_2 = winds_1_n                
            else:
                check = False
                winds_2 = False

    return [check, winds_2, ratio] 

def get_ratios_from_r1_and_beta(R1, beta, phi1, phi2): 
    A_pp = np.exp(2*beta*(phi1-phi2))*((1-10**(R1/20)*np.exp(-beta*np.pi))/(10**(R1/20)*np.exp(beta*np.pi)-1))
    R2_pp = 20*np.log10((1+A_pp)/(A_pp*np.exp(beta*np.pi)+np.exp(-beta*np.pi)))
    A_nn = np.exp(2*beta*(phi1-phi2))*((10**(R1/20)*np.exp(beta*np.pi)-1)/(1-10**(R1/20)*np.exp(-beta*np.pi)))
    R2_nn = 20*np.log10((1+A_nn)/(A_nn*np.exp(-beta*np.pi)+np.exp(beta*np.pi)))
    A_pn = np.exp(2*beta*(phi1-phi2))*((10**(R1/20)*np.exp(beta*np.pi)-1)/(1-10**(R1/20)*np.exp(-beta*np.pi)))
    R2_pn = 20*np.log10((1+A_pn)/(A_pn*np.exp(beta*np.pi)+np.exp(-beta*np.pi)))
    A_np = np.exp(2*beta*(phi1-phi2))*((1-10**(R1/20)*np.exp(-beta*np.pi))/(10**(R1/20)*np.exp(beta*np.pi)-1))
    R2_np = 20*np.log10((1+A_np)/(A_np*np.exp(-beta*np.pi)+np.exp(beta*np.pi)))
    
    ratios_possible = [R2_pp, R2_nn, R2_pn, R2_np]
    
    ratios_picked = []
    winds_picked = []
    for R2 in ratios_possible:
        [truth, winds_2, ratio_passed] = check(phi1, R1, beta, phi2, R2)        
        if truth:
            ratios_picked.append(ratio_passed)
            winds_picked.append(winds_2)
        
    return [ratios_picked, winds_picked]

def rand_between_2_5grid_ind_fuzzy(min_ratio, max_ratio, min_beta, max_beta):
    ratio_ind = random.uniform(min_ratio, max_ratio)
    beta_ind = random.uniform(min_beta, max_beta)
    return [ratio_ind,beta_ind]

def evalCosts_5grid_ind_fuzzy_xtra(individual, ratios, ratio_mins, ratio_maxs, phis):
#    print individual
    wind = np.nan
    R1 = individual[0]
    beta = individual[1]
    print [R1, beta]
#    betas[0,0], betas[0,2], betas[2,0], betas[2,2] = np.nan, np.nan, np.nan, np.nan 
#    divisor = np.count_nonzero(~np.isnan(betas)) - 1
    
#    get the other ratio values     
    phis_neighbours = np.array([phis[0,1],phis[1,2],phis[2,1],phis[1,0]])
    ratios_original = np.array([ratios[0,1],ratios[1,2],ratios[2,1],ratios[1,0]])
    
    ratios_possible, winds_possible, orig_ratio = [],[],[]
    for i, phi2 in enumerate(phis_neighbours):
        [ratios_neighbour, winds_picked] = get_ratios_from_r1_and_beta(R1, beta, phis[1,1], phi2)
        orig_ratio.append(len(ratios_neighbour)*[ratios_original[i]])
        ratios_possible.append(ratios_neighbour)
        winds_possible.append(winds_picked)
    
    flat_rats,flat_rats_orig,flat_winds = np.ravel(np.asarray(ratios_possible),'c'),np.ravel(np.asarray(orig_ratio),'c'),np.ravel(np.asarray(winds_possible),'c')
        
    # keep the better option     
    if np.asarray(ratios_possible).size != 8:
        J_dist = 40
    else:
        wind_choices = np.unique(flat_winds)         
        diffs_1 = np.abs(flat_rats[flat_winds==wind_choices[0]] - flat_rats_orig[flat_winds==wind_choices[0]])
        diffs_2 = np.abs(flat_rats[flat_winds==wind_choices[1]] - flat_rats_orig[flat_winds==wind_choices[0]])
        if np.sum(diffs_1) < np.sum(diffs_2):
            J_dist = 0
            wind = np.rad2deg(wind_choices[0])
            if ((flat_rats[flat_winds==wind_choices[0]]<(flat_rats_orig[flat_winds==wind_choices[0]]-3)) | (flat_rats[flat_winds==wind_choices[0]]>(flat_rats_orig[flat_winds==wind_choices[0]]+3))).any():
                J_dist = 10*np.sum(((flat_rats[flat_winds==wind_choices[0]]<(flat_rats_orig[flat_winds==wind_choices[0]]-3)) | (flat_rats[flat_winds==wind_choices[0]]>(flat_rats_orig[flat_winds==wind_choices[0]]+3))))                
        else:
            J_dist = 0
            wind = np.rad2deg(wind_choices[1])
            
    return (J_dist, np.abs(R1-ratios[1,1]), wind)
