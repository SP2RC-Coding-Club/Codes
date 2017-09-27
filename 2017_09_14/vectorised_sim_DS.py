#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 11:17:48 2017

Vectorised version of the doppler spectrum simulation program.

@author: rach
"""
from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.integrate as integrate
import time
import os
import sys
cwd = os.getcwd()
sys.path.insert(0, cwd + "/../Geometry_Functions/")
from scatter_geometry import ScatterGeometry_rotate 

#np.set_printoptions(threshold=np.nan)

# %% function section 

def sech(x):
    return np.cosh(x)**(-1)

def csch(x):
    return np.sinh(x)**(-1)

def mag_between_winds(ang_1,ang_2):
    """
        get the smallest magnitude between two angles 
        inputs:
            angle 1 - radians 
            angle 2 - radians
        Returns:
            ang_dist - angle between them in radians 
    """    
    ang_prelim = ang_1 - ang_2
    ang_dist = np.abs(np.mod(ang_prelim + np.pi, 2*np.pi) - np.pi)    
    return ang_dist

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """
    return np.isnan(y), lambda z: z.nonzero()[0]

#%%
def PMSpecII(theta, thetaStar, K, U, s):
    """
    Pierson Moskowitz Spectrum to calculate the power density of the spectrum at a certain frequency.
    A cosine directional spectrum is used to get the power of the PM spec in a certain direction 
        Inputs:
            theta - 'look' direction (rads)*
            thetaStar - wind direction (rads)*
            K - wave number of the wave
            U - wind speed m/s
            *doesn't matter on orientation as long as they are in the same frame of ref
        Returns:
            SpectrumValue - what it says on the tin. 
    """

#   calculate the normalised propotion of the power in the desired direction
    G = (np.abs(np.cos((mag_between_winds(theta,thetaStar))/2)))**(2*s)   
        
    DirectionIntegralpre = lambda theta: (np.cos(theta/2)+0j)**(2*s)
    DirectionIntegral = integrate.quad(DirectionIntegralpre,0,2*np.pi)[0]

#   calculate the power at the desired frequency (well wavenumber)    
    alpha, beta, g = 0.0081, 0.74, 9.807 # constants for the PM spec 
    SpectrumValue = ((alpha/(2*K**4))*np.exp(-((beta*g**2)/(K**2*U**4))))*(G/DirectionIntegral)
    
    return SpectrumValue

def electro_cc(k1,k2,theta,phi_bi, kB, k0):
    """
        Calculate the electromagnetic coupling coefficient - look in thesis for description
        ... too long
    """        
    k1x = k1*np.cos(theta)
    k2x = -kB*np.cos(phi_bi)-k1*np.cos(theta)
    k1_dot_a = -k1*np.cos(2*phi_bi+theta)
    k2_dot_a = kB*np.cos(phi_bi)+k1*np.cos(2*phi_bi+theta)

    A = -k1x*k2_dot_a - 2*np.cos(phi_bi)**2*( -k2**2 + 2*k0*k2_dot_a )
    B = -k2x*k1_dot_a - 2*np.cos(phi_bi)**2*( -k1**2 + 2*k0*k1_dot_a )
    b1 = np.sqrt(-k2**2 + 2*k0*k2_dot_a + 0j) + 0.5*(0.11-0.12j)
    b2 = np.sqrt(-k1**2 + 2*k0*k1_dot_a + 0j) + 0.5*(0.11-0.12j)

    electro = (1/(np.cos(phi_bi)**2))*np.sqrt(((A/b1)*(A/b1 + B/b2)))
    return electro

def hydro_cc(K1,K2,depth,eta,omegaB, kB, theta, phi_bi, m,md,L):
    """
        Calculate the (shallow water) hydrodynamic coupling coefficient - look in thesis for description
        ... too long
    """    
    g = 9.807
    etaTerm = (eta**2+omegaB**2)/(eta**2-omegaB**2)
    omega1 = np.sqrt(g*K1*np.tanh(K1*depth))
    omega2 = np.sqrt(g*K2*np.tanh(K2*depth))
    hydro = ((-2*1j)*( K1*np.tanh(K1*depth)+K2*np.tanh(K2*depth) - 
        (np.sqrt(K1)*(K2*np.tanh(K1*depth)*np.tanh(K2*depth)+kB*np.cos(theta+phi_bi)+K1))*
        (etaTerm)*(1/(L*np.sqrt(K2*np.tanh(K1*depth)*np.tanh(K2*depth)))) +
         (eta/g)*((m*omega1)**3*csch(K1*depth)**2+(md*omega2)**3*csch(K2*depth)**2)/(eta**2-omegaB**2)  ))
    return hydro

def jacobian_term(K1,K2,ystar,kB,theta,phi_bi,depth,m,md):
    """
        Calculate the jacobian (consult thesis) 
    """
    g = 9.807
    dydh = (1/abs(np.sqrt(g)*
             (  m*(np.sqrt(np.tanh(K1*depth)) + (K1*depth*(sech(K1*depth))**2)/(np.sqrt(np.tanh(K1*depth)))) 
                + (md/K2**1.5)*(ystar**3+ystar*kB*np.cos(theta+phi_bi))*(np.sqrt(np.tanh(K2*depth))+K2*depth*((sech(K2*depth))**2)/(np.sqrt(np.tanh(K2*depth))))  )  
            ))
    jacobian = (2*ystar**3)*dydh
    return jacobian
 
def spec(theta, td, thetastar, K1,K2,U,s,phi_bi, L, eta):
    """
        Get the contribution of the spectral terms in the integral   
    """
    index_range = np.arange(theta.shape[0]).reshape([K1.shape[0],1])
    idx_1 = np.array(index_range[(L==1) & (eta>0)])
    idx_2 = np.array(index_range[(L==1) & (eta<0)])
    idx_3 = np.array(index_range[(L==-1) & (eta>0)])
    idx_4 = np.array(index_range[(L==-1) & (eta<0)])
    
    SP1,SP2,SP3,SP4 = np.empty((theta.shape)),np.empty((theta.shape)),np.empty((theta.shape)),np.empty((theta.shape))
    
    SP1[idx_1,:] = PMSpecII(theta[idx_1,:], thetastar, K1[idx_1,:], U, s)
    SP2[idx_1,:] = PMSpecII(td[idx_1,:], thetastar, K2[idx_1,:], U, s)
    SP3[idx_1,:] = PMSpecII(-td[idx_1,:] - 2*phi_bi, thetastar, K2[idx_1,:], U, s)
    SP4[idx_1,:] = PMSpecII(-theta[idx_1,:] - 2*phi_bi, thetastar, K1[idx_1,:], U, s)
    
    SP1[idx_2,:] = PMSpecII(theta[idx_2,:], thetastar + np.pi, K1[idx_2,:], U, s)
    SP2[idx_2,:] = PMSpecII(td[idx_2,:], thetastar + np.pi, K2[idx_2,:], U, s)
    SP3[idx_2,:] = PMSpecII(-td[idx_2,:] - 2*phi_bi, thetastar + np.pi, K2[idx_2,:], U, s)
    SP4[idx_2,:] = PMSpecII(-theta[idx_2,:] - 2*phi_bi, thetastar + np.pi, K1[idx_2,:], U, s) 
    
    SP1[idx_3,:] = PMSpecII(theta[idx_3,:], thetastar + np.pi, K1[idx_3,:], U, s)
    SP2[idx_3,:] = PMSpecII(td[idx_3,:], thetastar, K2[idx_3,:], U, s)
    SP3[idx_3,:] = PMSpecII(-td[idx_3,:] - 2*phi_bi, thetastar, K2[idx_3,:], U, s)
    SP4[idx_3,:] = PMSpecII(-theta[idx_3,:] - 2*phi_bi, thetastar + np.pi, K1[idx_3,:], U, s) 
    
    SP1[idx_4,:] = PMSpecII(theta[idx_4,:], thetastar, K1[idx_4,:], U, s)
    SP2[idx_4,:] = PMSpecII(td[idx_4,:], thetastar + np.pi, K2[idx_4,:], U, s)
    SP3[idx_4,:] = PMSpecII(-td[idx_4,:] - 2*phi_bi, thetastar + np.pi, K2[idx_4,:], U, s)
    SP4[idx_4,:] = PMSpecII(-theta[idx_4,:] - 2*phi_bi, thetastar, K1[idx_4,:], U, s)  

    return SP1,SP2,SP3,SP4

# %%               
def vectorised_sim_DS_rotate(f, thetastar, s, U, distance, R, depth, zeta, phi,etapoints, NQ, current_speed, current_angle, noise_variance):
    """
    Inputs:
        f - frequency of the radar in Mhz
        thetastar - wind direction in degrees (0˚=N, 90º=E, etc)
        s - spreading value for a cosine model 
        U - Wind speed in m/s
        distance - distance between the Tx and Rx
        R - distance from the Rx to the scatter point (same units as distance)
        depth - depth of the ocean at the scatter point (m) ?????? 
        zeta - bearing of the Rx-Tx
        phi - bearing from Rx to the scatter point (0˚=N, 90º=E, etc)
        etapoints - number of points to calculate the RCS for on eta range
        NQ - number of points in the integral for each eta point
        current_speed
        current_angle
        noise_variance
    Returns:
        eta - frequency range in Hz 
        RCS - returned RCS in dB
    """
    
# %%   ------- 
    t = time.time()
    eta = np.array(np.linspace(-6,6,etapoints)) # range of etas to calc rcs for 

#   scatter geometry to return bistatic angle eta
    [bragg_bearing, phi_bi, wind_dir, tx_sp_bearing,x,y, theta_tx_rx_rot] = ScatterGeometry_rotate(zeta, phi, thetastar, distance, R)
#   rotate the wind to get the geometry on the correct axes
    rotation_amount = - tx_sp_bearing
    rotated_wind = np.rad2deg(np.mod(wind_dir + rotation_amount, 2*np.pi))
#   alter for if the Rx is in the bottom half of the plane     
    if (theta_tx_rx_rot< 2*np.pi and theta_tx_rx_rot > np.pi):
        rotated_wind = - rotated_wind 

    thetastar = np.deg2rad(rotated_wind) # convert to radians 
    
#   Calculate the radar wave number, k0, given the radar frequency 
    SpeedOfLight = 299.792458 # speed of light in units to allow input of MHz 
    lmbda = SpeedOfLight/f # wavelength of radio wave
    k0 = 2*np.pi/lmbda # radar wave number
    g = 9.807    
        
    omegaB = np.sqrt(2*g*k0*np.cos(phi_bi)*np.tanh(2*k0*np.cos(phi_bi)*depth))
    kB = 2*k0*np.cos(phi_bi)
        
#    BETTER WAY??
    if any(eta == omegaB) is False:        
        eta = np.append(eta,omegaB)    
    if any(eta == -omegaB) is False:        
        eta = np.append(eta,-omegaB)
    
    eta = np.sort(eta)
    etapoints=eta.size
    eta = eta.reshape(etapoints,1)

    crossing_point = 2*np.sqrt(k0*g*np.cos(phi_bi)*np.tanh(k0*np.cos(phi_bi)*depth))

    u = np.abs(eta) - omegaB   # For seeing if the eta value lies between the Bragg lines or outide

    L = np.empty_like(eta)
    m = -np.ones_like(eta)
    md = -np.ones_like(eta)

    L[u>0], L[u<0]  = 1, -1 

    m[np.logical_or(eta>omegaB, (L==-1) & (eta<0) )] = 1
    md[np.logical_or(eta>omegaB, (L==-1) & (eta>0) )] = 1    

    up = np.pi*np.ones((etapoints,1))
    etas_above_crossing = eta[abs(eta) > crossing_point]
    print '1st section: '+str(time.time() - t)
    
# %%   -------     
    t = time.time()
    a0 = etas_above_crossing**2/(4*g)
    f = lambda k1: k1*np.tanh(depth*k1) - a0    
    a = fsolve(f, a0)

    up[abs(eta) > crossing_point] = np.pi - np.arccos(k0*np.cos(phi_bi)/a)
    
    y0 = np.empty((etapoints,1))
    y0[L==1] = (np.sqrt(g)/(2*eta[L==1]*m[L==1]))*(eta[L==1]**2/g-kB)      # Outside the Bragg lines
    y0[L==-1] = (m[L==-1]*eta[L==-1]/np.sqrt(g)+np.sqrt(2*kB-(eta[L==-1]**2/g)))/2   # Inside the Bragg lines
             
    Thetas = np.empty((etapoints, NQ))
    print '2nd section: '+str(time.time() - t)
    
# %%   -------    
    t = time.time()
    for i in np.linspace(0,etapoints-1,etapoints, dtype = int):
        if L[i] == 1:
            Thetas[i, :] = np.linspace(-phi_bi, -phi_bi + up[i], NQ)
        else:
            Thetas[i,:] = np.linspace(np.pi - phi_bi, -phi_bi, NQ)
    print '3rd section: '+str(time.time() - t)
    
# %%   -------
    t = time.time()
    ystar = np.empty((etapoints,NQ))
    for i in np.linspace(0,etapoints-1,etapoints, dtype = int):
        for j in np.linspace(0,NQ-1,NQ, dtype = int):
            f = lambda y:eta[i]*m[i] - np.sqrt(g)*y - L[i]*np.sqrt(g)*(y**4+2*y**2*kB*np.cos(Thetas[i,j]+phi_bi)+kB**2)**0.25
            a0 = y0[i]
            ystar[i,j] = fsolve(f, a0) 
    print '4th section: '+str(time.time() - t)
    
# %%   -------    
    # wave numbers 
    t = time.time()
    K1 = np.array(ystar**2)
    K2 = np.array(np.sqrt(K1**2+2*kB*K1*np.cos(Thetas+phi_bi)+kB**2))

    # Coupling Coefficient
    gamma_EM = electro_cc(K1,K2,Thetas,phi_bi, kB, k0)
    gammaH = hydro_cc(K1,K2,depth,eta,omegaB, kB, Thetas, phi_bi, m,md,L)
    gammaL = (abs(gammaH + gamma_EM))**2
    
    # Jacobian
    dydh_v = jacobian_term(K1,K2,ystar,kB,Thetas,phi_bi,depth,m,md)

    # theta' - angle of the long wave vector 
    td = np.array(np.real(np.arccos((kB+K1*np.cos(Thetas+phi_bi))/K2) + np.pi - phi_bi))
    # Pierson Moskowitz spectra calculations 
    [SP1,SP2,SP3,SP4] = spec(Thetas, td, thetastar, K1,K2,U,s,phi_bi, L, eta)
    SpectrumContribution = SP1*SP2 + SP3*SP4
    print '5th section: '+str(time.time() - t)
    
# %%   -------    
    t = time.time()
    InnerIntegral = (gammaL*dydh_v*SpectrumContribution)
 
    TA = InnerIntegral[:,0].reshape(etapoints,1) 
    TB = InnerIntegral[:,-1].reshape(etapoints,1)  
                              
    # Use the trapezoid rule to calculate the integral and save to the F array 
    RCS = (2**5)*np.pi*k0**4*np.cos(phi_bi)**4*(((np.sum(InnerIntegral,1)).reshape(etapoints,1)  - TA/2 - TB/2)*up/(NQ - 1)) 
    # Calculate the first order peaks 
    RCS[eta==omegaB] = 2**6*np.pi*k0**4*(np.cos(phi_bi))**4*PMSpecII(np.pi - phi_bi, thetastar, k0, U, s)
    RCS[eta==-omegaB] = 2**6*np.pi*k0**4*(np.cos(phi_bi))**4*PMSpecII(-phi_bi, thetastar, k0, U, s)
    
    # calculate the noise at each point 
    sd = np.sqrt(noise_variance)
    noise = sd*np.random.randn(len(eta),1)

    RCS = 10*np.log10(RCS+noise) # add the gaussian noise & convert to dB
    RCS[abs(eta)<0.5] = np.nan # replace the centre values with nans 
    RCS[RCS<-70] = -70 # sort a noise floor 
   
    nans, x = nan_helper(RCS) # get locations of nans for interpolation
    RCS[nans] = np.interp(x(nans), x(~nans), RCS[~nans]) # interpolate! interpolate!
    print '6th section: '+str(time.time() - t) 
    
    return [eta, RCS] 

#%% Example

if __name__ == "__main__":
    etapoints = 256
    f = 16.5           # Radar frequency in MHz
    thetastar = 80     #% Direction of the wind i.e. predominant direction of the waves
    s = 4             #% Spreading value (ocean spectra directional dist_n)
    U = 15            # % Wind speed (ms^-1)
    distance = 17    #  % Kms between Tx and Rx (0 for monostatic)
    R = 20           #% Distance from the transmitter to the scatter patch 
    depth = 100
    zeta = 45
    phi = 0        #% Measurement angle (c'clockwise from x axis - degrees 0-180) 
    NQ = 36           # % # of points to go through for each eta value
    current_speed = 4
    current_angle = 45    #    % Measurement angle (c'clockwise from x axis - degrees 0-180) 
    
    t = time.time()
    eta, rcs = vectorised_sim_DS_rotate(f, thetastar, s, U, distance, R, depth, zeta, phi,etapoints, NQ, current_speed, current_angle, 0.00000000001)
    print time.time() - t
    
    eta = eta/(2*np.pi)
    
    plt.figure()
    plt.plot(eta,rcs)
    plt.show()
