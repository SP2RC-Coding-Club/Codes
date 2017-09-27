# Dealine 5/01/16  object: J005908.99-710648.6

import math
import numpy as np
import pylab as plt
import random
from numpy.random import normal
from scipy.stats import norm
import matplotlib.mlab as mlab

# ----------------------------------------------------------
# Data
# ----------------------------------------------------------
# Loads data from model SED
(lambda_model_nm,flux_model_CGS)= np.loadtxt('data_oats_py',usecols=(2,4),comments='#',unpack=True)
lambda_model_SI = lambda_model_nm*10**(-9) # convert nm to m
flux_model_SI = flux_model_CGS*10**(-3) #
lambda_model_micron = lambda_model_nm/1000
# Speed of light
c = 3*10**(8)# m/s
flux_lambda_model_SI = 4*c*flux_model_SI/(lambda_model_SI**2) #***# !!!!!!!!!!!!!!!!!!!!!!!!!!

# Collected data for catalogs 
# H band for kevin (4 or 5) chi^2 = 11-12, D = 1 N = band + 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(magnitude,magnitude_error,flux_zero_point,Lambda)= np.loadtxt('SED_data.dat',usecols=(0,1,2,3),comments='#',unpack=True)
#(magnitude,magnitude_error,flux_zero_point,Lambda)= np.loadtxt('\alex data\Kia.dat',usecols=(0,1,2,3),comments='#',unpack=True)
Lambda_SI = Lambda*10**(-6) # Convert microns to meters
lambda_nm = Lambda_SI*10**(9) # Convert m to nm
lambda_mircons = lambda_nm/1000 
flux_zero_point_SI = flux_zero_point*10**(-26) # convert Jy to W/(m^2 Hz)
flux_lambda_SI = flux_zero_point_SI*(c/Lambda_SI**2)*(10**(-0.4*magnitude)) # Watt/(m^2 m)
F_lambda_error_SI = 0.4*np.log(10)*magnitude_error*flux_lambda_SI
print flux_lambda_SI
# Loads data from redding law of the galaxy
(wl_Rv,Rv)= np.loadtxt('Rv.dat',usecols=(0,1),comments='#',unpack=True)
wl_Rv_SI = wl_Rv*1e-10 # A to m

index_bands = []
index_model = []
# Gets the index for each band and for the model
for i in range(len(Lambda_SI)):
    index_bands.append(np.argmin(np.abs(wl_Rv_SI - Lambda_SI[i])))
    index_model.append(np.argmin(np.abs(lambda_model_SI - Lambda_SI[i]))) 


# To chose which band to scale to. Band are in order of smallest CWL to largest
band = 6
band_name = ['U','B','V','I','J','H','K','w1','3.6']
print 'Scaled to', band_name[band], 'band.'
#%-------------------------------------------------------------------------------
#U_index = np.argmin(np.abs(wl_Rv_SI - Lambda_SI[0])) # band = 0
#B_index = np.argmin(np.abs(wl_Rv_SI - Lambda_SI[1])) # band = 1
#V_index = np.argmin(np.abs(wl_Rv_SI - Lambda_SI[2])) # band = 2
#I_index = np.argmin(np.abs(wl_Rv_SI - Lambda_SI[3])) # band = 3
#J_index = np.argmin(np.abs(wl_Rv_SI - Lambda_SI[4])) # band = 4
#H_index = np.argmin(np.abs(wl_Rv_SI - Lambda_SI[5])) # band = 5
#K_index = np.argmin(np.abs(wl_Rv_SI - Lambda_SI[6])) # band = 6
#w1_index =np.argmin(np.abs( wl_Rv_SI - Lambda_SI[7])) # band = 7
#three_six_index = np.argmin(np.abs(wl_Rv_SI - Lambda_SI[8])) # band = 8
#%------------------------------------------------------------------------------

#%------------------------------------------------------------------------------
# Chi Squared test
# 1/(n-d) sum ((O_i - C_i)/F_lambda_error_i)^2
# Index gives the position in the array of the closet value to wave length of selected band in the model
index = np.argmin(np.abs(lambda_model_SI-Lambda_SI[band]))
E_B_V_dist = []
min_chi_squared_dist = []
F_lambda_error_SI_random = []
# Set number of iteration for find the E(B-V) value
n = 1000
# Create an array for E(B-V)
E_B_V = np.linspace(0.0,1.0,n)
E_B_V = np.asarray(E_B_V)

N = band - 1.0 # Number of data points in Chi^2
D = 2.0 # number of dependent variables
#for x in range(10000):
for x in range(2):
    print 'run', x
    final_chi_squared = []
    F_lambda_error_SI_random = random.gauss(0.0,F_lambda_error_SI) 
    for E in np.linspace(0,1.0,n):
        chi_list = []
        chi_sum_1 = []
        C = []  
        scaled_flux_lambda_model_SI_C = []
        O = []
        Deredden_flux_lambda_SI = []
        Deredden_lambda_flux_lambda_SI = []
        chi = [] 
        for j in range(len(Lambda_SI)):    
#            Monte carlo
            Deredden_flux_lambda_SI.append(((flux_lambda_SI[j]+F_lambda_error_SI_random[j])*10**(E*0.4*Rv[index_bands[j]])))
#            Derreden flux without random error
#            Deredden_flux_lambda_SI.append(((flux_lambda_SI[j])*10**(E*0.4*Rv[index_bands[j]])))
            Deredden_lambda_flux_lambda_SI.append((Lambda_SI[j]*(flux_lambda_SI[j]+F_lambda_error_SI_random[j])*10**(E*0.4*Rv[index_bands[j]])))
            # To create an array once all photometric data has been dereddened        
            if j == len(Lambda_SI)-1:
                O = Deredden_flux_lambda_SI
                Deredden_lambda_flux_lambda_SI_O = Deredden_lambda_flux_lambda_SI
                # Scales observed data to certain band fromt he derredened model
                scaled_factor = O[band]/flux_lambda_model_SI[index]
                scaled_flux_lambda_model_SI = (flux_lambda_model_SI*scaled_factor) # scaled flux of the model 
                for i in range(band+1): 
                    scaled_flux_lambda_model_SI_C.append(scaled_flux_lambda_model_SI[index_model[i]])
                    # To create an array once the  all desired fluxes have been scaled               
                    if i == band:
                        C = scaled_flux_lambda_model_SI_C
#                        #To check the scaling is correct                   
#                        fig = plt.figure(figsize=(6,4), dpi= 120)
#                        ax =fig.add_subplot(111)
#                        fig.tight_layout()
#
#                        plt.plot((lambda_model_micron),(lambda_model_SI*scaled_flux_lambda_model_SI),lw='2', color= 'blue')
#                        plt.scatter((lambda_mircons),(Deredden_lambda_flux_lambda_SI_O), color= 'black')
#                        plt.scatter((lambda_mircons),((Lambda_SI*flux_lambda_SI)), color= 'red')
#                        ax.set_xscale('log')
#                        ax.set_yscale('log')
#                        ax.set_ylabel(r'$\lambda f(\lambda)$ (Watt/m^2)')
#                        ax.set_xlabel(r'$\lambda$ (micron)')
#                        plt.xlim((10**(-1),10**(2)))
#                        plt.ylim((10**(-16),10**(-13)))                                        
#                        plt.show()
        for k in range(band+1):
            # This calculates the chi^2
            Chi = (((N-D))**(-1))*((O[k]-C[k])/F_lambda_error_SI[k])**2
            chi.append(Chi)
            # Once chi array is to the length of the point before band scaled to, the array needs to be summed to get a chi^2 value
            if  k == band:      
                chi_squared = sum(chi)
                final_chi_squared.append(chi_squared)
                # Gives minimuim chi^2 values
                min_chi_squared =  min(final_chi_squared)               
                # Finds the min E(B-V) value in the Chi^2
                min_E_B_V = np.argmin(final_chi_squared)
                # This gives the E(B-V) for our model
                E_B_V_val = E_B_V[min_E_B_V]
                if E == 1:
                    # When E(B-V) reaches its last value it will save as below arrays allowing for plot distrbution of E(B-V) from the monte carlo simulation
                    min_chi_squared_dist.append(min_chi_squared)
                    E_B_V_dist.append(E_B_V_val)
print 'min chi^2 val', min_chi_squared
print 'E(B-V) value is:', E_B_V_val

## To plot Histogram
##-------------------------------------------------------------------------------------------
#fig = plt.figure(figsize=(6,4), dpi= 120)
#ax =fig.add_subplot(111)
#fig.tight_layout()
#    
#mean = np.mean(E_B_V_dist)
#std = np.std(E_B_V_dist)      
#x = np.linspace(0.4,.42,1000)
#plt.hist(E_B_V_dist, normed = True)
#plt.plot(x,mlab.normpdf(x,mean,std), linewidth= 3)
#plt.show
#-------------------------------------------------------------------------------------------------------
# Plot Chi Squared dist 
fig = plt.figure(figsize=(6,4), dpi= 120)
ax  = fig.add_subplot(111)
fig.tight_layout()
plt.plot(E_B_V,np.log10(final_chi_squared))
plt.xlabel('$E(B-V)$')
plt.ylabel('$log(\chi^2)$')


# ----------------------------------------------------------------------------------------------------------

# To derreden SED with the Calculated E(B-V).
deredden_Lambda_flux_lambda_SI = []
deredden_scaler_flux_lambda_SI = []
for i in range(len(Lambda_SI)):
    deredden_scaler_flux_lambda_SI.append(flux_lambda_SI[i]*10**(E_B_V_val*0.4*Rv[index_bands[i]]))    
    deredden_Lambda_flux_lambda_SI.append(Lambda_SI[i]*flux_lambda_SI[i]*10**(E_B_V_val*0.4*Rv[index_bands[i]]))

best_scaled_factor = deredden_scaler_flux_lambda_SI[band]/flux_lambda_model_SI[index]
best_scaled_flux_lambda_model_SI = (flux_lambda_model_SI*best_scaled_factor)

#------------------------------------------------------------------------------------------------------      
# Plots SED      
fig = plt.figure(figsize=(6,4), dpi= 120)
ax  = fig.add_subplot(111)
fig.tight_layout()


plt.plot((lambda_model_micron),(lambda_model_SI*best_scaled_flux_lambda_model_SI),lw='2', color= 'blue', label='SED Model') # Model SED
plt.scatter((lambda_mircons),(Lambda_SI*flux_lambda_SI), color= 'red', label='Photometric Fluxes') # Spectra
plt.scatter((lambda_mircons),(deredden_Lambda_flux_lambda_SI), color= 'black', label='Photometric Fluxes Deredden') # Deredden spectra
plt.scatter((lambda_mircons[band]),(deredden_Lambda_flux_lambda_SI[band]), color= 'yellow', label='Photometric Fluxes Deredden') # Highlight scaled band
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r'$\lambda f(\lambda)$ (Watt/m^2)')
ax.set_xlabel(r'$\lambda$ (micron)') # fix lable after !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# micron range
plt.ylim((10**(-16),10**(-13)))
plt.xlim((10**(-1),10**(2)))

#------------------------------------------------------------------------------------------
# Calculation of basic stellar parameters
lambda_model_SI_int = []
best_scaled_flux_lambda_model_SI_int = []
for i in range(index_model[band]+1):
	lambda_model_SI_int.append(lambda_model_SI[i])
	best_scaled_flux_lambda_model_SI_int.append(best_scaled_flux_lambda_model_SI[i])
 # Intergrates under curves from selected points using trapsuim rule
flux_star = np.trapz(best_scaled_flux_lambda_model_SI_int,lambda_model_SI_int)
print 'flux_star', flux_star

Lambda_SI_int_L_TE = []
best_scaled_flux_lambda_SI_int_L_TE = []
for i in range(band+1,len(Lambda_SI)):
	Lambda_SI_int_L_TE.append(Lambda_SI[i])
	best_scaled_flux_lambda_SI_int_L_TE.append(deredden_scaler_flux_lambda_SI[i])
flux_L_TE = np.trapz(best_scaled_flux_lambda_SI_int_L_TE,Lambda_SI_int_L_TE)
print 'flux_L_TE', flux_L_TE



# distance from Earth to SMC
D_SMC_SI = (2*10**(5)*(9.4605284*10**15))# Kpc to pc to m
print 'D_SMC_SI', D_SMC_SI
T_eff = 5500 # K
L_solar = 3.846*10**26  #W
boltzmann_constant = 5.670*10**(-8)  # W m^−2 K^−4

# Calculation of lumnosity
# L = sigma(T_eff^4)pi(R^2)
luminosity_star = (4*np.pi*flux_star*(D_SMC_SI**2)) # W
luminosity_star_solar = (4*np.pi*flux_star*(D_SMC_SI**2))*((L_solar)**(-1)) # L_sun
print 'luminosity of star in soalr units is:'
print luminosity_star_solar
# R = sqrt(L /sigma(T_eff^4)pi)
#value for radius
Radius_SI = ((luminosity_star)/(4*np.pi*boltzmann_constant*(T_eff**4)))**(0.5) # m
Radius_SI_km = Radius_SI/1000 # km
Radius_SI_mega_km = Radius_SI_km*(10**(-6)) # mega km
Raduis_solar = Radius_SI_km*(695500**(-1))
print 'Radius in mega km is:'
print Radius_SI_mega_km
print 'Radius in solar is:'
print Raduis_solar

# Calculation of dust tmepature
# Dust tempature = (2.9/lambda_max) [mm][K]/[mm] from Weins Law
lambda_peak = lambda_nm[11]/(1000000) # mm
T_d = 2.9/lambda_peak # K
print 'tempature of dust in Kelvin is:'
print T_d
# Distance of the inner shell from the star
# D_shell = (0.5)*R_star*(T_eff/T_d)^2
D_shell = (0.5)*Radius_SI*((T_eff/T_d)**2)
print 'The distance from of the inner shell from the star in meters is:' 
print D_shell
# Infrared Lumnosity 
# L_TE = 4*pi*D_smc^2*f_measured_TE
L_TE = 4*np.pi*(D_SMC_SI**2)*flux_L_TE
print 'lumnnosity of thermal emmission in W:'
print L_TE  