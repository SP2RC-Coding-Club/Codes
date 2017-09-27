#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 15:06:59 2017

To find the wind direction and spreading values amongst noisy data 

Use fuzzy ratios to find the value for beta and r1 which minimises the change in the surrounding 
ratio values

@author: rach
"""

#%% imports 
from __future__ import division
import numpy as np
import os 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cPickle as pickle
import datetime 
import sys
cwd = os.getcwd()
sys.path.insert(0, cwd + "/../../../Programming/Read_data/")
sys.path.insert(0, cwd + "/../../../Programming/Geometry_Functions/")
sys.path.insert(0, cwd + "/../../../Programming/Beamforming/")
sys.path.insert(0, cwd + "/../../../Programming/Extra_functions/")
from read_data_functions import Winds_from_ASCAT, Winds_from_buoy, geo_to_polar 
from functions_for_winds import alphas_from_betas, return_chosen_winds_spiral, lat_plus_metres_to_lat, long_plus_metres_to_long, get_min_beta 
import time
from functions_for_fuzzy import get_ratios_from_r1_and_beta
from scipy.optimize import minimize

#%% cost function

def cost_function_fuzzy_GA_nm(individual, ratios, phis, r1_index_x,r1_index_y, THRESHOLD=4):
    """
    Inputs:
        individual - [r1, beta] the values to be assessed
        ratios - the original array of measured ratios
        phis - the original array of bistatic angles
        THRESHOLD - the amount that renders a valid output *this may be obsolete depending on choices 
        r1_index_x,r1_index_y - the index of the desired point to assess r1/beta choice 
    """
    
    np.reshape(np.asarray(individual),(2,1))

    try:
        R1, beta = individual[0][0], individual[0][1]
    except IndexError:
          R1, beta = individual[0], individual[1]
     
    min_beta = get_min_beta(R1)
    if beta< min_beta:
        return np.abs(beta-min_beta)*10000
    
    # get the original values - just for reference at this point
    orig_phi = phis[r1_index_x, r1_index_y]
    
    J = np.reshape(np.array([0.0,0.0]),(1,2)) # initialise the cost counter
    counter = np.reshape(np.array([0,0]),(1,2)) # initialise the counter for keeping track of the costs
    
#    get the ratios and phis around the desired point
    ratios_pt = ratios[r1_index_x-1:r1_index_x+2,r1_index_y-1:r1_index_y+2]
    phis_pt = phis[r1_index_x-1:r1_index_x+2,r1_index_y-1:r1_index_y+2]
       
    idxs = np.array([[0,1],[1,2],[2,1],[1,0]]) # the indices of the nearest 4 (easier to control this way)
        
#    select the 4 values nearest to the desired point - both phis and ratios 
    phis_neighbours1 = np.array([phis_pt[idxs[0,0],idxs[0,1]],phis_pt[idxs[1,0],idxs[1,1]],phis_pt[idxs[2,0],idxs[2,1]],phis_pt[idxs[3,0],idxs[3,1]]])
    ratios_neighbours1 = np.array([ratios_pt[idxs[0,0],idxs[0,1]],ratios_pt[idxs[1,0],idxs[1,1]],ratios_pt[idxs[2,0],idxs[2,1]],ratios_pt[idxs[3,0],idxs[3,1]]])
    
#    loop over the neighbour values in the order specified in idxs
    for i, phi1 in enumerate(phis_neighbours1): 
        orig_rat_1 = ratios_neighbours1[i] # original value of ratio to calc cost
        # for the centre point, specified by R1 and beta, find the values the neighbour
        # value must take to pattern fit correctly given their phi values        
        [ratios_picked_1, winds_picked_1] = get_ratios_from_r1_and_beta(R1, beta, orig_phi, phi1)
        
        if i==0:
            winds_assert_shape=np.reshape(np.asarray(winds_picked_1[:]),(1,-1))
            
            if winds_assert_shape.size==1:
                wind_choices = np.empty_like(J)
                wind_choices[0,0] = np.mod(winds_assert_shape,2*np.pi)
            else:
                wind_choices = np.mod(winds_assert_shape,2*np.pi)
                        
        # get a new set of neighbour data to check the fit of the parameters further 
        # new centre point index (row/column = x/y)
        [neighbour_idx_x, neighbour_idx_y] = [r1_index_x,r1_index_y] + (idxs[i] - 1)
        
        # get the 9grid first
        ratios_neighbour_new = ratios[neighbour_idx_x-1:neighbour_idx_x+2, neighbour_idx_y-1:neighbour_idx_y+2]
        phis_neighbour_new = phis[neighbour_idx_x-1:neighbour_idx_x+2, neighbour_idx_y-1:neighbour_idx_y+2]
        # select the 4 nearest in the 9grid
        
        phis_neighbours2 = np.array([phis_neighbour_new[idxs[0,0],idxs[0,1]], phis_neighbour_new[idxs[1,0],idxs[1,1]], phis_neighbour_new[idxs[2,0],idxs[2,1]], phis_neighbour_new[idxs[3,0],idxs[3,1]]])
        ratios_neighbours_2 = np.array([ratios_neighbour_new[idxs[0,0],idxs[0,1]], ratios_neighbour_new[idxs[1,0],idxs[1,1]], ratios_neighbour_new[idxs[2,0],idxs[2,1]], ratios_neighbour_new[idxs[3,0],idxs[3,1]]])
             
        # loop over the (two) ratios picked from the pf algorithm 
        for i_rat, R_i in enumerate(ratios_picked_1):
            wind_i = np.mod(winds_picked_1[i_rat],2*np.pi) # corresponding wind direction to the chosen ratio 

            J[wind_choices==wind_i] += 2*np.abs(orig_rat_1-R_i) # calculate cost of the neighbour
            counter[wind_choices==wind_i]+=2 # keep track of the counter 

            # loop through the neighbour phis and get the matching ratio values
            for i_2, phi2_2 in enumerate(phis_neighbours2): 
                orig_rat_2 = ratios_neighbours_2[i_2] # original ratio value for comparing             
                [ratios_picked_2, winds_picked_2] = get_ratios_from_r1_and_beta(R_i, beta, phi1, phi2_2)
                try:
                    ratio_val_for_wind = np.asarray(ratios_picked_2)[np.isclose(np.mod(winds_picked_1,2*np.pi), np.mod(wind_i,2*np.pi), rtol=1e-05, atol=1e-04)]
                    J[wind_choices==wind_i]+= np.abs(orig_rat_2-ratio_val_for_wind[0])
                    counter[wind_choices==wind_i]+=1
                except:
#                    print 'random error 2'
                    
                    return 1000
                    
    idx_min = np.nanargmin(J)
    
    if (counter==0).all():
        J_return = 1000
    else:
        J_return = J[0,idx_min]
    return J_return

#%% get data 
# range of loaded data
RANGES = range(3,16)
BEAMS = range(120)

# input file names
#file_name_radar = '20140712.SORT'
file_name_buoy = 'CANDHIS_export_pem_08302_Base.csv'
file_name_ASCAT = 'OASWC12_20140712_084200_09415_M01.nc'

# get the data directories 
dir_link = os.getcwd() # get the name of the directory this program is in
data_dir = dir_link + "/../Data/" # the data should be in a subfolder called 'Data' in the folder the program's run from 
data_dir_radar = data_dir + 'Radar/'
data_dir_buoy = data_dir + 'French_buoy_data/'
data_dir_ASCAT = data_dir + 'EUMETSAT/JULY_2014/'

#%% get radar data
Rs_r = pickle.load(open('Rs_r.pickle','rb'))
Rs_ps = pickle.load(open('Rs_ps.pickle','rb'))
Rs_mv = pickle.load(open('Rs_mv.pickle','rb'))
xs = pickle.load(open('xs.pickle','rb'))
ys = pickle.load(open('ys.pickle','rb'))
thetas = pickle.load(open('thetas.pickle','rb'))
phibis = pickle.load(open('phibis.pickle','rb'))
braggs = pickle.load(open('braggs.pickle','rb'))
omegas = pickle.load(open('omegas.pickle','rb'))
time_measurement = pickle.load(open('time_measure.pickle','rb'))
                      
print 'radar data read'
print 'time: ' + str(time_measurement)

#%% get buoy winds
[buoy_time, buoy_dir, buoy_dir_mean, buoy_hm0] = Winds_from_buoy(file_name_buoy, data_dir_buoy, time_measurement)
buoy_dir+=180
buoy_dir_polar = geo_to_polar(buoy_dir) # 
print 'buoy data read'

# %% get ASCAT winds 
[lons_ASCAT, lats_ASCAT, times_ASCAT, wind_dir_ASCAT, wind_dir_model_ASCAT] = Winds_from_ASCAT(file_name_ASCAT, data_dir_ASCAT)
wind_dir_ASCAT_np = np.array(wind_dir_ASCAT)
plot_ASCAT_dirs = geo_to_polar(wind_dir_ASCAT) # 
print 'ascat data read'

# %% get french map data and et locations of buoy/rx/tx

#limits for the quiver plot (lat/lon)
x1,x2 = 6,6.5 #longitudes
y1,y2 = 42.7,43.1 #latitudes

#get the basemap data for the french coast 

#france_map = Basemap(resolution='f',projection='merc',area_thresh = 0.1, llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_ts=(x1+x2)/2)
#c (crude), l (low), i (intermediate), h (high), f (full)
#pickle.dump(france_map, open('fr_map.pickle','wb'),-1) # pickle for speed 
france_map = pickle.load(open('fr_map.pickle','rb'))

#tx and rx coordinates and transformation into map coords 
lonloc_Tx,latloc_Tx = 6.204189,42.983084
lon_Tx, lat_Tx = france_map(lonloc_Tx,latloc_Tx)

lonloc_Rx,latloc_Rx = 6.357316666,43.091933333
lon_rx, lat_rx = france_map(lonloc_Rx,latloc_Rx)
# wave buoy data 
buoy_lon,buoy_lat = 6.280556, 42.96667
lon_buoy, lat_buoy = france_map(buoy_lon,buoy_lat)

x_ascat,y_ascat = france_map(lons_ASCAT,lats_ASCAT) # ascat data -> french map coords
# get the latitude and longitudes of the scatter points from the radar data set
lat = [lat_plus_metres_to_lat(latloc_Tx, change) for change in ys]
lon = [long_plus_metres_to_long(lonloc_Tx, change, latloc_Rx) for change in xs ]
# convert the radar data points into the french map coords 
x,y = france_map(lon,lat)
# Note that this is from the TX because the scatter geometry functions returns the sp from a coordinate system with Tx at (0,0)

# %% time check of the ASCAT data - get the time of the data point nearest to the radar coverage area 

target_lat = 43
target_lon = 6.35
diff = np.inf # set as upper limit for initial comparison 
for i in range(wind_dir_ASCAT.shape[0]):
    for j in range(wind_dir_ASCAT.shape[1]):
        lat = lats_ASCAT[i,j]
        lon = lons_ASCAT[i,j]
        lat_d = target_lat - lat
        lon_d = target_lon - lon
        diff_n = np.sqrt(lat_d**2 + lon_d**2)
        if diff_n < diff:
            diff = diff_n
            time_ix = i
            space_ix = j

time_measurement = times_ASCAT[time_ix,space_ix]
julian_time = datetime.datetime(1990,1,1,0,0,0)
time_of_measure = julian_time + datetime.timedelta(seconds=int(time_measurement))

#%% choose radar data and mask ratios where necessary - also where smoothing may take place 

# reshape arrays for the rest of the program (choosing certain sections) 
ratios_r = np.reshape(Rs_ps, (len(RANGES),len(BEAMS)))[:,20:100:5] 
phis = np.mod(np.reshape(braggs, (len(RANGES),len(BEAMS))), 2*np.pi)[:,20:100:5]
XS = np.reshape(x, (len(RANGES),len(BEAMS)))[:,20:100:5]
YS = np.reshape(y, (len(RANGES),len(BEAMS)))[:,20:100:5]

R = np.ma.masked_where(np.isnan(ratios_r),ratios_r) # mask the nans of R 

#%% Minimise the cost function to fit each beta and ratio value 
print 'pattern fitting on fuzzy ratios about to begin'

THRESHOLD = 4
MIN_BETA,MAX_BETA = 0.05,10

BETA_RANGE = np.arange(0,10,.1)

NO_X,NO_Y = 9,12
lon_points,lat_points,phi_picked = np.empty((NO_X,NO_Y)), np.empty((NO_X,NO_Y)), np.empty((NO_X,NO_Y))
alphas_picked_nm, beta_points_nm = np.empty((NO_X,NO_Y)),np.empty((NO_X,NO_Y))

x0 = [[0.0, 2]]
# loop through the specified range to get beta and ratio values and hence wind direction 
counter_r = 0
for range_idx in np.arange(NO_X)+2: 
    counter_b = 0
    for beam_idx in np.arange(NO_Y)+2:
        # store position values for plotting on the map 
        lon_points[counter_r,counter_b] = XS[range_idx,beam_idx] 
        lat_points[counter_r,counter_b] = YS[range_idx,beam_idx]
        phi_picked[counter_r,counter_b] = phis[range_idx,beam_idx] # bragg angle
        print [range_idx, beam_idx]
        t_nm = time.time()
        # find the minimum of the cost function using Nelder Mead 
        x0[0][0] = R[range_idx,beam_idx] # an initial guess for the minimising algorithm 
        res = minimize(cost_function_fuzzy_GA_nm, x0[0], args=(R, phis, range_idx,beam_idx), method='nelder-mead', options={'xtol': 1e-8, 'disp': False})
        print 'Nelder Mead time: '+str(time.time()-t_nm)+' secs'
        
        nm_r1,nm_beta,j_nm = res.x[0], res.x[1], res.fun # solutions and the cost (NB the cost is for plotting)
        print j_nm
        
        nan
        
        counter_b+=1
    counter_r+=1

# pick the most consistent winds from the 2 options to plot 
chosen_winds_nm = return_chosen_winds_spiral(phi_picked+alphas_picked_nm, phi_picked-alphas_picked_nm, beta_points_nm.shape[0], beta_points_nm.shape[1], [0,0])[0]     

# %% plot the final wind vectors on the french map 

fig = plt.figure()
france_map.drawcountries(linewidth=0.5)
france_map.drawcoastlines(linewidth=0.5)
france_map.drawmapboundary(fill_color='lightcyan')
france_map.fillcontinents(color='mediumseagreen',lake_color='skyblue')
france_map.plot(lon_rx, lat_rx, marker='D',color='m')
france_map.plot(lon_buoy, lat_buoy, marker='o',color='g')
france_map.quiver(lon_buoy, lat_buoy, np.cos(np.deg2rad(chosen_winds_nm)), np.sin(np.deg2rad(chosen_winds_nm)),scale = 10, edgecolor='b', facecolor='black')
france_map.plot(lon_Tx, lat_Tx, marker='D',color='b')
txt = plt.title('Beta Background and Wind Vectors: Nelder Mead')
france_map.plot(xs,ys,marker='o',ls='none')
        
betas_masked_nm = np.ma.masked_where(np.isnan(beta_points_nm),beta_points_nm) # mask beta for plotting
france_map.pcolormesh(lon_points, lat_points, betas_masked_nm, cmap='hot', vmin=0, vmax=10)
x,y = france_map(lons_ASCAT,lats_ASCAT)
france_map.quiver(x,y,np.cos(np.deg2rad(plot_ASCAT_dirs)),np.sin(np.deg2rad(plot_ASCAT_dirs)),edgecolor='g', scale=20, facecolor='red', linewidth=.5)

france_map.quiver(lon_points, lat_points, np.cos(chosen_winds_nm), np.sin(chosen_winds_nm), scale=20, facecolor='blue', linewidth=.5)
france_map.quiver(lon_points, lat_points, np.cos(chosen_winds_nm), np.sin(chosen_winds_nm), scale=20, facecolor='blue', linewidth=.5)

plt.colorbar()
plt.show()
