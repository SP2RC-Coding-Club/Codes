# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:00:14 2017

@author: Matt
"""

import numpy as np
import slab_functions as sf
import move_seed_points as msp
from pysac.plot.mayavi_seed_streamlines import SeedStreamline
from mayavi import mlab
import mayavi_plotting_functions as mpf
import img2vid as i2v

###############################################################################

# Set the Alfven speed distribution - varying in the x-direction
vA2 = 1.
vA1 = 0.5
def vA_func(x):
    return (vA2 - vA1)/4 * x + vA1


# Which angle shall we view from?
view_options = ['front', 'front-parallel', 'top', 'top-parallel' 'front-top',
                'front-side', 'front-top-side']
#view = 'front'
#view = 'front-parallel'
#view = 'top'
#view = 'top-parallel'
#view = 'front-top'
#view = 'front-side'
view = 'front-top-side'

# Uniform lighting?
#uniform_light = True
uniform_light = False

show_mag = False
show_mag_scale = False
show_mag_fade = False
show_mag_vec = False
show_vel_top = False
show_disp_top = False
show_axes = False
show_axis_labels = False
show_mini_axis = False

# Uncomment the parametrer you would like to see
# No density perturbations or vel/disp pert for alfven modes.
show_mag = True
#show_mag_scale = True #must also have show_mag = True
#show_mag_fade = True
#show_mag_vec = True
#show_vel_top = True
#show_disp_top = True
show_axes = True
#show_axis_labels = True
show_mini_axis = True

# Video resolution
#res = (1920,1080)
res = tuple(101 * np.array((16,9)))
#res = tuple(51 * np.array((16,9)))
#res = tuple(21 * np.array((16,9)))

number_of_frames = 50

fps = 20

#save_images = False
save_images = True

#make_video = False
make_video = True

mode = 'alfven-mixed-driver'

#
##
###
####
#####
######
#######
########
#########

    
print('Starting ' + mode)
    
# Specify oscillation parameters
K = 2.
W = vA1
    
# Dependent variables:
# x = k*x
# y = k*y
# z = k*z
# W = omega/k
# K = k*x_0
# t = omega*t
    
#################################################################################
    
# Grid range
xmin = 0.
xmax = 4.
ymin = 0.
ymax = 4.
zmin = 0.
zmax = 2*np.pi
        
# You can change ny but be careful changing nx, nz.
nx = 100 #100
ny = 100 #100#20 #100
nz = 100 #100
nt = number_of_frames

# t (time) range.
t_start = 0.
t_end = 2*zmax
        
t = t_start

xvals = np.linspace(xmin, xmax, nx)
yvals = np.linspace(ymin, ymax, ny)
zvals = np.linspace(zmin, zmax, nz, endpoint=False)

x_spacing = max(nx, ny, nz) / nx
y_spacing = max(nx, ny, nz) / ny
z_spacing = max(nx, ny, nz) / nz
spacing =  np.array([x_spacing, z_spacing, y_spacing])
        
# For masking points
mod = int(3 * nx / 100)
mod_y = int(np.ceil(mod / y_spacing))

# Displacement data
xixvals_t = np.real(np.repeat(sf.xix_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
xizvals_t = np.real(np.repeat(sf.xiz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
xiyvals_t = np.real(np.repeat(sf.xiy_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))

# Velocity field data
if show_vel_top == True:
    vxvals_t = np.real(np.repeat(sf.vx_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
    vzvals_t = np.real(np.repeat(sf.vz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
    vyvals_t = np.real(np.repeat(sf.vy_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
        
# Magnetic field data
bxvals_t = np.real(np.repeat(sf.bx_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
byvals_t = np.real(np.repeat(sf.by_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
bz_eq3d = np.repeat(sf.bz_eq_amd(xvals, zvals)[:, :, np.newaxis], ny, axis=2)
bzvals_t = np.real(np.repeat(-sf.bz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2) +
                   bz_eq3d)

#Masking points - Have to do this manually due to Mayavi bug
if show_mag_vec == True:
    bxvals_mask_front_t, byvals_mask_front_t, bzvals_mask_front_t = mpf.mask_points(bxvals_t, byvals_t, bzvals_t, 
                                                                       'front', mod, mod_y)
if show_disp_top == True:    
    xixvals_mask_top_t, xiyvals_mask_top_t, xizvals_mask_top_t = mpf.mask_points(xixvals_t, xiyvals_t, xizvals_t, 
                                                                                 'top', mod, mod_y)
if show_vel_top == True:  
    vxvals_mask_top_t, vyvals_mask_top_t, vzvals_mask_top_t = mpf.mask_points(vxvals_t, vyvals_t, vzvals_t, 
                                                                              'top', mod, mod_y)
                                                                                    
# Initialise figure
fig = mlab.figure(size=res) # (1920, 1080) for 1080p , tuple(101 * np.array((16,9))) #16:9 aspect ratio for video upload

xgrid, zgrid, ygrid = np.mgrid[0:nx:(nx)*1j,
                               0:nz:(nz)*1j,
                               0:ny:(ny)*1j]

# Create magnetic field visualisation
field = mlab.pipeline.vector_field(bxvals_t, bzvals_t, byvals_t, name="B field", 
                                   figure=fig, scalars=zgrid)
field.spacing = spacing
field_ms = field.mlab_source

scalefactor = 4. # scale factor for direction field vectors

# Velocity direction field at the top of the domain
if show_vel_top == True:
    vdirfield_top = mlab.pipeline.vector_field(vxvals_mask_top_t, np.zeros_like(vxvals_mask_top_t),
                                                vyvals_mask_top_t, name="V field top",
                                                figure=fig)
    vdirfield_top.spacing = spacing
    mpf.vector_cut_plane(vdirfield_top, 'top', ny, nz, y_spacing, scale_factor=scalefactor)

#Set viewing angle
mpf.view_position(fig, view, nx, ny, nz)

if show_axes == True:
    mpf.axes(field, 'false', view)
        
if show_mini_axis == True:
    oa = mlab.orientation_axes(xlabel='x', ylabel='z', zlabel='y')
    oa.marker.set_viewport(0,0,0.25,0.25) # minx, miny, maxx, maxy

if uniform_light == True:
    #uniform lighting, but if we turn shading of volumes off, we are ok without
    mpf.uniform_lighting(fig)

#Black background
mpf.background_colour(fig, (0., 0., 0.))

for t_ind in range(nt):
    if t_ind != 0:
        # Update datasets at each time step
        xixvals_t = np.real(np.repeat(sf.xix_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
        xizvals_t = np.real(np.repeat(sf.xiz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
        xiyvals_t = np.real(np.repeat(sf.xiy_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
        
        if show_vel_top == True:
            vxvals_t = np.real(np.repeat(sf.vx_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
            vzvals_t = np.real(np.repeat(sf.vz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
            vyvals_t = np.real(np.repeat(sf.vy_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
                
        bxvals_t = np.real(np.repeat(sf.bx_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
        byvals_t = np.real(np.repeat(sf.by_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
        bz_eq3d = np.repeat(sf.bz_eq_amd(xvals, zvals)[:, :, np.newaxis], ny, axis=2)
        bzvals_t = np.real(np.repeat(-sf.bz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2) +
                           bz_eq3d)
        
        #Masking points
        if show_mag_vec == True:                            
            mpf.mask_points(bxvals_t, byvals_t, bzvals_t, 'front', mod, mod_y)
            
        if show_disp_top == True:    
            mpf.mask_points(xixvals_t, xiyvals_t, xizvals_t, 'top', mod, mod_y)
        
        if show_vel_top == True: 
            mpf.mask_points(vxvals_t, vyvals_t, vzvals_t, 'top', mod, mod_y)
        
    #
    ##
    ###
    ####
    #####
    ####
    ###
    ##
    #

    # Update mag field visualisation
    if t_ind != 0:
        field_ms.set(u=bxvals_t, v=bzvals_t, w=byvals_t)
    
    # Create field lines
    if show_mag == True:
        # Create an array of points for which we want mag field seeds
        nx_seed = 9
        ny_seed = 13
        start_x = 2.
        end_x = nx+1 - start_x
        start_y = 1.
        if ny == 20:
            end_y = ny - 1
        elif ny == 100:
            end_y = ny - 2
        else:
            end_y = ny - 1
        seeds=[]
        dx_res = (end_x - start_x) / (nx_seed-1)
        dy_res = (end_y - start_y) / (ny_seed-1)
        for j in range(0,ny_seed):
            for i in range(0,nx_seed):
                x = start_x + (i * dx_res) * x_spacing
                y = start_y + (j * dy_res) * y_spacing
                z = 1.# + (t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
                seeds.append((x,z,y))
        
        for i in range(nx_seed):
            del seeds[0]
            del seeds[-1]
        
        # I have managed to get the move_seeds function to work here.
        # Seed points are moved by the velocity field.
        seeds = msp.move_seeds(seeds, xixvals_t, xiyvals_t, xizvals_t, 
                               [xmin, ymin, zmin], [xmax, ymax, zmax])
        if t_ind != 0:
            field_lines.remove()
        field_lines = SeedStreamline(seed_points=seeds)
        field_lines.stream_tracer.integration_direction='both'
        field_lines.streamline_type = 'tube'        
        field_lines.stream_tracer.maximum_propagation = 500.
        field_lines.tube_filter.number_of_sides = 20
        field_lines.tube_filter.radius = 0.7
        field_lines.tube_filter.capping = True
        field_lines.actor.property.opacity = 1.0
                
        field.add_child(field_lines)
        module_manager = field_lines.parent
        
        # Colormap of magnetic field strength plotted on the field lines
        if show_mag_scale == True:
            module_manager.scalar_lut_manager.lut_mode = 'coolwarm'
            module_manager.scalar_lut_manager.data_range=[7,18]
        else:
            mag_lut = module_manager.scalar_lut_manager.lut.table.to_array()
            mag_lut[:,0] = [220]*256
            mag_lut[:,1] = [20]*256
            mag_lut[:,2] = [20]*256
            module_manager.scalar_lut_manager.lut.table = mag_lut
        if show_mag_fade == True:
            mpf.colormap_fade(module_manager, fade_value=20)
            
        # Direction field of Magntic field vectors
        if show_mag_vec == True:
            bdirfield_front = mlab.pipeline.vector_field(bxvals_mask_front_t, bzvals_mask_front_t,
                                                         byvals_mask_front_t, name="B field front",
                                                         figure=fig)
            bdirfield_front.spacing = spacing
            mpf.vector_cut_plane(bdirfield_front, 'front', ny, nz, y_spacing, scale_factor=scalefactor)
            
            
    if show_vel_top == True:
        vxvals_mask_top_t, vyvals_mask_top_t, vzvals_mask_top_t = mpf.mask_points(vxvals_t, vyvals_t, vzvals_t, 
                                                                                                  'top', mod, mod_y)                                                
        vdirfield_top.mlab_source.set(u=vxvals_mask_top_t, v=vzvals_mask_top_t, w=vyvals_mask_top_t)
        
    if show_disp_top == True:
        xidirfield_top = mlab.pipeline.vector_field(xixvals_mask_top_t, np.zeros_like(xixvals_mask_top_t),
                                                    xiyvals_mask_top_t, name="Xi field top",
                                                    figure=fig)
        xidirfield_top.spacing = spacing
        mpf.vector_cut_plane(xidirfield_top, 'top', ny, nz, y_spacing, scale_factor=scalefactor)

    #Set viewing angle
    mpf.view_position(fig, view, nx, ny, nz)

    if save_images == True:
        prefix = 'bouncing_field_test_amd_' + mode + '_' + view + '_vel_top_bounce_test'
        mlab.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_animations\\'
                 + prefix + str(t_ind+1) + '.png')
        
    t = t + (t_end - t_start) / nt
    
if make_video == True:
    mlab.close(fig)
    i2v.image2video(prefix=prefix, output_name=prefix + '_video', 
                    out_extension='mp4', fps=fps, n_loops=1, delete_images=True,
                    delete_old_videos=True, cover_page=False, res=res[1])
print('Finished ' + mode)
    