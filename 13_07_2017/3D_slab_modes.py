
#import pdb # pause code for debugging at pdb.set_trace()
import numpy as np
import toolbox as tool
import slab_functions as sf
from pysac.plot.mayavi_seed_streamlines import SeedStreamline
import matplotlib.pyplot as plt
from mayavi import mlab
import gc
#import move_seed_points as msp
import mayavi_plotting_functions as mpf
import dispersion_diagram
import img2vid as i2v

###############################################################################

# What mode do you want? OPTIONS:
mode_options = ['slow-kink-surf', 'slow-saus-surf', 'slow-saus-body-3',
                'slow-kink-body-3', 'slow-saus-body-2', 'slow-kink-body-2', 
                'slow-saus-body-1', 'slow-kink-body-1', 'fast-saus-body-1',
                'fast-kink-body-1', 'fast-saus-body-2', 'fast-kink-body-2',
                'fast-saus-body-3', 'fast-kink-body-3', 'fast-kink-surf',
                'fast-saus-surf', 'shear-alfven', 'shear-alfven-broadband']

# Create sublists of modes
alfven_mode_options = []
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
    if 'alfven' in mode:
        alfven_mode_options.append(mode)
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

# Which angle shall we view from? OPTIONS:
view_options = ['front', 'front-parallel', 'top', 'top-parallel', 'front-top',
                'front-side', 'front-top-side']

# Uniform lighting?
#uniform_light = True
uniform_light = False

show_density = False
show_density_pert = False
show_mag = False
show_mag_scale = False
show_mag_fade = False
show_mag_vec = False
show_vel_front = False
show_vel_front_pert = False
show_vel_top = False
show_vel_top_pert = False
show_disp_top = False
show_disp_front = False
show_axes = False
show_axis_labels = False
show_mini_axis = False
show_boundary = False

# Uncomment the parametrer you would like to see
# No density perturbations or vel/disp pert for alfven modes.

#show_density = True
#show_density_pert = True
show_mag = True
#show_mag_scale = True #must also have show_mag = True
#show_mag_fade = True
#show_mag_vec = True
#show_vel_front = True
#show_vel_front_pert = True
#show_vel_top = True
show_vel_top_pert = True
#show_disp_top = True
#show_disp_front = True
show_axes = True
#show_axis_labels = True
show_mini_axis = True
show_boundary = True

# Visualisation modules in string form for file-names
vis_modules = [show_density, show_density_pert, show_mag, show_mag_scale, 
               show_mag_fade, show_mag_vec, show_vel_front, show_vel_front_pert,
               show_vel_top, show_vel_top_pert, show_disp_top, show_disp_front]
vis_modules_strings = ['show_density', 'show_density_pert', 'show_mag', 'show_mag_scale', 
                       'show_mag_fade', 'show_mag_vec', 'show_vel_front', 'show_vel_front_pert',
                       'show_vel_top', 'show_vel_top_pert', 'show_disp_top', 'show_disp_front']
vis_mod_string = ''
for i in range(len(vis_modules)):
    if vis_modules[i] == True:
        vis_mod_string = vis_mod_string + vis_modules_strings[i][5:] + '_'


# Set to True if you would like the dispersion diagram with chosen mode highlighted.
show_dispersion = False
#show_dispersion = True

# Wanna see the animation? Of course you do
#show_animation = False
show_animation = True

# Basic plot to see which eigensolutions have been found.
show_quick_plot = False
#show_quick_plot = True

# Video resolution
#res = (1920,1080) # There is a problem with this resolution- height must be odd number - Mayavi bug apparently
res = tuple(101 * np.array((16,9)))
#res = tuple(51 * np.array((16,9)))
#res = tuple(21 * np.array((16,9)))

number_of_frames = 50

# Frames per second of output video
fps = 20

#save_images = False
save_images = True

#make_video = False
make_video = True

#
##
###
####
#####
######
#######
########
#########

# Variable definitions:
# x = k*x
# y = k*y
# z = k*z
# W = omega/k
# K = k*x_0
# t = omega*t

if np.array([show_density, show_vel_front_pert, show_vel_top_pert]).any() == True and mode in alfven_mode_options:
    raise NameError('Cannot show density or vel/disp pert for this mode')

# Loop through selected modes
for mode_ind in [0,1]:#range(8,14): # for all others. REMEMBER SBB pparameters
#for mode_ind in [14,15]: #for fast body surf. REMEMBER SBS parameters
#for mode_ind in [16, 17]:
#for mode_ind in [13]: #for an individual mode
#for mode_ind in range(2,14): 
    if mode_ind not in range(len(mode_options)):
        raise NameError('Mode not in mode_options')
        
    # (note that fast surface modes, i.e. 14 and 15, can only be 
    # found with SBS parameters in slab_functions...)
    mode = mode_options[mode_ind]    
    
    # Specify oscillation parameters
    if mode in slow_surf_mode_options + alfven_mode_options:
        K = 2.
    elif mode in slow_body_mode_options:
        K = 8.
    elif mode in fast_body_1_mode_options:
        K = 8.
    elif mode in fast_body_2_mode_options:
        K = 15.
    elif mode in fast_body_3_mode_options:
        K = 22.
    elif mode in fast_surf_mode_options:
        K = 8.
    else:
        raise NameError('Mode not found')

# Specify density ratio R1 := rho_1 / rho_0
#    R1 = 1.5 # Higher denisty on left than right
#    R1 = 1.8
#    R1 = 1.9 # Disp_diagram will only work for R1=1.5, 1.8, 2.0       
    R1 = 2. # Symmetric slab
    
    # Reduce number of variables in dispersion relation
    def disp_rel_asym_2var(W, K):
        return sf.disp_rel_asym(W, K, R1)
    
    # find eigenfrequencies W (= omega/k) within the range Wrange for the given parameters.
    Wrange1 = np.linspace(0., sf.cT, 11)
    Wrange2 = np.linspace(sf.cT, sf.c0, 401)
    Wrange3 = np.linspace(sf.c0, sf.c2, 11)
    
    Woptions_slow_surf = np.real(tool.point_find(disp_rel_asym_2var, np.array(K), Wrange1, args=None).transpose())
    Woptions_slow_body = np.real(tool.point_find(disp_rel_asym_2var, np.array(K), Wrange2, args=None).transpose())
    Woptions_fast = np.real(tool.point_find(disp_rel_asym_2var, np.array(K), Wrange3, args=None).transpose())
    
    # Remove W values that are very close to characteristic speeds - these are spurious solutions
    tol = 1e-2
    indices_to_rm = []
    for i in range(len(Woptions_slow_surf)):
        w = Woptions_slow_surf[i]
        if min(abs(np.array([w, w - sf.c0, w - sf.c1(R1), w - sf.c2, w - sf.vA]))) < tol or w < 0 or w > sf.cT:
            indices_to_rm.append(i)
    Woptions_slow_surf = np.delete(Woptions_slow_surf, indices_to_rm)
    Woptions_slow_surf.sort()
    
    indices_to_rm = []
    for i in range(len(Woptions_slow_body)):
        w = Woptions_slow_body[i]
        if min(abs(np.array([w, w - sf.c0, w - sf.c1(R1), w - sf.c2, w - sf.vA]))) < tol or w < sf.cT or w > sf.c0:
            indices_to_rm.append(i)
    Woptions_slow_body = np.delete(Woptions_slow_body, indices_to_rm)
    Woptions_slow_body.sort()
    
    indices_to_rm = []
    for i in range(len(Woptions_fast)):
        w = Woptions_fast[i]
        if min(abs(np.array([w, w - sf.c0, w - sf.c1(R1), w - sf.c2, w - sf.vA]))) < tol or w < sf.c0 or w > min(sf.c1, sf.c2):
            indices_to_rm.append(i)
    Woptions_fast = np.delete(Woptions_fast, indices_to_rm)
    Woptions_fast.sort()
    
    # remove any higher order slow body modes - we only want to do the first 3 saus/kink
    if len(Woptions_slow_body) > 6:
        Woptions_slow_body = np.delete(Woptions_slow_body, range(len(Woptions_slow_body) - 6))
    Woptions = np.concatenate((Woptions_slow_surf, Woptions_slow_body, Woptions_fast))
    
    # set W to be the eigenfrequency for the requested mode
    if mode in ['fast-saus-body-1', 'fast-saus-body-2', 'fast-saus-body-3', 'fast-kink-surf']:
        W = Woptions_fast[-2]
    elif mode in ['fast-kink-body-1', 'fast-kink-body-2', 'fast-kink-body-3', 'fast-saus-surf']:
        W = Woptions_fast[-1]
    elif mode in slow_surf_mode_options:
        W = Woptions_slow_surf[mode_ind]
    elif mode in slow_body_mode_options:
        W = Woptions_slow_body[mode_ind-2]
    
    if mode in alfven_mode_options:
        W = sf.vA
    else:
        W = np.real(W)
        
    # Quick plot to see if we are hitting correct mode    
    if show_quick_plot == True:
        plt.plot([K] * len(Woptions), Woptions, '.')
        plt.plot(K+0.5, W, 'go')
        plt.xlim([0,23])
        plt.show()

    #################################################################################
    
    # Plot dispersion diagram
    if show_dispersion == True:
        if mode in alfven_mode_options:
            raise NameError('Disperion plot requested for an alfven mode. Cant do it')
        
        dispersion_diagram.dispersion_diagram(mode_options, mode, 
                                              disp_rel_asym_2var, K, W, R1)
#        plt.tight_layout() # seems to make it chop the sides off with this
        plt.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_dispersion_diagrams\\'
                    + 'R1_' + str(R1) + '_' + mode + '.png')   
        plt.close()
    
    ##############################################################################
    
    if show_animation == True:
        print('Starting ' + mode)
        
        # set grid parameters
        xmin = -2.*K
        xmax = 2.*K
        ymin = 0.
        ymax = 4.
        zmin = 0.
        zmax = 2*np.pi
        
        # You can change ny but be careful changing nx, nz.
        nx = 300#100 #100 #300 gives us reduced bouncing of field lines for the same video size, but there is significant computational cost.
        ny = 300#100 #100 #100#20 #100
        nz = 300#100 #100
        nt = number_of_frames
        
        if nz % nt != 0:
            print("nt doesnt divide nz so there may be a problem with chopping in z direction for each time step")
        
        t_start = 0.
        t_end = zmax
        
        t = t_start
        
        xvals = np.linspace(xmin, xmax, nx)
        yvals = np.linspace(ymin, ymax, ny)
        zvals = np.linspace(zmin, zmax, nz, endpoint=False) # A fudge to give the height as exactly one wavelength

        x_spacing = max(nx, ny, nz) / nx
        y_spacing = max(nx, ny, nz) / ny
        z_spacing = max(nx, ny, nz) / nz        
        
        # For masking points for plotting vector fields- have to do it manually due to Mayavi bug
        mod = int(4 * nx / 100)
        mod_y = int(np.ceil(mod / y_spacing))
        
        # Get the data xi=displacement, v=velocity, b=mag field
        if show_disp_top == True or show_disp_front == True:
            xixvals = np.real(np.repeat(sf.xix(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            xizvals = np.real(np.repeat(sf.xiz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            xiyvals = np.real(np.repeat(sf.xiy(mode, xvals, zvals, t, W, K)[:, :, np.newaxis], ny, axis=2))
        
        if show_vel_front == True or show_vel_top == True:
            vxvals = np.real(np.repeat(sf.vx(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vzvals = np.real(np.repeat(sf.vz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vyvals = np.real(np.repeat(sf.vy(mode, xvals, zvals, t, K)[:, :, np.newaxis], ny, axis=2))
        
        if show_vel_front_pert == True or show_vel_top_pert == True:
            vxvals = np.real(np.repeat(sf.vx_pert(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vzvals = np.real(np.repeat(sf.vz_pert(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vyvals = np.zeros_like(vxvals)
        
        # Axis is defined on the mag field so we have to set up this data
        bxvals = np.real(np.repeat(sf.bx(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
        byvals = np.real(np.repeat(sf.by(mode, xvals, zvals, t, K)[:, :, np.newaxis], ny, axis=2))
        bz_eq3d = np.repeat(sf.bz_eq(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2)
        bzvals = np.real(np.repeat(-sf.bz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2) +
                         bz_eq3d)
        
        # displacement at the right and left boundaries
        if show_boundary == True:
            xix_boundary_r_vals = np.real(np.repeat(K + sf.xix_boundary(mode, zvals, t, W, K, R1, boundary='r')[:, np.newaxis], ny, axis=1))
            xix_boundary_l_vals = np.real(np.repeat(-K + sf.xix_boundary(mode, zvals, t, W, K, R1, boundary='l')[:, np.newaxis], ny, axis=1))
        
        if show_density == True:
            rho_vals = np.real(np.repeat(sf.rho(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
        
        if show_density_pert == True:
            rho_vals = np.real(np.repeat(sf.rho_pert(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
        
        bxvals_t = bxvals
        byvals_t = byvals
        bzvals_t = bzvals
        
        if show_disp_top == True or show_disp_front == True:
            xixvals_t = xixvals
            xiyvals_t = xiyvals
            xizvals_t = xizvals
            
        if np.array([show_vel_top, show_vel_top_pert, show_vel_front, show_vel_front_pert]).any() == True:
            vxvals_t = vxvals
            vyvals_t = vyvals
            vzvals_t = vzvals 
        
        if show_boundary == True:
            xix_boundary_r_vals_t = xix_boundary_r_vals
            xix_boundary_l_vals_t = xix_boundary_l_vals
            
        if show_density == True or show_density_pert == True:
            rho_vals_t = rho_vals
            
            
            
            
            
            
            
            
        zgrid_zy, ygrid_zy = np.mgrid[0:nz:(nz)*1j,
                          0:ny:(ny)*1j]
        
        #
        ##
        ###
        ####
        #####
        ####
        ###
        ##
        #
        
        fig = mlab.figure(size=res) # (1920, 1080) for 1080p , tuple(101 * np.array((16,9))) #16:9 aspect ratio for video upload

        # Spacing of grid so that we can display a visualisation cube without having the same number of grid points in each dimension
        spacing =  np.array([x_spacing, z_spacing, y_spacing])                     
                                             
        if show_density == True or show_density_pert == True:
            # Scalar field density   
            rho = mlab.pipeline.scalar_field(rho_vals_t, name="density", figure=fig)
            rho.spacing = spacing
            mpf.volume_red_blue(rho, rho_vals_t)
        
        #Masking points
        if show_mag_vec == True:
            bxvals_mask_front_t, byvals_mask_front_t, bzvals_mask_front_t = mpf.mask_points(bxvals_t, byvals_t, bzvals_t, 
                                                                                            'front', mod, mod_y)
        if show_disp_top == True:    
            xixvals_mask_top_t, xiyvals_mask_top_t, xizvals_mask_top_t = mpf.mask_points(xixvals_t, xiyvals_t, xizvals_t, 
                                                                                         'top', mod, mod_y)
        if show_disp_front == True:
            xixvals_mask_front_t, xiyvals_mask_front_t, xizvals_mask_front_t = mpf.mask_points(xixvals_t, xiyvals_t, xizvals_t, 
                                                                                               'front', mod, mod_y)
        if show_vel_top == True or show_vel_top_pert == True:  
            vxvals_mask_top_t, vyvals_mask_top_t, vzvals_mask_top_t = mpf.mask_points(vxvals_t, vyvals_t, vzvals_t, 
                                                                                      'top', mod, mod_y)                                                
        if show_vel_front == True or show_vel_front_pert == True:
            vxvals_mask_front_t, vyvals_mask_front_t, vzvals_mask_front_t = mpf.mask_points(vxvals_t, vyvals_t, vzvals_t, 
                                                                                            'front', mod, mod_y)        
        
        
        xgrid, zgrid, ygrid = np.mgrid[0:nx:(nx)*1j,
                                       0:nz:(nz)*1j,
                                       0:ny:(ny)*1j]
            
        field = mlab.pipeline.vector_field(bxvals_t, bzvals_t, byvals_t, name="B field", 
                                           figure=fig, scalars=zgrid)
        field.spacing = spacing        
        

        
        if show_axes == True:
            mpf.axes_no_label(field)
                
        if show_mini_axis == True:
            mpf.mini_axes()

        if uniform_light == True:
            #uniform lighting, but if we turn shading of volumes off, we are ok without
            mpf.uniform_lighting(fig)
        
        #Black background
        mpf.background_colour(fig, (0., 0., 0.))        
        
            #contours = mlab.pipeline.iso_surface(magnitude,
            #                                        contours=range(2, 14, 3),
            #                                        transparent=True,
            #                                        opacity=0.4,
            #                                        colormap='YlGnBu',
            #                                        vmin=0, vmax=14)
        
        
        scalefactor = 8. * nx / 100. # scale factor for direction field vectors
        
        # Set up visualisation modules
        if show_mag_vec == True:
            bdirfield_front = mlab.pipeline.vector_field(bxvals_mask_front_t, bzvals_mask_front_t,
                                                         byvals_mask_front_t, name="B field front",
                                                         figure=fig)
            bdirfield_front.spacing = spacing
            mpf.vector_cut_plane(bdirfield_front, 'front', nx, ny, nz, 
                                 y_spacing, scale_factor=scalefactor)
        
        if show_vel_top == True or show_vel_top_pert == True:
            vdirfield_top = mlab.pipeline.vector_field(vxvals_mask_top_t, np.zeros_like(vxvals_mask_top_t),
                                                       vyvals_mask_top_t, name="V field top",
                                                       figure=fig)
            vdirfield_top.spacing = spacing
            mpf.vector_cut_plane(vdirfield_top, 'top', nx, ny, nz, 
                                 y_spacing, scale_factor=scalefactor)
            
        if show_vel_front == True or show_vel_front_pert == True:
            vdirfield_front = mlab.pipeline.vector_field(vxvals_mask_front_t, vzvals_mask_front_t,
                                                         vyvals_mask_front_t, name="V field front",
                                                         figure=fig)
            vdirfield_front.spacing = spacing
            mpf.vector_cut_plane(vdirfield_front,'front', nx, ny, nz, 
                                 y_spacing, scale_factor=scalefactor)
        
        if show_disp_top == True:
            xidirfield_top = mlab.pipeline.vector_field(xixvals_mask_top_t, np.zeros_like(xixvals_mask_top_t),
                                                        xiyvals_mask_top_t, name="Xi field top",
                                                        figure=fig)
            xidirfield_top.spacing = spacing
            mpf.vector_cut_plane(xidirfield_top, 'top', nx, ny, nz, 
                                 y_spacing, scale_factor=scalefactor)
            
        if show_disp_front == True:
            xidirfield_front = mlab.pipeline.vector_field(xixvals_mask_front_t, xizvals_mask_front_t,
                                                         xiyvals_mask_front_t, name="Xi field front",
                                                         figure=fig)
            xidirfield_front.spacing = spacing
            mpf.vector_cut_plane(xidirfield_front, 'front', nx, ny, nz, 
                                 y_spacing, scale_factor=scalefactor)
        
        # Loop through time
        for t_ind in range(nt):
            if t_ind == 0:
                bxvals_t = bxvals
                byvals_t = byvals
                bzvals_t = bzvals
                
                if show_disp_top == True or show_disp_front == True:
                    xixvals_t = xixvals
                    xiyvals_t = xiyvals
                    xizvals_t = xizvals
                    
                if np.array([show_vel_top, show_vel_top_pert, show_vel_front, show_vel_front_pert]).any() == True:
                    vxvals_t = vxvals
                    vyvals_t = vyvals
                    vzvals_t = vzvals
                
                if show_boundary == True:
                    xix_boundary_r_vals_t = xix_boundary_r_vals
                    xix_boundary_l_vals_t = xix_boundary_l_vals
                    
                if show_density == True or show_density_pert == True:
                    rho_vals_t = rho_vals
                     
            else:
                
                bxvals = np.real(np.repeat(sf.bx(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
                byvals = np.real(np.repeat(sf.by(mode, xvals, zvals, t, K)[:, :, np.newaxis], ny, axis=2))
                bz_eq3d = np.repeat(sf.bz_eq(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2)
                bzvals = np.real(np.repeat(-sf.bz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2) +
                                 bz_eq3d)
                                 
                bxvals_t = bxvals
                byvals_t = byvals
                bzvals_t = bzvals
                
                # Update mag field data
                field.mlab_source.set(u=bxvals_t, v=bzvals_t, w=byvals_t)
                
                # Update mag field visualisation module
                if show_mag_vec == True:
                    bxvals_mask_front_t, byvals_mask_front_t, bzvals_mask_front_t = mpf.mask_points(bxvals_t, byvals_t, bzvals_t, 
                                                                                                    'front', mod, mod_y)
                    bdirfield_front.mlab_source.set(u=bxvals_mask_front_t, v=bzvals_mask_front_t, w=byvals_mask_front_t)                                                                                        
                
                # Update displacement field data
                if show_disp_top == True or show_disp_front == True:            
                    xixvals_split = np.split(xixvals, [nz - (nz / nt) * t_ind], axis=1)
                    xiyvals_split = np.split(xiyvals, [nz - (nz / nt) * t_ind], axis=1)
                    xizvals_split = np.split(xizvals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    xixvals_t = np.concatenate((xixvals_split[1], xixvals_split[0]), axis=1)
                    xiyvals_t = np.concatenate((xiyvals_split[1], xiyvals_split[0]), axis=1)
                    xizvals_t = np.concatenate((xizvals_split[1], xizvals_split[0]), axis=1)
                    
                    # Update displacement field visualisation module
                    if show_disp_top == True:
                        xixvals_mask_top_t, xiyvals_mask_top_t, xizvals_mask_top_t = mpf.mask_points(xixvals_t, xiyvals_t, xizvals_t, 
                                                                                                     'top', mod, mod_y)
                        xidirfield_top.mlab_source.set(u=xixvals_mask_top_t, v=np.zeros_like(xixvals_mask_top_t), w=xiyvals_mask_top_t)
                    if show_disp_front == True:
                        xixvals_mask_front_t, xiyvals_mask_front_t, xizvals_mask_front_t = mpf.mask_points(xixvals_t, xiyvals_t, xizvals_t, 
                                                                                                           'front', mod, mod_y)
                        xidirfield_front.mlab_source.set(u=xixvals_mask_front_t, v=xizvals_mask_front_t, w=xiyvals_mask_front_t)
                
                # Update velocity field data
                if np.array([show_vel_top, show_vel_top_pert, show_vel_front, show_vel_front_pert]).any() == True:            
                    vxvals_split = np.split(vxvals, [nz - (nz / nt) * t_ind], axis=1)
                    vyvals_split = np.split(vyvals, [nz - (nz / nt) * t_ind], axis=1)
                    vzvals_split = np.split(vzvals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    vxvals_t = np.concatenate((vxvals_split[1], vxvals_split[0]), axis=1)
                    vyvals_t = np.concatenate((vyvals_split[1], vyvals_split[0]), axis=1)
                    vzvals_t = np.concatenate((vzvals_split[1], vzvals_split[0]), axis=1)
                    
                    # Update velocity field visualisation module
                    if show_vel_top == True or show_vel_top_pert == True:
                        vxvals_mask_top_t, vyvals_mask_top_t, vzvals_mask_top_t = mpf.mask_points(vxvals_t, vyvals_t, vzvals_t, 
                                                                                                  'top', mod, mod_y)                                                
                        vdirfield_top.mlab_source.set(u=vxvals_mask_top_t, v=np.zeros_like(vxvals_mask_top_t), w=vyvals_mask_top_t)
                    if show_vel_front == True or show_vel_front_pert == True:
                        vxvals_mask_front_t, vyvals_mask_front_t, vzvals_mask_front_t = mpf.mask_points(vxvals_t, vyvals_t, vzvals_t, 
                                                                                                        'front', mod, mod_y) 
                        vdirfield_front.mlab_source.set(u=vxvals_mask_front_t, v=vzvals_mask_front_t, w=vyvals_mask_front_t)
                
                # Update boundary displacement data
                if show_boundary == True:
                    xix_boundary_r_vals_split = np.split(xix_boundary_r_vals, [nz - (nz / nt) * t_ind], axis=0)
                    xix_boundary_l_vals_split = np.split(xix_boundary_l_vals, [nz - (nz / nt) * t_ind], axis=0)
        
                    xix_boundary_r_vals_t = np.concatenate((xix_boundary_r_vals_split[1], xix_boundary_r_vals_split[0]), axis=0)
                    xix_boundary_l_vals_t = np.concatenate((xix_boundary_l_vals_split[1], xix_boundary_l_vals_split[0]), axis=0)
                
                # Update density data
                if show_density == True or show_density_pert == True:            
                    rho_vals_split = np.split(rho_vals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    rho_vals_t = np.concatenate((rho_vals_split[1], rho_vals_split[0]), axis=1)    

                    rho.mlab_source.set(scalars=rho_vals_t)                     
            
            # Boundary data - Letting mayavi know where to plot the boundary
            if show_boundary == True:
                ext_min_r = ((nx) * (xix_boundary_r_vals_t.min() - xmin) / (xmax - xmin)) * x_spacing
                ext_max_r = ((nx) * (xix_boundary_r_vals_t.max() - xmin) / (xmax - xmin)) * x_spacing
                
                ext_min_l = ((nx) * (xix_boundary_l_vals_t.min() - xmin) / (xmax - xmin)) * x_spacing
                ext_max_l = ((nx) * (xix_boundary_l_vals_t.max() - xmin) / (xmax - xmin)) * x_spacing
                                                   
 
            #Make field lines
            if show_mag == True:
                # move seed points up with phase speed. - Bit of a fudge.
                # Create an array of points for which we want mag field seeds
                nx_seed = 9
                ny_seed = 13
                start_x = 30. * nx / 100.
                end_x = nx+1 - start_x
                start_y = 1.
                if ny == 20: # so that the lines dont go right up to the edge of the box
                    end_y = ny - 1.
                elif ny == 100:
                    end_y = ny - 2.
                elif ny == 300:
                    end_y = ny - 6.
                else:
                    end_y = ny - 1
                seeds=[]
                dx_res = (end_x - start_x) / (nx_seed-1)
                dy_res = (end_y - start_y) / (ny_seed-1)
                for j in range(0,ny_seed):
                    for i in range(0,nx_seed):
                        x = start_x + (i * dx_res) * x_spacing
                        y = start_y + (j * dy_res) * y_spacing
                        z = 1. + (t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
                        seeds.append((x,z,y))
                
                if mode in alfven_mode_options:
                    for i in range(nx_seed):
                        del seeds[0]
                        del seeds[-1]
                    
                # Remove previous field lines - field lines cannot be updated, just the data that they are built from
                if t_ind != 0:
                    field_lines.remove() # field_lines is defined in first go through loop
                field_lines = SeedStreamline(seed_points=seeds)
                
                # Field line visualisation tinkering
                field_lines.stream_tracer.integration_direction='both'
                field_lines.streamline_type = 'tube'
                field_lines.stream_tracer.maximum_propagation = nz * 2
                field_lines.tube_filter.number_of_sides = 20
                field_lines.tube_filter.radius = 0.7 * max(nx, ny, nz) / 100.
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

            
            # Which views do you want to show? Options are defined at the start
            views_selected = [2,3]#[0,1,4,5,6] #range(7) #[2,3]
            for view_ind in range(len(views_selected)):
                view = view_options[views_selected[view_ind]]
                
                # Display boundary - cannot be updated each time
                if show_boundary == True:
                    # Boundaries should look different depending on view
                    if view == 'front-parallel':
                        #remove previous boundaries
                        if t != 0 or view_ind != 0:
                            boundary_r.remove()
                            boundary_l.remove()
                        
                        # Make a fading colormap by changing opacity at ends
                        lut = np.reshape(np.array([150, 150, 150, 255]*256), (256,4))
                        fade_value = 125
                        lut[:fade_value,-1] = np.linspace(0, 255, fade_value)
                        lut[-fade_value:,-1] = np.linspace(255, 0, fade_value)
                        
                        # Set up boundary visualisation
                        boundary_r = mlab.mesh(xix_boundary_r_vals_t, zgrid_zy, ygrid_zy,
                                                     extent=[ext_min_r, ext_max_r, 1, nz, 0, (ny-1) * y_spacing],
                                                     opacity=1., representation='wireframe',
                                                     line_width=12., scalars=zgrid_zy)
                        boundary_l = mlab.mesh(xix_boundary_l_vals_t, zgrid_zy, ygrid_zy,
                                                     extent=[ext_min_l, ext_max_l, 1, nz, 0, (ny-1) * y_spacing],
                                                     opacity=1., representation='wireframe',
                                                     line_width=12., scalars=zgrid_zy)
                        
                        # Boundary color and other options
                        boundary_r.module_manager.scalar_lut_manager.lut.table = lut
                        boundary_l.module_manager.scalar_lut_manager.lut.table = lut
                        boundary_r.actor.property.lighting = False
                        boundary_r.actor.property.shading = False
                        boundary_l.actor.property.lighting = False
                        boundary_l.actor.property.shading = False
                                                     
                    else:
                        #remove previous boundaries
                        if t != 0 or view_ind != 0:
                            boundary_r.remove()
                            boundary_l.remove()
                            
                        # Make a fading colormap by changing opacity at ends
                        lut = np.reshape(np.array([150, 150, 150, 255]*256), (256,4))
                        fade_value = 20
                        lut[:fade_value,-1] = np.linspace(0, 255, fade_value)
                        lut[-fade_value:,-1] = np.linspace(255, 0, fade_value)
                        
                        # Set up boundary visualisation
                        boundary_r = mlab.mesh(xix_boundary_r_vals_t, zgrid_zy, ygrid_zy,
                                               extent=[ext_min_r, ext_max_r, 1, nz, 0, (ny-1) * y_spacing],
                                               opacity=0.7, scalars=zgrid_zy)
        
                        boundary_l = mlab.mesh(xix_boundary_l_vals_t, zgrid_zy, ygrid_zy,
                                               extent=[ext_min_l, ext_max_l, 1, nz, 0, (ny-1) * y_spacing],
                                               opacity=0.7, scalars=zgrid_zy)
                        
                        # Boundary color and other options
                        boundary_r.module_manager.scalar_lut_manager.lut.table = lut
                        boundary_l.module_manager.scalar_lut_manager.lut.table = lut
                        boundary_r.actor.property.lighting = False
                        boundary_r.actor.property.shading = False
                        boundary_l.actor.property.lighting = False
                        boundary_l.actor.property.shading = False       

                # Set viewing angle - For some unknown reason we must redefine the camera position each time.
                # This is something to do with the boundaries being replaced each time.
                mpf.view_position(fig, view, nx, ny, nz)
                
                if save_images == True:
                    prefix = 'R1_'+str(R1) + '_' + mode + '_' + vis_mod_string + view + '_'# + '_norho_'
                    mlab.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_animations\\'
                                 + prefix + str(t_ind+1) + '.png')
                if t_ind == nt - 1:
                    if make_video == True:
                        i2v.image2video(prefix=prefix, output_name=prefix+'video', 
                                        out_extension='mp4', fps=fps, n_loops=4, 
                                        delete_images=True, delete_old_videos=True, res=res[1])
            
            # Keep us updated with progress
            if t_ind % 5 == 4:
                print('Finished frame number ' + str(t_ind + 1) + ' out of ' + str(number_of_frames))

            #Release some memory after each time step
            gc.collect()
            
            #step t forward
            t = t + (t_end - t_start) / nt
            
        # Close Mayavi window each time if we cant to make a video
        if make_video == True:
            mlab.close(fig)
        print('Finished ' + mode)
                            