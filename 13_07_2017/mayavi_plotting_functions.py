# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 16:16:47 2017

@author: Matt
"""

# Functions that will initiate mayavi plotting modules

import numpy as np
from mayavi import mlab

def uniform_lighting(scene):
    scene.scene.light_manager.number_of_lights = 6
    
    camera_light1 = scene.scene.light_manager.lights[0]
    camera_light1.activate = True
    camera_light1.intensity = 0.7
    camera_light1.elevation = 90.
    camera_light1.azimuth = 0.

    camera_light2 = scene.scene.light_manager.lights[1]
    camera_light2.activate = True
    camera_light2.intensity = 0.7
    camera_light2.elevation = -90.
    camera_light2.azimuth = 0.

    camera_light3 = scene.scene.light_manager.lights[2]
    camera_light3.activate = True
    camera_light3.intensity = 0.7
    camera_light3.elevation = 0.
    camera_light3.azimuth = -90

    camera_light4 = scene.scene.light_manager.lights[3]
    camera_light4.activate = True
    camera_light4.intensity = 0.7
    camera_light4.elevation = 0.
    camera_light4.azimuth = 0.

    camera_light5 = scene.scene.light_manager.lights[4]
    camera_light5.activate = True
    camera_light5.intensity = 0.7
    camera_light5.elevation = 0.
    camera_light5.azimuth = 90.

    camera_light6 = scene.scene.light_manager.lights[5]
    camera_light6.activate = True
    camera_light6.intensity = 0.7
    camera_light6.elevation = 0.
    camera_light6.azimuth = 180.

def background_colour(scene, colour):
    """
    colour = tuple colour code
    e.g (0.,0.,0.) for black
    """
    scene.scene.background = colour

def mini_axes():
    oa = mlab.orientation_axes(xlabel='x', ylabel='z', zlabel='y')
    oa.marker.set_viewport(0,0,0.25,0.25) # minx, miny, maxx, maxy

def axes(scene, show_axis_labels, view):
    axes = mlab.axes(scene, nb_labels=1, line_width=3)
    axes.axes.label_format = ''
    if show_axis_labels == True:
        axes.axes.x_label = 'x'
        if view == 'front-parallel':
            axes.axes.y_label = 'z'
            axes.axes.z_label = ''
        elif view == 'top-parallel':
            axes.axes.y_label = ''
            axes.axes.z_label = 'y'
        else:
            axes.axes.y_label = 'z'
            axes.axes.z_label = 'y'
    else:
        axes.axes.x_label = ''
        axes.axes.y_label = ''
        axes.axes.z_label = ''
        
def axes_no_label(scene):
    axes = mlab.axes(scene, nb_labels=1, line_width=3)
    axes.axes.label_format = ''
    axes.axes.x_label = ''
    axes.axes.y_label = ''
    axes.axes.z_label = ''
        
def view_position(scene, view, nx, ny, nz):
    #Set viewing angle
    if view == 'front-parallel':
        scene.scene.z_plus_view()
        scene.scene.parallel_projection = True
        scene.scene.camera.zoom(1.65) # Parallel projection zoom is done in this way, different to perspective projection
    if view == 'front':
#        scene.scene.z_plus_view()
#        scene.scene.parallel_projection = True
#        scene.scene.camera.zoom(1.65)
##        scene.scene.camera.view_angle = 21.
        scene.scene.parallel_projection = False
        scene.scene.camera.position = [50.5 * nx / 100., 50.5 * nz / 100., 382.37955413300307  * ny / 100.]
        scene.scene.camera.focal_point = [50.5 * nx / 100., 50.5 * nz / 100., 50.0 * ny / 100.]
        scene.scene.camera.view_angle = 21.0
        scene.scene.camera.view_up = [0.0, 1.0, 0.0]
        scene.scene.camera.clipping_range = [200, 1000]#229.55575859167305 * 2, 462.61524744499809 * 2]

    if view == 'top':
        scene.scene.parallel_projection = False
        scene.scene.camera.position = [53.10 * nx / 100., 523.36 * nz / 100., 50.95 * ny / 100.]
        scene.scene.camera.focal_point = [50.82 * nx / 100., 50.41 * nz / 100., 50.16 * ny / 100.]
        scene.scene.camera.view_angle = 14.
        scene.scene.camera.view_up = [-0, 0, -1]
        scene.scene.camera.clipping_range = [368.83 * nx / 100., 605.15 * nx / 100.]
    if view == 'top-parallel':
        scene.scene.parallel_projection = True
        scene.scene.camera.zoom(1.65)
        scene.scene.camera.position = [53.11 * nx / 100., 523.36 * nz / 100., 50.95 * ny / 100.]
        scene.scene.camera.focal_point = [50.82 * nx / 100., 50.41 * nz / 100., 50.16 * ny / 100.]
#            scene.scene.camera.view_angle = 14.
        scene.scene.camera.view_up = [0, 0, -1]
        scene.scene.camera.clipping_range = [368.83 * nx / 100., 605.15 * nx / 100.]
    if view == 'front-top':
        scene.scene.parallel_projection = False
        scene.scene.camera.position = [48.76 * nx / 100., 223.65 * nz / 100., 498.62 * ny / 100.]
        scene.scene.camera.focal_point = [50.82 * nx / 100., 46. * nz / 100., 50.16 * ny / 100.]
        scene.scene.camera.view_angle = 16.0
        scene.scene.camera.view_up = [-0.002418791139063777, 0.93281530024654913, -0.36034672896443193]
        scene.scene.camera.clipping_range = [345.98 * nx / 100., 650.72 * nx / 100.]
     
    if view == 'front-side':
        scene.scene.parallel_projection = False
        scene.scene.camera.position = [126.6 * nx / 100., 60.5 * nz / 100., 524.8 * ny / 100.]
        scene.scene.camera.focal_point = [50.8 * nx / 100., 50.4 * nz / 100., 50.2 * ny / 100.]
        scene.scene.camera.view_angle = 14.
        scene.scene.camera.view_up = [-0.01695, 0.999686, -0.0184509]
        scene.scene.camera.clipping_range = [366.21 * nx / 100., 631.08 * nx / 100.]
        
    if view == 'front-top-side':
        scene.scene.parallel_projection = False
        scene.scene.camera.position = [400.87 * nx / 100., 181.10 * nz / 100., 495.97 * ny / 100.]
        scene.scene.camera.focal_point = [50.80 * nx / 100., 43. * nz / 100., 50.20 * ny / 100.]
        scene.scene.camera.view_angle = 13.0
        scene.scene.camera.view_up = [-0.14482357103326407, 0.9744012643898321, -0.17195438124301951]
        scene.scene.camera.clipping_range = [418.37 * nx / 100., 789.31 * nx / 100.]
        

def mask_points(var_x, var_y, var_z, front_or_top, mod, mod_y):
    if front_or_top == 'front':
        var_x_mask = np.copy(var_x)
        var_y_mask = np.copy(var_y)
        var_z_mask = np.copy(var_z)
        
        for i in range(var_x_mask.shape[0]):
            for j in range(var_x_mask.shape[1]):
                for k in range(var_x_mask.shape[2]):
                    if (i%mod) != 1 or (j%mod) != 1:
                        var_x_mask[i,j,k] = 0.
                        var_y_mask[i,j,k] = 0.
                        var_z_mask[i,j,k] = 0.
    elif front_or_top == 'top':
        var_x_mask = np.copy(var_x)
        var_y_mask = np.copy(var_y)
        var_z_mask = np.copy(var_z)
        
        for i in range(var_x_mask.shape[0]):
            for j in range(var_x_mask.shape[1]):
                for k in range(var_x_mask.shape[2]):
                    if (i%mod) != 1 or (k%mod_y) != 1:
                        var_x_mask[i,j,k] = 0.
                        var_y_mask[i,j,k] = 0.
                        var_z_mask[i,j,k] = 0.
    else:
        raise ValueError("front_or_top can be only 'front' or 'top'")
    return [var_x_mask, var_y_mask, var_z_mask]

def vector_cut_plane(vec_field, front_or_top, nx, ny, nz, y_spacing, scale_factor):
#    scale_factor = 4. # scale factor for direction field vectors
    if front_or_top == 'front':
        vector_cut_plane_front = mlab.pipeline.vector_cut_plane(vec_field, 
                                                                scale_factor=scale_factor)
        vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
        vector_cut_plane_front.implicit_plane.widget.origin = np.array([ nx / 2., nz / 4., (ny-1)*y_spacing])
        vector_cut_plane_front.glyph.color_mode = 'no_coloring'
        vector_cut_plane_front.implicit_plane.widget.enabled = False
        vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
        vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
    
    if front_or_top == 'top':
        vector_cut_plane_top = mlab.pipeline.vector_cut_plane(vec_field, 
                                                              scale_factor=scale_factor)
        vector_cut_plane_top.implicit_plane.widget.normal_to_y_axis = True
        vector_cut_plane_top.glyph.color_mode = 'no_coloring'
        vector_cut_plane_top.implicit_plane.widget.origin = np.array([ nx / 2.,nz-0.1, nx / 2.])
        vector_cut_plane_top.implicit_plane.widget.enabled = False
        vector_cut_plane_top.glyph.glyph_source.glyph_source = vector_cut_plane_top.glyph.glyph_source.glyph_dict['arrow_source']
        vector_cut_plane_top.glyph.glyph_source.glyph_position = 'center'
    
    
def volume_red_blue(scalar_field, scalar_data):
    minr = scalar_data.min()
    maxr = scalar_data.max()
    
    #Volume for high pressure
    rvmin1 = minr + 0.5 * (maxr - minr)
    rvmax1 = minr + 1. * (maxr - minr)
    rvol1 = mlab.pipeline.volume(scalar_field, vmin=rvmin1, vmax=rvmax1)
    
    # Changing the ctf:
    from tvtk.util.ctf import ColorTransferFunction
    ctf1 = ColorTransferFunction()
    ctf1.add_rgb_point(rvmin1, 1., 1., 0.5)
    ctf1.add_rgb_point(rvmin1 + 0.4 * (rvmax1 - rvmin1), 1, 0.3, 0.1)
    ctf1.add_rgb_point(rvmax1, 1., 0., 0.)
    # ...
    rvol1._volume_property.set_color(ctf1)
    rvol1._ctf = ctf1
    rvol1.update_ctf = True
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    from tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(rvmin1, 0)
    otf.add_point(rvmin1 + (rvmax1-rvmin1)*0.2, 0.012)
    otf.add_point(rvmin1 + (rvmax1-rvmin1)*0.5, 0.05)
    otf.add_point(rvmax1, 0.15)
    ##vol1._otf = otf
    rvol1._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol1.volume_property.shade = False
    rvol1.volume_property.scalar_opacity_unit_distance = 2.0
    
    
    #Volume for low pressure
    rvmin2 = minr + 0. * (maxr - minr)
    rvmax2 = minr + 0.5 * (maxr - minr)
    rvol2 = mlab.pipeline.volume(scalar_field, vmin=rvmin2, vmax=rvmax2)
    
    # Changing the ctf:
    ctf2 = ColorTransferFunction()
    ctf2.add_rgb_point(rvmin2, 0., 0.5, 1.)
    ctf2.add_rgb_point(rvmin2 + 0.6 * (rvmax2 - rvmin2), 0.1, 0.7, 1.)
    ctf2.add_rgb_point(rvmax2, 0.5, 1., 1.)
    # ...
    rvol2._volume_property.set_color(ctf2)
    rvol2._ctf = ctf2
    rvol2.update_ctf = True

    #Changing the opacity of the volume vol2
    ## Changing the otf:
    otf = PiecewiseFunction()
    otf.add_point(rvmax2, 0)
    otf.add_point(rvmax2 - (rvmax2-rvmin2)*0.2, 0.012)
    otf.add_point(rvmax2 - (rvmax2-rvmin2)*0.5, 0.05)
    otf.add_point(rvmin2, 0.15)
    ##vol1._otf = otf
    rvol2._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol2.volume_property.shade = False
    rvol2.volume_property.scalar_opacity_unit_distance = 2.0
    
def colormap_fade(module_manager, fade_value=20):
    lut = module_manager.scalar_lut_manager.lut.table.to_array()
    lut[:fade_value,-1] = np.linspace(0, 255, fade_value)
    lut[-fade_value:,-1] = np.linspace(255, 0, fade_value)
    module_manager.scalar_lut_manager.lut.table = lut
