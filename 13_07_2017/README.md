# 3D visualisations of MHD waves in asymmetric an asymmetric slab

## Prerequisites:
[Python 2](https://www.python.org/download/releases/2.7.2/),
[Mayavi](http://mayavi.sourceforge.net/),
[ffmpeg](https://ffmpeg.org/).

## 3D_slab_modes.py:
Builds an animation of the magnetoacoustic eigenmodes of an isolated asymmetric magnetic slab with possible visualisations of:
* magnetic field lines,
* magentic field strength,
* magnetic field vectors,
* velocity (eulerian and lagrangian),
* displacement,
* density (eularian and lagrangian),
* perturbed slab boundary.

## dispersion_diagram.py:
Functions that plot dispersion diagrams (omega/k against kx_0).

## img2vid.py:
Uses *ffmpeg* command line operations to make a set of images into a video.

## alfven_mixed_driver.py:
Builds an animation of a shear alfven wave in a plasma with spacially variable alfven speed (and otherwise uniform plasma).

We use the analytical solutions given in our paper which can be found [here](http://link.springer.com/article/10.1007/s11207-017-1054-y).