import toolbox_extended as t
import numpy as np
import scipy.optimize as sp
from functools import partial

# Attempts to find the roots of a defined function in a box defined by two
# lists defined by numpy.linspace, using scipy.optimize.newton.
# Arguments:
# func - The function to be solved. This is must be a function of AT LEAST
# two variables, but it may have any number.
# x_range, y_range - (numpy.linspace) Lists defining the box.
# args - The values of the arguments of func. If this is omitted, the point
# finder assumes that func is only dependent on two variables, i.e. that 
# it is of the form f(a,b). The point finder will iterate over 'b' and find 
# the values of 'a' for which f(a,b) = 0. The point finder will ALWAYS find the
# roots of the FIRST variable in func. If the function depends on more than two
# variables, the variables must be imputed IN ORDER. This means that if you
# have a function f(a, b, c), the point finder will always use its root solver
# on 'a', and depending on how you input 'b' and 'c', it will keep one constant
# and iterate over the other. If, for example, your input is args=(None, 1), 
# the point finder will find the values of 'a' for each value of 'b' in x_range
# around the points in y_range, while keeping 'c' constant, with c=1. 
# If, on the other hand, your input is args=(1, None), it will find 'a', while
# keeping 'b' constant, with b=1, and iterate over 'c', where 'c' will take the
# values in y_range.

# Example:
# Suppose you have a function f(a, b, c, d) = 0.
# The point finder will attempt to use Newton's method to find values of the
# variable 'a' such that f(a) = 0 for some set values of 'b', 'c', and 'd',
# where one of these is iterated over. 
# Say you want b=1, c to be iterated over, and d=3. You would input this as
# y_pts = point_finder_scipy(func, x_range, y_range, args=(1, None, 3))
# Your output will be a numpy array called y_pts.
# You may then use plot(x_range, y_pts).

def point_find(func, x_range, y_range, args=None):
    y_pts = []
    if args is None:
        args = [None]
    if x_range.shape == ():
        out = []
        for a in args:
            out.append(float(x_range))
        root_temp = []
        for y in y_range:
            try:
                root = sp.newton(func, y, args = tuple(out), tol = 1e-20)
                root_temp = t.check(root, root_temp, 10**(-6))
            except RuntimeError:
                pass
        y_pts.append(root_temp)
    else:
        for x in x_range:
            out = []
            for a in args:
                if a is None:
                    a = x
                out.append(a)
            root_temp = []
            for y in y_range:
                try:
                    root = sp.newton(func, y, args = tuple(out), tol = 1e-20)
                    root_temp = t.check(root, root_temp, 10**(-6))
                except RuntimeError:
                    pass
            y_pts.append(root_temp)
    y_array = t.list_to_array(y_pts)
    return y_array

# Attempts to find a root that is close to a root inputted into the finder.
# Arguments:
# func - The function. Must depend on a single variable.
# root_list - List of roots of the function found so far. Can be empty.
# x - The x coordinate of the starting point.
# y - The y coordinate of the starting point.
# step_size - The step size to be used in calculating the next point.
# x_values -

def root_find(func, root_list, x, y, step_size, x_values):
    if len(root_list) == 0:
        root = sp.newton(func, y, maxiter = 2000)
        root_list.append(root)
        x_values.append(x)
    elif len(root_list) == 1:
        root = sp.newton(func, root_list[-1], maxiter = 2000)
        root_list.append(root)
        x_values.append(x)
    elif len(root_list) == 2:
        grad = ((root_list[-1] - root_list[-2])*
        np.abs(step_size/(x_values[-1] - x_values[-2])))
        root = sp.newton(func, root_list[-1] + grad, maxiter = 2000)
        if np.abs(root-root_list[-1]) > 0.1:
            raise Exception('largegrad at {}'.format(x_values[-1]))
        else:
            root_list.append(root)
            x_values.append(x)
    else:
        grad = ((root_list[-1] - root_list[-2]) + \
        (root_list[-1] - 2*root_list[-2] + root_list[-3]))* \
        np.abs(step_size/(x_values[-1] - x_values[-2]))
        root = sp.newton(func, root_list[-1] + grad, maxiter = 2000)
        if np.abs(root-root_list[-1]) > 0.1:
            raise Exception('largegrad at {}'.format(x_values[-1]))
        else:
            root_list.append(root)
            x_values.append(x)
    return x_values, root_list

# Attempts to trace a line of roots of a function.
# Arguments:
# func - The function. This can have any number of variables, 
# but the first one is always the one that the root finder will be used on.
# x - The x coordinate of the starting point.
# y - The y coordinate of the starting point.
# step_size - The step size to be used in calculating the next point.
# x_end_left, x_end_right - These define the limits of the interval in which
# line is to be traced.

def line_trace(func, x, y, step_size, x_end_left, x_end_right, args=None):
    if args is None:
        args = [None]
    x_values = []
    root_list = []
    while x > x_end_left:
        out = []
        for a in args:
            if a is None:
                a = x
            out.append(a)
        flip = lambda *args: func(*args[::-1])
        flip_f = partial(flip, *out[::-1])
        x_values, root_list = root_find(flip_f, root_list, x, y, -step_size, x_values)
        x -= step_size
    x = x_values[-1] + step_size
    x_values_new = []
    root_list_new = []
    while x <= x_end_right:
        out = []
        for a in args:
            if a is None:
                a = x
            out.append(a)
        flip = lambda *args: func(*args[::-1])
        flip_f = partial(flip, *out[::-1])
        x_values, root_list = root_find(flip_f, root_list, x, y, step_size, x_values)
        x_values_new.append(x_values[-1])
        root_list_new.append(root_list[-1])
        x += step_size
    root_array = np.array(root_list_new)
    return x_values_new, root_array
