
import numpy as np
import scipy as sc
from scipy.optimize import fsolve

#Different modes can be found with different characteristic speed orderings
# Choose wisely!
#
## SBB
## Define the sound speeds and alfven speeds.
#c2 = 1.2
#c0 = 1.
#vA = 0.9

#
## for xi of x slow surface and maybe others
#c2 = 0.7
#c0 = 1.
#vA = 0.4


# SBS
# Define the sound speeds and alfven speed and tube speed.
c2 = 1.2
c0 = 1.
vA = 1.3
cT = sc.sqrt(c0**2 * vA**2*(c0**2 + vA**2)**(-1))

mode_options = ['slow-kink-surf', 'slow-saus-surf', 'slow-saus-body-3',
                'slow-kink-body-3', 'slow-saus-body-2', 'slow-kink-body-2', 
                'slow-saus-body-1', 'slow-kink-body-1', 'fast-saus-body-1',
                'fast-kink-body-1', 'fast-saus-body-2', 'fast-kink-body-2',
                'fast-saus-body-3', 'fast-kink-body-3', 'fast-kink-surf',
                'fast-saus-surf', 'shear-alfven', 'shear-alfven-broadband']

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

def alfven_shear_width(mode, K):
    if mode == 'shear-alfven':
        return [-8*K/16, -5*K/16]#[-14*K/16, -13*K/16]
    elif mode == 'shear-alfven-broadband':
        return [-K, K]
    
alfven_amplitude = 0.3 # Amplitude of alfven modes

R2 = 2. # rho2/rho0

# To maintain equilibrium pressure balance the sound speed on one side is:
def c1(R1):
    return c2 * sc.sqrt(R2 / R1)

# Equilibrium magnetic field strength within the slab:
B0 = 1.

# variables:
# x = k*x
# y = k*y
# z = k*z
# W = omega/k
# K = k*x_0
# t = omega*t

# Helpful functions fr defining the eigenfunctions
def m0(W):
    return sc.sqrt((c0**2 - W**2)*(vA**2 - W**2) / 
                   ((c0**2 + vA**2)*(cT**2 - W**2)))
                   
def m00(W):
    return sc.sqrt(1 - W**2/c0**2)
    
def m1(W, R1):
    return sc.sqrt(1 - W**2*(c2*sc.sqrt(R2/R1))**(-2))
    
def m2(W):
    return sc.sqrt(1 - W**2/c2**2)

def lamb0(W):
    return -(vA**2-W**2)*1.j/(m0(W)*W)

def lamb00(W):
    return W*1.j/m00(W)

def lamb1(W, R1):
    return R1*W*1.j/m1(W, R1)
    
def lamb2(W):
    return R2*W*1.j/m2(W)

# The displacement required to make the mode look good. 
# Some modes have much smaller amplitude than others for the same parameters 
# so this accounts for that.
def required_xi(mode, K):
    if mode in slow_surf_mode_options:
        return K / 3.
    elif mode in slow_body_1_mode_options:
        return K / 130. #30.
    elif mode in slow_body_2_mode_options:
        return K / 140. #90.
    elif mode in slow_body_3_mode_options:
        return K / 250.
    elif mode in fast_body_1_mode_options:
        return K / 80.
    elif mode in fast_body_2_mode_options:
        return K / 250. #180.
    elif mode in fast_body_3_mode_options:
        return K / 1000. #400. #250.
    elif mode in fast_surf_mode_options:
        return K / 30. #40.
    else:
        print('Not a recognised mode')

def const(mode, W, K, R1):
    return 0.6
#    const_val_r = (W * required_xi(mode, K) / (constB_dash(mode, W, K, R1)*sc.cosh(m0(W)*K) +
#                                     constC_dash(mode, W, K, R1)*sc.sinh(m0(W)*K)))
#    const_val_l = (W * required_xi(mode, K) / (constB_dash(mode, W, K, R1)*sc.cosh(m0(W)*-K) +
#                                     constC_dash(mode, W, K, R1)*sc.sinh(m0(W)*-K)))
#    return min(const_val_r, const_val_l)
#print('WARNING: Limited const in slab_functions')
#    if mode in slow_surf_mode_options:
#        return 0.05
#    elif mode in slow_body_1_mode_options:
#        return 0.23
#    elif mode in slow_body_2_mode_options:
#        return 0.1
#    elif mode in slow_body_3_mode_options:
#        return 0.1
#    elif mode in fast_body_1_mode_options:
#        return 0.9

# Dispersion relation for an asymmetric magnetic slab
def disp_rel_asym(W, K, R1):
    return ((W**4*m0(W)**2*R1*R2 + (vA**2 - W**2)**2*m1(W, R1)*m2(W) -
            0.5*m0(W)*W**2*(vA**2 - W**2)*(R2*m1(W, R1) + R1*m2(W))*
            (sc.tanh(m0(W)*K) + (sc.tanh(m0(W)*K))**(-1))) /
            (vA**2 - W**2)*(c0**2 - W**2)*(cT**2 - W**2))

# Equilibrium magnetic field - 1 dimensional
def bz_eq_1d(mode, x, z, t, W, K, R1):
    truth = np.array((x <= (K + xix_boundary(mode, z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                 (x >= (-K + xix_boundary(mode, z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
    indices = np.where(truth == True)
    bz_eqfunction = np.zeros(len(x))
    for i in indices:
        bz_eqfunction[i] = B0
    return bz_eqfunction
    
# Equilibrium magnetic field - 2 dimensional
def bz_eq(mode, x, z, t, W, K, R1):
    bz_eq_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bz_eq_array[:,i] = bz_eq_1d(mode, x, z[i], t, W, K, R1)
    return bz_eq_array


##############################################

# constants in the eigenfunctions
def constB(mode, W, K, R1):
    if mode in kink_mode_options:
        return const(mode, W, K, R1)
    elif mode in saus_mode_options:
        return const(mode, W, K, R1)*((lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)) /
                               (lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)))

def constC(mode, W, K, R1):
    if mode in kink_mode_options:
        return const(mode, W, K, R1)*((lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)) /
                            (lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)))
    elif mode in saus_mode_options:
        return const(mode, W, K, R1)
    
def constA(mode, W, K, R1):
    return ((constB(mode, W, K, R1)*sc.cosh(m0(W)*K) - constC(mode, W, K, R1)*sc.sinh(m0(W)*K)) /
            (sc.cosh(m1(W, R1)*K) - sc.sinh(m1(W, R1)*K)))

def constD(mode, W, K, R1):
    return ((constB(mode, W, K, R1)*sc.cosh(m0(W)*K) + constC(mode, W, K, R1)*sc.sinh(m0(W)*K)) /
            (sc.cosh(m2(W)*K) - sc.sinh(m2(W)*K)))

# constants without const(mode). Used for finding appropriate value for const(mode)
def constB_dash(mode, W, K, R1):
    if mode in kink_mode_options:
        return 1.
    elif mode in saus_mode_options:
        return ((lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)) /
                               (lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)))

def constC_dash(mode, W, K, R1):
    if mode in kink_mode_options:
        return ((lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)) /
                            (lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)))
    elif mode in saus_mode_options:
        return 1.

# velocity amplitude
def vxhat(mode, x, W, K, R1):
    if mode in alfven_mode_options:
        vxhatfunction = np.zeros_like(x)
    else:
        if type(x) == float or type(x) == np.float64:
            if np.abs(x) <= K:
                vxhatfunction = (constB(mode, W, K, R1)*sc.cosh(m0(W)*x) + 
                                constC(mode, W, K, R1)*sc.sinh(m0(W)*x))
            elif x < -K:
                vxhatfunction = constA(mode, W, K, R1)*(sc.cosh(m1(W, R1)*x) + 
                                                      sc.sinh(m1(W, R1)*x))
            elif x > K:
                vxhatfunction = constD(mode, W, K, R1)*(sc.cosh(m2(W)*x) - 
                                                      sc.sinh(m2(W)*x))
        else:
            truth = np.array(np.abs(x) <= K*np.ones(len(x)))
            indices = np.where(truth == True)
            vxhatfunction = np.zeros(len(x), dtype=complex)
            for i in indices:
                vxhatfunction[i] = (constB(mode, W, K, R1)*sc.cosh(m0(W)*x[i]) + 
                                    constC(mode, W, K, R1)*sc.sinh(m0(W)*x[i]))
            truth2 = np.array(x < -K*np.ones(len(x)))
            indices2 = np.where(truth2 == True)
            for i in indices2:
                vxhatfunction[i] = constA(mode, W, K, R1)*(sc.cosh(m1(W, R1)*x[i]) + 
                                                          sc.sinh(m1(W, R1)*x[i]))
            truth3 = np.array(x > K*np.ones(len(x)))
            indices3 = np.where(truth3 == True)
            for i in indices3:
                vxhatfunction[i] = constD(mode, W, K, R1)*(sc.cosh(m2(W)*x[i]) - 
                                                          sc.sinh(m2(W)*x[i]))
    return vxhatfunction

def vyhat(mode, x, K):
    if mode in alfven_mode_options:
        vyhat_alfven = np.zeros_like(x)
        for i in range(len(x)):
            if x[i] > alfven_shear_width(mode, K)[0] and x[i] < alfven_shear_width(mode, K)[1]:
                vyhat_alfven[i] = alfven_amplitude
        return vyhat_alfven
    else:
        if type(x) == np.float64 or type(x) == float:
            return 0. + 0.j
        else:
            return np.zeros_like(x, dtype=complex)

def vzhat(mode, x, W, K, R1):
    if mode in alfven_mode_options:
        vzhat_function = np.zeros_like(x)
    else:
        if type(x) == float or type(x) == np.float64:
            if np.abs(x) <= K:
                vzhat_function = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*x) +
                                           constC(mode, W, K, R1)*sc.cosh(m0(W)*x))
            elif x < -K:
                vzhat_function = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*x) + 
                                                                     sc.cosh(m1(W, R1)*x))
            elif x > K:
                vzhat_function = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD(mode, W, K, R1)*(sc.sinh(m2(W)*x) -
                                                                 sc.cosh(m2(W)*x))
        else:
            truth = np.array(np.abs(x) <= K*np.ones(len(x)))
            indices = np.where(truth == True)
            vzhat_function = np.zeros(len(x), dtype=complex)
            for i in indices:
                vzhat_function[i] = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*x[i]) +
                                               constC(mode, W, K, R1)*sc.cosh(m0(W)*x[i]))
            truth2 = np.array(x < -K*np.ones(len(x)))
            indices2 = np.where(truth2 == True)
            for i in indices2:
                vzhat_function[i] = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) + 
                                                                         sc.cosh(m1(W, R1)*x[i]))
            truth3 = np.array(x > K*np.ones(len(x)))
            indices3 = np.where(truth3 == True)
            for i in indices3:
                vzhat_function[i] = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD(mode, W, K, R1)*(sc.sinh(m2(W)*x[i]) -
                                                                     sc.cosh(m2(W)*x[i]))
    return vzhat_function

# velocity
def vx(mode, x, z, t, W, K, R1):
    if type(vxhat(mode, x, W, K, R1)) == np.complex128:
        return vxhat(mode, x, W, K, R1) * np.exp(1j*(z-t))
    else:
        return np.outer(vxhat(mode, x, W, K, R1), np.exp(1j*(z-t)))

def vy(mode, x, z, t, K):
    return np.outer(vyhat(mode, x, K), np.exp(1j*(z-t)))

def vz(mode, x, z, t, W, K, R1):
    if type(vzhat(mode, x, W, K, R1)) == np.complex128:
        return vzhat(mode, x, W, K, R1) * np.exp(1j*(z-t))
    else:
        return np.outer(vzhat(mode, x, W, K, R1), np.exp(1j*(z-t)))

# velocity - lagrangian
def vx_pert(mode, x, z, t, W, K, R1):
    vx_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    vx_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix(mode, r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz(mode, r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                vx_hat_vals[i,j] = (constB(mode, W, K, R1)*sc.cosh(m0(W)*sol[0]) + 
                            constC(mode, W, K, R1)*sc.sinh(m0(W)*sol[0]))
            elif sol[0] < -K:
                vx_hat_vals[i,j] = constA(mode, W, K, R1)*(sc.cosh(m1(W, R1)*sol[0]) + 
                                                  sc.sinh(m1(W, R1)*sol[0]))
            elif sol[0] > K:
                vx_hat_vals[i,j] = constD(mode, W, K, R1)*(sc.cosh(m2(W)*sol[0]) - 
                                                  sc.sinh(m2(W)*sol[0]))
            vx_vals[i,j] = vx_hat_vals[i,j] * np.exp(1j*(z[j]-t))
    return vx_vals

def vz_pert(mode, x, z, t, W, K, R1):
    vz_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    vz_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix(mode, r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz(mode, r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                vz_hat_vals[i,j] = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*sol[0]) +
                                       constC(mode, W, K, R1)*sc.cosh(m0(W)*sol[0]))
            elif sol[0] < -K:
                vz_hat_vals[i,j] = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*sol[0]) + 
                                                                 sc.cosh(m1(W, R1)*sol[0]))
            elif sol[0] > K:
                vz_hat_vals[i,j] = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD(mode, W, K, R1)*(sc.sinh(m2(W)*sol[0]) -
                                                             sc.cosh(m2(W)*sol[0]))
            vz_vals[i,j] = vz_hat_vals[i,j] * np.exp(1j*(z[j]-t))
    return vz_vals    

# Displacement amplitude
def xix_hat(mode, x, W, K, R1):
    return vxhat(mode, x, W, K, R1) / W
    
def xiy_hat(mode, x, W, K):
    return vyhat(mode, x, K) / W
    
def xiz_hat(mode, x, W, K, R1):
    return vzhat(mode, x, W, K, R1) / W 
   
# Displacement
def xix(mode, x, z, t, W, K, R1):
    if type(vxhat(mode, x, W, K, R1)) == np.complex128:
        return (1j * vxhat(mode, x, W, K, R1) / W) * np.exp(1j*(z-t))
    else:
        return np.outer(1j * vxhat(mode, x, W, K, R1) / W, np.exp(1j*(z-t)))

def xiy(mode, x, z, t, W, K):
    if type(vyhat(mode, x, K)) == complex or type(vyhat(mode, x, K)) == np.complex128:
        return (1j * vyhat(mode, x, K) / W) * np.exp(1j*(z-t))
    else:
        return np.outer(1j * vyhat(mode, x, K) / W, np.exp(1j*(z-t)))

def xiz(mode, x, z, t, W, K, R1):
    if type(vzhat(mode, x, W, K, R1)) == np.complex128:
        return (1j * vzhat(mode, x, W, K, R1) / W) * np.exp(1j*(z-t))
    else:
        return np.outer(1j * vzhat(mode, x, W, K, R1) / W, np.exp(1j*(z-t)))

# Displacement amplitude at boundary
def xixhat_boundary(mode, W, K, R1, boundary='r'):
    if mode in alfven_mode_options:
        xixhat = 0.
    else:
        if boundary == 'r' or boundary == 'right':
            xixhat = (1j / W) * (constB(mode, W, K, R1)*sc.cosh(m0(W)*K) + 
                                     constC(mode, W, K, R1)*sc.sinh(m0(W)*K))
        if boundary == 'l' or boundary == 'left':
            xixhat = (1j / W) * (constB(mode, W, K, R1)*sc.cosh(m0(W)*-K) + 
                                     constC(mode, W, K, R1)*sc.sinh(m0(W)*-K))
    return xixhat

def xix_boundary(mode, z, t, W, K, R1, boundary='r'):           
    return xixhat_boundary(mode, W, K, R1, boundary) * np.exp(1j*(z-t))

# Mag field amplitude
def bxhat(mode, x, z, t, W, K, R1):
    if mode in alfven_mode_options:
        bxhat_function = np.zeros(len(x), dtype=complex)
    else:
        truth = np.array((x <= (K + xix_boundary(mode, z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                         (x >= (-K + xix_boundary(mode, z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
        indices = np.where(truth == True)
        bxhat_function = np.zeros(len(x), dtype=complex)
        for i in indices:
            bxhat_function[i] = (-B0/W)*(constB(mode, W, K, R1)*sc.cosh(m0(W)*x[i]) +
                                       constC(mode, W, K, R1)*sc.sinh(m0(W)*x[i]))
    return bxhat_function

def byhat(mode, x, K):
    if mode in alfven_mode_options:
        byhat_function = -(B0/vA) * vyhat(mode, x, K)
    else:
        byhat_function = np.zeros(len(x), dtype=complex)
    return byhat_function

def bzhat(mode, x, z, t, W, K, R1):
    if mode in alfven_mode_options:
        bzhat_function = np.zeros(len(x), dtype=complex)
    else:
        truth = np.array((x <= (K + xix_boundary(mode, z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                         (x >= (-K + xix_boundary(mode, z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
        indices = np.where(truth == True)
        bzhat_function = np.zeros(len(x), dtype=complex)
        for i in indices:
            bzhat_function[i] = (-1j*B0/W)*m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*x[i]) +
                                       constC(mode, W, K, R1)*sc.cosh(m0(W)*x[i]))
    return bzhat_function

# Mag field
def bx(mode, x, z, t, W, K, R1):
    bx_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bx_array[:,i] = bxhat(mode, x, z[i], t, W, K, R1) * np.exp(1j*(z[i]-t))
    return bx_array
    
def by(mode, x, z, t, K):
    by_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        by_array[:,i] = byhat(mode, x, K) * np.exp(1j*(z[i]-t))
    return by_array
    
def bz(mode, x, z, t, W, K, R1):
    bz_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bz_array[:,i] = bzhat(mode, x, z[i], t, W, K, R1) * np.exp(1j*(z[i]-t))
    return bz_array
    
# Density amplitude
def rho_hat(mode, x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    rho_hatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        rho_hatfunction[i] = m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*x[i]) +
                             constC(mode, W, K, R1)*sc.cosh(m0(W)*x[i])) * lamb00(W) / (c0**2 * m00(W))
#        if mode in slow_surf_mode_options + slow_body_1_mode_options + slow_body_2_mode_options:
#            rho_hatfunction[i] = -rho_hatfunction[i]
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        rho_hatfunction[i] = constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) +
                             sc.cosh(m1(W, R1)*x[i])) * lamb1(W, R1) / c1(R1)**2
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        rho_hatfunction[i] = constD(mode, W, K, R1)*(sc.sinh(m2(W)*x[i]) -
                             sc.cosh(m2(W)*x[i])) * lamb2(W) / c2**2
    return rho_hatfunction

# Density
def rho(mode, x, z, t, W, K, R1):
    rho_func = np.outer(rho_hat(mode, x, W, K, R1), np.exp(1j*(z-t)))
    return rho_func

# Denisty - lagrangian
def rho_pert(mode, x, z, t, W, K, R1):
    rho_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    rho_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix(mode, r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz(mode, r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                rho_hat_vals[i,j] = m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*sol[0]) +
                                constC(mode, W, K, R1)*sc.cosh(m0(W)*sol[0])) * lamb00(W) / (c0**2 * m00(W))
#                if mode in slow_surf_mode_options + slow_body_1_mode_options + slow_body_2_mode_options:
#                    rho_hat_vals[i,j] = -rho_hat_vals[i,j]
            elif sol[0] < -K:
                rho_hat_vals[i,j] = constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*sol[0]) +
                                sc.cosh(m1(W, R1)*sol[0])) * lamb1(W, R1) / c1(R1)**2
            elif sol[0] > K:
                rho_hat_vals[i,j] = constD(mode, W, K, R1)*(sc.sinh(m2(W)*sol[0]) -
                                     sc.cosh(m2(W)*sol[0])) * lamb2(W) / c2**2
            
            rho_vals[i,j] = rho_hat_vals[i,j] * np.exp(1j*(z[j]-t))
    return rho_vals
    
##############################################################################
# Alfven mixed driver functions    

# Frequency of the wave is dependant on x
def W_amd(x, vA_func):
    return vA_func(x)

def vxhat_amd(x):
    return np.zeros(len(x), dtype=complex)
    
def vyhat_amd(x, vA_func):
    return alfven_amplitude * vA_func(x)    

def vzhat_amd(x):
    return np.zeros(len(x), dtype=complex)
    
def vzhat_dash_amd(x):
    return np.zeros(len(x), dtype=complex)
    
def vx_amd(x, z, t, vA_func):
    vxhat_amd_md = vxhat_amd(x) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(vxhat_amd_md, np.exp(1j*z))

def vy_amd(x, z, t, vA_func):
    vyhat_amd_md = vyhat_amd(x, vA_func) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(vyhat_amd_md, np.exp(1j*z))

def vz_amd(x, z, t, vA_func):
    vzhat_amd_md = vzhat_amd(x) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(vzhat_amd_md, np.exp(1j*z))
    
def xix_hat_amd(x, vA_func):
    return vxhat_amd(x) / W_amd(x, vA_func)
    
def xiy_hat_amd(x, vA_func):
    return vyhat_amd(x, vA_func) / W_amd(x, vA_func)
    
def xiz_hat_amd(x, vA_func):
    return vzhat_amd(x) / W_amd(x, vA_func)

def xix_amd(x, z, t, vA_func):
    xix_hat_amd_md = xix_hat_amd(x, vA_func) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(1j*xix_hat_amd_md, np.exp(1j*z))

def xiy_amd(x, z, t, vA_func):
    xiy_hat_amd_md = xiy_hat_amd(x, vA_func) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(1j*xiy_hat_amd_md, np.exp(1j*z))

def xiz_amd(x, z, t, vA_func):
    xiz_hat_amd_md = xiz_hat_amd(x, vA_func) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(1j*xiz_hat_amd_md, np.exp(1j*z))
    
def bxhat_amd(x, vA_func):
    return -(B0/vA_func(x)) * vxhat_amd(x)

def byhat_amd(x, vA_func):
    return -(B0/vA_func(x)) * vyhat_amd(x, vA_func)

def bzhat_amd(x, vA_func):
    return -1j * (B0/W_amd(x, vA_func)) * vzhat_dash_amd(x)
    
def bx_amd(x, z, t, vA_func):
    bxhat_amd_md = bxhat_amd(x, vA_func) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(bxhat_amd_md, np.exp(1j*z))
    
def by_amd(x, z, t, vA_func):
    byhat_amd_md = byhat_amd(x, vA_func) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(byhat_amd_md, np.exp(1j*z))
    
def bz_amd(x, z, t, vA_func):
    bzhat_amd_md = bzhat_amd(x, vA_func) * np.exp(-1j*W_amd(x, vA_func)*t)
    return np.outer(bzhat_amd_md, np.exp(1j*z))

# equilibrium magnetic field is uniform
def bz_eq_amd(x, z):
    return B0 * np.ones((len(x), len(z)), dtype=complex)
    