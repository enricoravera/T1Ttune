#! /usr/bin/env python3

import numpy as np
import klassez as kz
import scipy.constants


def d(r=1.02e-10, nuc1='1H', nuc2='15N'):
    """
    Function to calculate dipole-dipole coupling constant d
    
    Parameters:
    -----------
    r : float
        distance between the two nuclei in meters (default is 1.02e-10 m, typical for NH bond)
    nuc1 : str
        first nucleus (default is '1H')
    nuc2 : str
        second nucleus (default is '15N')
    
    Returns:
    --------
    d_val : float
        dipole-dipole coupling constant in rad/s
        
    """
    gamma1 = kz.sim.gamma[nuc1]*2*np.pi*1e6
    gamma2 = kz.sim.gamma[nuc2]*2*np.pi*1e6
    d_val = scipy.constants.mu_0 * scipy.constants.h * gamma1 * gamma2 / (r)**3 / (16 * np.pi**2)
    #print("d value: ", d_val / (2*np.pi*1e3), " kHz")
    return d_val

def c(Deltasigma=-160, nuc='15N'):
    """
    function to calculate the CSA constant c
    
    Parameters:
    -----------
    Deltasigma: float
        CSA value in ppm (default is -160)
    nuc: str
        nucleus (default is '15N')
    
    Returns:
    --------
    c_val: float
        CSA constant in rad T^-1 s^-1
    """
    gamma = kz.sim.gamma[nuc]*2*np.pi*1e6 # in rad T^-1 s^-1
    ds = Deltasigma * 1e-6 # dimensionless
    c_val = -(ds * gamma) / 3# in rad T^-1 s^-1
    #print("c value: ", c_val / (2*np.pi*1e3), " kHz T^-1")
    return c_val

def omega(B, nuc='1H'):
    """
    Function to calculate the Larmor frequency omega
    
    Parameters:
    -----------
    B: float
        magnetic field strength in Tesla
    nuc: str, optional
        nucleus (default is '1H')
    
    Returns:
    --------
    w: float
        Larmor frequency in rad/s
    """
    gamma = kz.sim.gamma[nuc]*2*np.pi*1e6
    w = gamma * B
    #print("Larmor frequency for ", nuc, " at ", B, " T: ", w / (2*np.pi*1e6), " MHz")
    return w

def J_omega_tau_iso(w, tau):
    r"""
    Function to calculate the spectral density J for isotropic tumbling, given the Larmor frequency w and the correlation time :math:`\tau`
    
    Parameters:
    -----------
    w: float
        Larmor frequency in rad/s
    tau: float
        correlation time in seconds
    
    Returns:
    --------
    J: float
        isotropic spectral density
    """
    
    J = (2/5) * (tau / (1 + (w * tau)**2))
    return J

def J_omega_tau_aniso(w, D_comp, r_pdb):
    """
    Function to calculate the anisotropic spectral density J, given the tensor D and the orientation r
    
    Parameters:
    -----------
    w: float
        Larmor frequency in rad/s
    D_comp: array_like
        components of the diffusion tensor (Dxx, Dyy, Dzz, Dxy, Dxz, Dyz)
    r_pdb: array_like
        unit vector along the NH bond in the PDB frame
    
    Returns:
    --------
    J: float
        anisotropic spectral density
    """
    r_pdb /= np.linalg.norm(r_pdb)
    D_mat = [[D_comp[0], D_comp[3], D_comp[4]],
             [D_comp[3], D_comp[1], D_comp[5]],
             [D_comp[4], D_comp[5], D_comp[2]]]
    Deig, R = np.linalg.eig(D_mat)
    r = R.T @ r_pdb
    Diso = (Deig[0] + Deig[1] + Deig[2]) / 3
    Dq = (Deig[0] * Deig[1] + Deig[0] * Deig[2] + Deig[1] * Deig[2]) / 3
    Dd = 6*np.sqrt(Diso**2-Dq)
    D = np.zeros(5)
    A = np.zeros(5)
    D[0] = 4 * Deig[0] + Deig[1] + Deig[2]
    D[1] = 4 * Deig[1] + Deig[0] + Deig[2]
    D[2] = 4 * Deig[2] + Deig[0] + Deig[1]
    D[3] = 6 * Diso + Dd
    D[4] = 6 * Diso - Dd
    A[0] = 3 * r[1]**2 * r[2]**2
    A[1] = 3 * r[0]**2 * r[2]**2
    A[2] = 3 * r[0]**2 * r[1]**2
    d0 = (D[0] - Diso) / np.sqrt(Diso**2 - Dq)
    d1 = (D[1] - Diso) / np.sqrt(Diso**2 - Dq)
    d2 = (D[2] - Diso) / np.sqrt(Diso**2 - Dq)

    A[3] = (1/4) * (3*(np.sum(r**4)) - 1) - (1/12) * (d0 * (3 * r[0]**4 + 6 * r[1]**2 * r[2]**2 - 1) + d1 * (3 * r[1]**4 + 6 * r[0]**2 * r[2]**2 - 1) + d2 * (3 * r[2]**4 + 6 * r[0]**2 * r[1]**2 - 1))
    A[4] = (1/4) * (3*(np.sum(r**4)) - 1) + (1/12) * (d0 * (3 * r[0]**4 + 6 * r[1]**2 * r[2]**2 - 1) + d1 * (3 * r[1]**4 + 6 * r[0]**2 * r[2]**2 - 1) + d2 * (3 * r[2]**4 + 6 * r[0]**2 * r[1]**2 - 1))

    Dlor = D / (D**2 + np.ones(5)*w**2)
    J = (2/5) * np.sum(A * Dlor)
 
    return J

def LS_iso(w, S2s, taus):
    """
    Isotropic reorientation extended model-free Lipari-Szabo density function, given the Larmor frequency w, the order parameters S2s and the correlation times taus. S2s should be one less than taus, as the last term corresponds to the global tumbling. The function calculates the weights for each term based on the S2s and then sums up the contributions from each term using the J_omega_tau_iso function.
    
    Parameters:
    -----------
    w: float
        Larmor frequency in rad/s
    S2s: single value or list of floats
        order parameters one less than the taus
    taus: list of floats
        correlation times
    
    Returns:
    --------
    LS: float
        Lipari-Szabo spectral density
    """
    if isinstance(S2s, (float, int)):
        S2s = [S2s]
    if len(S2s) != len(taus)-1:
        raise ValueError("S2s must be one less than the taus ")
    weights = np.zeros(len(taus))
    prod = 1.0
    for i, s in enumerate(S2s):
        weights[i] = s * prod
        prod *= (1 - s)
    weights[-1] = prod
    LS = 0
    for i in range(len(taus)):
        LS += weights[i] * J_omega_tau_iso(w, taus[i])
    return LS

def LS_aniso(w, S2s, taus, D_comp, r_pdb):
    """
    Anisotropic reorientation extended model-free Lipari-Szabo density function, given the Larmor frequency w, the order parameters S2s, the correlation times taus, the diffusion tensor components D_comp and the orientation r_pdb. S2s should be as long as taus, as the first term corresponds to the anisotropic contribution and the rest to the isotropic contributions. The function calculates the weights for each term based on the S2s and then sums up the contributions from each term using the J_omega_tau_aniso and J_omega_tau_iso functions.
    
    Parameters:
    -----------
    w: float
        Larmor frequency in rad/s
    S2s: single value or list of floats
        order parameters as long as the taus
    taus: list of floats
        correlation times
    D_comp: array_like
        diffusion tensor components
    r_pdb: array_like
        orientation of the NH bond in the PDB frame
    
    Returns:
    --------
    LS: float
        Lipari-Szabo spectral density
    """
    if isinstance(S2s, (float, int)):
        S2s = [S2s]
    if len(S2s) != len(taus):
        raise ValueError("S2s must be as long as taus ")
    weights = np.zeros(len(taus)+1)
    prod = 1.0
    for i, s in enumerate(S2s):
        weights[i] = s * prod
        prod *= (1 - s)
    weights[-1] = prod
    LS = weights[0] * J_omega_tau_aniso(w, D_comp, r_pdb)
    for i in range(len(taus)):
        LS += weights[i] * J_omega_tau_iso(w, taus[i])
    return LS

def J_Freed(w, d, D_target, D_cosolute, tau_1=1e-9, tau_2=None):
    """
    Freed spectral density function for outer sphere relaxation. Equations 6.42, 6.48 and 6.50 in `Bertini et al. 2016`_.
 
    .. _Bertini et al. 2016: https://www.sciencedirect.com/science/chapter/monograph/pii/B9780444634368000065

    Parameters
    -----------
    w: float
        Larmor frequency in rad/s
    d: float
        distance of closest approach of the paramagnetic center and the nucleus in meters
    D_target: float
        diffusion coefficient of the target molecule in m^2/s
    D_cosolute: float
        diffusion coefficient of the cosolute in m^2/s
    tau_1: float
        longitudinal correlation time of the electron in seconds
    tau_2: float or None
        transverse correlation time of the electron in seconds (default None, will be set to tau_1 if not provided)

    Returns
    --------
    J: float
        spectral density
    """
 
    tau_D = d**2 / (D_target + D_cosolute) # equation 6.42
    if tau_2 is None:
        tau_2 = tau_1
    z = np.sqrt(2 * np.abs(w) * tau_D + tau_D / tau_1) # equation 6.50
    
    J = (2/5) * (1 + 5 * z / 8 + z**2 / 8) / (1 + z + z**2 / 2 + z**3 / 6 + 4 * z**4 / 81 + z**5 / 81 + z**6 / 648)  #equation 6.48

    return J

def J(J_func, w, f_args):
    """
    General function to calculate the spectral density J using a given J_func. Supports:
    - Isotropic Lipari-Szabo: J_func = :func:`LS_iso`, f_args = (S2s, taus)
    - Anisotropic Lipari-Szabo: J_func = :func:`LS_aniso`, f_args = (S2s, taus, D_comp, r_pdb)
    - Freed outer sphere: J_func = :func:`J_Freed`, f_args = (d, D_target, D_cosolute, tau_1, tau_2)
        
    Parameters:
    -----------
    J_func: function 
        spectral density function to use (:func:`LS_iso` or :func:`LS_aniso`)
    w: float
        Larmor frequency in rad/s
    f_args: list or dict
        arguments for the J_func
    
    Returns
    --------
    J: float
        spectral density
    """
    if isinstance(f_args, dict):
        return J_func(w, **f_args)

    return J_func(w, *f_args)

# H-X relaxation terms R1, R2, NOE. everything defaults to NH. 
def R1R2nOe(B, r=1.02e-10, nuc1='1H', nuc2='15N', Deltasigma=-160, func=LS_iso, f_args=([0.9], [1e-9, 1e-11]), Rex=0):
    r"""
    H-X relaxation rates R1, R2 and NOE enhancement factor eta, given the magnetic field strength B, the distance r between the nuclei, the types of nuclei nuc1 and nuc2, the CSA value Deltasigma, the function to calculate the spectral density func and its arguments f_args, and an optional Rex contribution to R2. The function calculates the dipole-dipole coupling constant d and the CSA constant c, then uses them to calculate R1, R2 and eta based on the formulas reported by `Fushman`_:
    
    .. _Fushman: https://pmc.ncbi.nlm.nih.gov/articles/PMC4361738/
    
    .. math:: R_1 = d^2 (J(\omega_H - \omega_X) + 6 J(\omega_H + \omega_X)) + 3 (c^2 B^2 + d^2) J(\omega_X)
    
    .. math:: R_2 = R_{ex} + \frac{1}{2} d^2 (J(\omega_H - \omega_X) + 6 J(\omega_H + \omega_X) + 6 J(\omega_H)) + \frac{1}{2} (c^2 B^2 + d^2) (4 J(0) + 3 J(\omega_X))
    
    .. math:: \eta = 1 - \frac{d^2 \gamma_1}{\gamma_2} \frac{6 J(\omega_H + \omega_X) - J(\omega_H - \omega_X)}{R_1}
    
    Parameters:
    -----------
    B: float
        magnetic field strength in Tesla
    r: float
        distance between the nuclei in meters (default is 1.02e-10   m, typical for NH bond)
    nuc1: str
        first nucleus (default is '1H')
    nuc2: str
        second nucleus (default is '15N')
    Deltasigma: float
        CSA value in ppm (default is -160)
    func: function
        function to calculate the spectral density (default is LS_iso)
    f_args: list
        arguments for the spectral density function (default is ([0.9], [1e-9, 1e-11]))
    Rex: float
        optional exchange contribution to R2 (default is 0)
    
    Returns:
    --------
    R1: float
        longitudinal relaxation rate in s^-1
    R2: float
        transverse relaxation rate in s^-1
    eta: float
        NOE enhancement factor (dimensionless)
    """
    d_val = d(r=r, nuc1=nuc1, nuc2=nuc2)
    c_val = c(Deltasigma=Deltasigma, nuc=nuc2)
    gamma1overgamma2 = np.abs(kz.sim.gamma[nuc1] / kz.sim.gamma[nuc2])
    J_0 = J(func, 0, f_args)
    J_nuc1 = J(func, omega(B, nuc1), f_args)
    J_nuc2 = J(func, omega(B, nuc2), f_args)
    J_sum = J(func, omega(B, nuc1) + omega(B, nuc2), f_args)
    J_diff = J(func, abs(omega(B, nuc1) - omega(B, nuc2)), f_args)
    R1 = (d_val**2 * (J_diff + 6 * J_sum)) + 3 * (c_val**2 * B**2 + d_val**2) * J_nuc2
    R2 = Rex + 0.5 * d_val**2 * (J_diff + 6 * J_sum + 6 * J_nuc1) + 0.5 * (c_val**2 * B**2 + d_val**2) * (4 * J_0 + 3 * J_nuc2)
    eta = 1 - (d_val**2 / R1) * gamma1overgamma2 * (6 * J_sum - J_diff)
    return R1, R2, eta

def eta_z_eta_xy(B, r=1.02e-10, nuc1='1H', nuc2='15N', Deltasigma=-160, theta=17*np.pi/180, func=LS_iso, f_args=([0.9], [1e-9, 1e-11])):
    r"""
    Compute cross-correlated relaxation rates eta_z and eta_xy for a given magnetic field strength B, distance r between the nuclei, types of nuclei nuc1 and nuc2, CSA value Deltasigma, angle theta between the NH bond and the magnetic field, function to calculate the spectral density func and its arguments f_args. The function calculates the dipole-dipole coupling constant d and the CSA constant c, then uses them to calculate eta_z and eta_xy based on the formulas reported by `Salvi`_:
    
    .. _Salvi: https://www.sciencedirect.com/science/article/pii/S0079656517300213
    
    .. math:: \eta_z = \frac{1}{15} P_2(\cos\theta) d c B 6 J(\omega_X)
    
    .. math:: \eta_{xy} = \frac{1}{15} P_2(\cos\theta) d c B (4 J(0) + 3 J(\omega_X))
    
    Parameters:
    -----------
    B: float
        magnetic field strength in Tesla
    r: float
        distance between the nuclei in meters (default is 1.02e-10 m, typical for NH bond)
    nuc1: str
        first nucleus (default is '1H')
    nuc2: str
        second nucleus (default is '15N')
    Deltasigma: float
        CSA value in ppm (default is -160)
    theta: float
        angle between the NH bond and the magnetic field in radians (default is 17 degrees)
    func: function
        function to calculate the spectral density (default is :func:`LS_iso`)
    f_args: list
        arguments for the spectral density function (default is ``([0.9], [1e-9, 1e-11])``)
    
    Returns:
    --------
    eta_z and eta_xy: float
        cross-correlated relaxation rates in s^-1
    
    """
    d_val = d(r=r, nuc1=nuc1, nuc2=nuc2)
    c_val = c(Deltasigma=Deltasigma, nuc=nuc2)
    P2_costheta = 0.5 * (3 * np.cos(theta)**2 - 1)
    J_0 = J(func, 0, f_args)
    J_nuc2 = J(func, omega(B, nuc2), f_args)
    eta_z = P2_costheta * d_val * c_val * B * 6 * J_nuc2
    eta_xy = P2_costheta * d_val * c_val * B * (4 * J_0 + 3 * J_nuc2)
    return eta_z, eta_xy 
