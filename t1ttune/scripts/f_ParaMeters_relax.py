import numpy as np
import scipy.constants as constants
from . import fun_hetrelax_models
import klassez as kz

def compute_taue(B, tauv, S=0.5, deltat=None, AMe=None, I=None, g=None):     
    """
    Compute the electron relaxation time taue using the Bloembergen-Morgan model for transient zero-field splitting (ZFS) contributions.
    
    Parameters
    ----------
    B: float
        magnetic field strength in Tesla
    S: float
        electron spin quantum number (default 0.5, one unpaired electron)
    tauv: float
        correlation time for transient coordination sphere fluctuations in seconds (default None, no transient ZFS)                 
    deltat: float or None
        transient ZFS parameter in cm^-1 (default None, no transient ZFS)
    AMe: float or None
        molecular alignment parameter (default None)
    I: float or None
        nuclear spin quantum number (default None)
    g: float or None
        electron g-factor (default None, uses free electron g-factor)

    Returns
    -------
    T1e: float
        longitudinal electron relaxation time in seconds
    T2e: float
        transverse electron relaxation time in seconds
    """

    if S == 0.5:
        T1e = 1./rotational_taue(g, B, tauv, AMe, I)[0]
        T2e = 1./rotational_taue(g, B, tauv, AMe, I)[1]
    else:
        T1e = 1./(transient_zfs(deltat, B, tauv, S)[0])
        T2e = 1./(transient_zfs(deltat, B, tauv, S)[1])
            
    return T1e, T2e    
            

  
def contactrelax(A, B, taue=1e-9, tau_M=np.inf, S=0.5, tauv=None, deltat=None, AMe=None, I=None, g=None):
    """
    Calculate the contact relaxation rates R1 and R2 using a simple contact relaxation model, with optional transient zero-field splitting (ZFS) contributions.
    
    Parameters
    ----------
    A: float
        contact coupling constant in Hz
    B: float
        magnetic field strength in Tesla
    taue: float
        electron relaxation time in seconds (default 1e-9 s)
    tau_M: float
        exchange correlation time in seconds (default np.inf)
    deltat: float or None
        transient ZFS parameter in cm^-1 (default None, no transient ZFS)
    tauv: float or None
        correlation time for transient ZFS fluctuations in seconds (default None, no transient ZFS)
    S: float
        electron spin quantum number (default 0.5, one unpaired electron)
        
    Returns
    -------
    R1_con: float
        longitudinal relaxation rate in s^-1
    R2_con: float
        transverse relaxation rate in s^-1
    """

    k_con = (5./3.)* S * (S+1) * (A)**2
    #Warning: The 2/5 factor is included in the spectral density function.    
    omegas = -B * constants.physical_constants['electron gyromag. ratio'][0]

 
 
    if tauv is not None:
        T1e, T2e = compute_taue(B, tauv, S=S, deltat=deltat, AMe=AMe, I=I, g=g)
    else:
        T1e = taue
        T2e = taue
    tau1 = 1./(1./T1e + 1./tau_M)
    tau2 = 1./(1./T2e + 1./tau_M)

    J_0_1 = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, 0, (tau1,))
    J_omegas_2 = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, omegas, (tau2,))
    
    
    R1_con = k_con * J_omegas_2
    R2_con = 0.5 * (k_con) * (J_omegas_2 + J_0_1)

    return R1_con, R2_con

def curie(B, r, S=0.5, nuc='1H', T=298, tau_r=1e-9, tau_M=np.inf, chi=None, sigma=None):
    r"""
    Curie-spin relaxation is implemented with two set of equations:
    - If the :math:`\chi` tensor is not provided the Curie-spin relaxation is computed using equations 4.30 and 4.31 in `Bertini et al. 2016`_, which assume an isotropic susceptibility.
    - If the :math:`\chi` tensor is provided, the Curie-spin relaxation is computed using equations 17-19 in `Suturina et al. 2018`_, which account for anisotropy in the susceptibility.
    
    .. _Bertini et al. 2016: https://www.sciencedirect.com/science/chapter/monograph/pii/B9780444634368000041
    .. _Suturina et al. 2018: https://pubs.rsc.org/en/content/articlelanding/2018/cp/c8cp01332b
    
    Parameters
    ----------
    B: float
        magnetic field strength in Tesla
    r: float or np.ndarray
        distance between the nucleus and the paramagnetic center in Angstroms as norm or as vector
    S: float
        electron spin quantum number (default 0.5, one unpaired electron)
    nuc: str
        nucleus type (default '1H')
    T: float
        temperature in Kelvin (default 298 K)
    chi: np.ndarray, optional
        susceptibility tensor (default None, isotropic)
    sigma: float, optional
        chemical shielding anisotropy (default None, not used)

    Returns
    -------
    R1_curie: float
        longitudinal relaxation rate in s^-1
    R2_curie: float
        transverse relaxation rate in s^-1
    """
    tau = 1./(1./tau_r + 1./tau_M)
    gamma_nuc = kz.sim.gamma[nuc]*(2*np.pi)*1e6 # in rad/s/T
    omegai = gamma_nuc * B
    J_0_tau = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, 0, (tau,))
    J_omegai_tau = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, omegai, (tau,))


    if chi is None:
        if isinstance(r, np.ndarray):
            r = np.linalg.norm(r)
        gamma_e = constants.physical_constants['electron gyromag. ratio'][0] # in rad/s/T
        k_B = constants.Boltzmann
        h = constants.Planck
        k_curie = (constants.mu_0/(4 * np.pi ))**2 * B**2 * gamma_nuc**2 * ((h * gamma_e/(2*np.pi))**4) * (S*(S+1.))**2 / ((3 * k_B * T)**2)
        r6 = r**6
        R1_curie = (k_curie/r6) * 3 * J_omegai_tau
        R2_curie = 0.5 * (k_curie/r6) * (3 * J_omegai_tau + 4 * J_0_tau)
        return R1_curie, R2_curie
    else:
        if isinstance(r, float) or isinstance(r, int):
            r_vec = np.array([r, 0, 0]) # assume r along x-axis if only distance is given    
        else:
            r_vec = r
            r = np.linalg.norm(r_vec)
        J_3omegai_tau = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, 3*omegai, (tau,))
        sigma_tot = sigma if sigma is not None else np.zeros(3)
        sigma_tot += (-1/(4*np.pi)) * (1/r**3) * chi @ (3 * np.outer(r_vec, r_vec)/r**2 - np.eye(3)) # equation 6 and 8 in Suturina et al. 2018
        bigLambdasquared = (sigma_tot[0,1] - sigma_tot[1,0])**2 + (sigma_tot[0,2] - sigma_tot[2,0])**2 + (sigma_tot[1,2] - sigma_tot[2,1])**2 # equation 18 in Suturina et al. 2018
        bigDeltasquared = sigma_tot[0,0]**2 + sigma_tot[1,1]**2 + sigma_tot[2,2]**2 - sigma_tot[0,0]*sigma_tot[1,1] - sigma_tot[0,0]*sigma_tot[2,2] - sigma_tot[1,1]*sigma_tot[2,2] + (3/4) * ((sigma_tot[0,1] + sigma_tot[1,0])**2 + (sigma_tot[0,2] + sigma_tot[2,0])**2 + (sigma_tot[1,2] + sigma_tot[2,1])**2) # equation 19 in Suturina et al. 2018
        R1_curie = (5./12.) * (omegai**2) * bigLambdasquared * J_3omegai_tau + (1./3.) * (omegai**2) * bigDeltasquared * J_omegai_tau
        R2_curie = (5./24.) * (omegai**2) * bigLambdasquared * J_3omegai_tau + (1./18.) * (omegai**2) * bigDeltasquared * (4 * J_0_tau + 3 * J_omegai_tau)
        return R1_curie, R2_curie

def SBM(r, B, taue=1e-9, taur=1e-9, tau_M=np.inf, S=0.5, tauv=None, deltat=None, AMe=None, I=None, g=None, nuc='1H'):
    """
    Calculate the dipolar relaxation rates R1 and R2 using the Solomon-Bloembergen-Morgan (SBM) equations, with optional transient zero-field splitting (ZFS) contributions.
    
    Parameters
    ----------
    
    r: float
        distance of between the nucleus and the paramagnetic center in Angstroms
    B: float
        magnetic field strength in Tesla
    g: float, or list of two floats, or 3x3 numpy array
        electron g-factor (default is the free electron g-factor)
    taue: float
        electron relaxation time in seconds (default 1e-9 s)
    taur: float
        rotational correlation time in seconds (default 1e-9 s)
    tau_M: float
        exchange correlation time in seconds (default np.inf)
    S: float
        electron spin quantum number (default 0.5, one unpaired electron)
    deltat: float or None
        transient ZFS parameter in cm^-1 (default None, no transient ZFS)
    tauv: float or None
        correlation time for transient ZFS fluctuations in seconds (default None, no transient ZFS)

    Returns
    -------
    R1_SBM: float
        longitudinal relaxation rate in s^-1
    R2_SBM: float
        transverse relaxation rate in s^-1
    """
    r6 = r**6
    gamma_nuc = kz.sim.gamma[nuc]*(2*np.pi)*1e6 # in rad/s/T
    if isinstance(g, (int, float, np.integer, np.floating)):
        pass
    elif isinstance(g, list) and len(g) == 2 and all(isinstance(x, (int, float, np.integer, np.floating)) for x in g):
    # Case 2: List of two elements (and both numeric)
        g = g[0] # use g_iso for the SBM calculation
    elif isinstance(g, np.ndarray) and g.shape == (3, 3) and np.issubdtype(g.dtype, np.number):
        gvals, gvecs = np.linalg.eig(g)
        g = np.mean(gvals) # use the g_iso for the SBM calculation
    elif g is None:
        g = constants.physical_constants['electron g factor'][0] # use free electron g-factor if g is not provided
    else:
        raise ValueError("g must be a float, a list of two floats, a 3x3 numpy array, or None")
        
    k_sbm = (1./3.) * (constants.mu_0/(4 * np.pi ))**2 * constants.physical_constants['Bohr magneton'][0]**2 * gamma_nuc**2 * g**2 *(S*(S+1.))
    #Warning: The 2/5 factor is included in the spectral density functions.
    omegai = gamma_nuc * B
    omegas = -B * constants.physical_constants['electron gyromag. ratio'][0]

 
    if deltat is not None:
        T1e, T2e = compute_taue(B, tauv, S=S, deltat=deltat, AMe=AMe, I=I, g=g)
    else:
        T1e = taue
        T2e = taue
    tau1 = 1./(1./T1e + 1./taur + 1./tau_M)
    tau2 = 1./(1./T2e + 1./taur + 1./tau_M)
    
    
    J_0_1 = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, 0, (tau1,))
    J_omegas_2 = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, omegas, (tau2,))
    J_omegai_1 = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, omegai, (tau1,))

    #print(f"J_0_1: {J_0_1}, J_omegas_2: {J_omegas_2}, J_omegai_1: {J_omegai_1}")

    R1_SBM = (k_sbm/r6)*(7. * J_omegas_2 + 3. * J_omegai_1)
    R2_SBM = 0.5 * (k_sbm/r6)*(13. * J_omegas_2 + 3. * J_omegai_1+4. * J_0_1)

    return R1_SBM, R2_SBM

def transient_zfs(deltat, B, tauv, S):
    """
    Calculate the transient zero-field splitting (ZFS) contributions to the electron relaxation rates T1e and T2e using the Bloembergen-Morgan model.
    
    Parameters
    ----------
    deltat: float
        transient ZFS parameter in cm^-1
    B: float
        magnetic field strength in Tesla
    tauv: float
        correlation time for transient ZFS fluctuations in seconds
    S: float
        electron spin quantum number
    
    Returns
    -------
    R1e: float
        electron longitudinal relaxation rate in s^-1
    R2e: float
        electron transverse relaxation rate in s^-1
    
    """

    deltatf = deltat*(2.998e10*2.*np.pi)
    
    omegas = -B * constants.physical_constants['electron gyromag. ratio'][0]
    J_omegas_v = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, omegas, (tauv,))
    J_2omegas_v = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, 2*omegas, (tauv,))
    J_0_v = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, 0, (tauv,))
    
    k = ((deltatf**2)/10.)*(4.*S*(S+1.)-3.)
    
    R1e = k * (J_omegas_v + 2* J_2omegas_v)    
    R2e = 0.5 * k * (3 * J_0_v + 5 * J_omegas_v + J_2omegas_v)

    return R1e, R2e

def rotational_taue(g, B, tauv, A=None, I=None):
    """
    Calculate the rotational contribution to the electron relaxation times.
    
    Parameters
    ----------
    g: float or array-like
        Electron g-factor or g-tensor
    B0: float
        Magnetic field strength in Tesla
    tauv: float
        Correlation time for rotational fluctuations in seconds
    A: float or array-like, optional
        Hyperfine coupling constant or tensor
    I: float or array-like, optional
        Nuclear spin quantum number or tensor
    
    Returns
    -------
    R1e: float
        Electron longitudinal relaxation rate in s^-1
    R2e: float
        Electron transverse relaxation rate in s^-1
    
    """

    ge = constants.physical_constants['electron g factor'][0]
    if isinstance(g, (int, float, np.integer, np.floating)):
        bigDeltag = np.abs(g - ge)
        smalldeltag = 0
    elif isinstance(g, list) and len(g) == 2 and all(isinstance(x, (int, float, np.integer, np.floating)) for x in g):
    # Case 2: List of two elements (and both numeric)
        bigDeltag = g[0]
        smalldeltag = g[1]
    elif isinstance(g, np.ndarray) and g.shape == (3, 3) and np.issubdtype(g.dtype, np.number):
    # Case 3: 3x3 NumPy array (numeric and correct shape)
        gvals, gvecs = np.linalg.eig(g)
    # find the entry in gvals that differs the most from the average of the other two
        avg_others = [(gvals[1] + gvals[2]) / 2, (gvals[0] + gvals[2]) / 2, (gvals[0] + gvals[1]) / 2]
        diff_others = [(gvals[1] - gvals[2]), (gvals[0] - gvals[2]), (gvals[0] - gvals[1])]
        distances = gvals - avg_others
        whichone = np.argmax(distances**2)
        bigDeltag = distances[whichone]
        smalldeltag = diff_others[whichone]/2
    else:
    #if none of the above, print an error message and exit
        raise ValueError("Invalid g input: must be a number, a list of two numbers, or a 3x3 numeric array.")

    omegas = B * constants.physical_constants['electron gyromag. ratio'][0]
    J_0_v = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, 0, (tauv,))
    J_omegas_v = fun_hetrelax_models.J(fun_hetrelax_models.J_omega_tau_iso, omegas, (tauv,))
    if A is None:
        ka1list = [0]
        ka2list = [0]
        kclist = [0]
    else:
        if I is None:
            raise ValueError("I missing, cannot compute contact relaxation without nuclear spin quantum number I.")
        else:
            if isinstance(A, (int, float, np.integer, np.floating)):
                bigDeltaA = A
                smalldeltaA = 0
            elif isinstance(A, list) and len(A) == 2 and all(isinstance(x, (int, float, np.integer, np.floating)) for x in A):
                bigDeltaA = A[0]
                smalldeltaA = A[1]
            elif isinstance(A, np.ndarray) and A.shape == (3, 3) and np.issubdtype(A.dtype, np.number):
                Avals, Avecs = np.linalg.eig(A)
                avg_others = [(Avals[1] + Avals[2]) / 2, (Avals[0] + Avals[2]) / 2, (Avals[0] + Avals[1]) / 2]
                diff_others = [(Avals[1] - Avals[2]), (Avals[0] - Avals[2]), (Avals[0] - Avals[1])]
                distances = Avals - avg_others
                whichone = np.argmax(distances**2)
                bigDeltaA = distances[whichone]
                smalldeltaA = diff_others[whichone]/2
            else:
                raise ValueError("Invalid A input: must be a number, a list of two numbers, or a 3x3 numeric array.")
            bigDeltaA*=constants.hbar
            smalldeltaA*=constants.hbar
            ka1list = [(1./8.)*(bigDeltaA**2 + smalldeltaA**2)*(7*I*(I+1)-mI**2) for mI in np.arange(-I, I+1)]
            ka2list = [(1./8.)*(bigDeltaA**2 + smalldeltaA**2)*(3*I*(I+1)-5*mI**2) for mI in np.arange(-I, I+1)]
            kclist = [(bigDeltag*bigDeltaA + 3*smalldeltag*smalldeltaA)*constants.physical_constants['Bohr magneton'][0]*B*mI for mI in np.arange(-I, I+1)]
 
    kg = (constants.physical_constants['Bohr magneton'][0]**2)*(B**2)*((bigDeltag**2)/3+smalldeltag**2)
    R1e = 0
    R2e = 0
    for mI, ka in enumerate(ka1list):
        k1 = (2./5.)*(1/constants.hbar**2)*(kg + ka - kclist[mI])
        k2 = (2./5.)*(1/constants.hbar**2)*((4./3.)*kg + ka2list[mI] -(4./3.)*kclist[mI])
        R1e += k1*J_omegas_v
        R2e += 0.5*(k2*J_0_v+k1*J_omegas_v)
        #print(k1, k2, kg, ka, kclist[mI])
    #print(1./R1e, 1./R2e)
    #exit()
    if len(ka1list)>1:
        R1e /= (2*I+1)
        R2e /= (2*I+1)

    return R1e,R2e

def OuterSphere(B, c=1, d=3.6e-10, D_target = 1e-10, D_cosolute = 2.6e-10, f=0.5, taue=1e-9, tauv=2.6e-11, deltat=0.014, AMe=None, I=None, g=None, S=3.5, nuc='1H'):
    """
    Calculate the outer sphere relaxation rates R1 and R2 using the Freed model, with optional transient zero-field splitting (ZFS) contributions.
    Default values are for 1 mM Gd-DOTA and a protein of 10 kDa at room temperature. Values for Gd-DOTA are taken from `Li et al. 2002`_.
    
    .. _Li et al. 2002: https://pubs.acs.org/doi/full/10.1021/ic0200390
    
    
    Parameters
    ----------
    B: float
        magnetic field strength in Tesla
    c: float
        concentration of the paramagnetic cosolute in mM (default 1 mM)
    d: float
        distance of closest approach between the nucleus and the paramagnetic center in Angstroms
    D_target: float
        diffusion coefficient of the target molecule in m^2/s (default 1e-10 m^2/s, corresponding to a protein of ~10 kDa at room temperature)
    D_cosolute: float
        diffusion coefficient of the paramagnetic cosolute in m^2/s (default 2.6e-10 m^2/s, corresponding to DOTA at room temperature)
    f: float
        fraction of the sphere of accessibility for the cosolute (default 0.5)
    taue: float
        electron relaxation time in seconds (default 1e-9 s)
    tauv: float or None
        correlation time for transient ZFS fluctuations in seconds (default 2.6e-11 s, for Gd^3+ complexes)
    deltat: float or None
        transient ZFS parameter in cm^-1 (default 0.014 cm^-1 for Gd^3+ complexes, None for no transient ZFS)
    AMe: float or None
        Hyperfine coupling constant to the metal center in Hz (default None, not used for Gd-DOTA)
    I: float or None
        nuclear spin quantum number (default None, not used for Gd-DOTA)
    g: float or None
        electron g-factor (default None, uses free electron g-factor, not used for Gd-DOTA)
    S: float
        electron spin quantum number (default 3.5, for Gd^3+)
    nuc: str
        nuclear spin quantum number (default '1H')

    Returns
    -------
    R1_outer: float
        longitudinal relaxation rate in s^-1
    R2_outer: float
        transverse relaxation rate in s^-1
    """
    
    h = constants.h
    N_A = constants.N_A
    mu_0 = constants.mu_0
    gamma_nuc = kz.sim.gamma[nuc]*(2*np.pi)*1e6 # in rad/s/T
    gamma_e = -1 * constants.physical_constants['electron gyromag. ratio'][0]
    k_outer = (16./81.) * np.pi * (mu_0/(4*np.pi))**2 * (h**2) * (gamma_nuc**2) * (gamma_e**2 / (4*np.pi**2)) * S*(S+1.) * N_A * c / (d * (D_target + D_cosolute))
    
    omegai = gamma_nuc * B
    omegas = gamma_e * B

 
    if tauv is not None:
        T1e, T2e = compute_taue(B, tauv, S=S, deltat=deltat, AMe=AMe, I=I, g=g)
    else:
        T1e = taue
        T2e = taue
    
    print("\nParameters for Outer Sphere Relaxation Calculation:")
    print("-" * 50)
    print("\n")
    print(f"Magnetic field strength (B): {B} T, corresponding to nuclear Larmor frequency of {omegai/(2*np.pi*1e6)} MHz for {nuc} and electron Larmor frequency of {omegas/(2*np.pi*1e9)} GHz.")
    print(f"\u03c4_e: {taue}")
    if tauv is not None:
        print("Bloembergen-Morgan pseudorotational contributions is used.") 
        print(f"\u03c4_v: {tauv}, S: {S}")
        print(f"\u0394_t: {deltat} cm^-1, A: {AMe}, I: {I}, g: {g}") 
    print(f"T1e: {T1e}, T2e: {T2e}")
    print("\n")
    print(f"Diffusion parameters:")
    print("-" * 20)
    print("\n")
    print(f"Concentration of the paramagnetic cosolute (c): {c} mM")
    print(f"Fraction of the sphere of accessibility for the cosolute (f): {f}")  
    print(f"distance of closest approach: {d*1e10} \u212B, D_target: {D_target} m^2/s, D_cosolute: {D_cosolute} m^2/s")
    
    J_0_outer = fun_hetrelax_models.J(fun_hetrelax_models.J_Freed, 0, {'d': d, 'D_target': D_target, 'D_cosolute': D_cosolute, 'tau_1': T1e, 'tau_2': T2e})
    J_omegai_outer = fun_hetrelax_models.J(fun_hetrelax_models.J_Freed, omegai, {'d': d, 'D_target': D_target, 'D_cosolute': D_cosolute, 'tau_1': T1e, 'tau_2': T2e})    
    J_omegas_outer = fun_hetrelax_models.J(fun_hetrelax_models.J_Freed, omegas, {'d': d, 'D_target': D_target, 'D_cosolute': D_cosolute, 'tau_1': T1e, 'tau_2': T2e})
    
    
    R1_outer = k_outer * (7 * J_omegas_outer + 3 * J_omegai_outer)
    R2_outer = 0.5 * k_outer * (13 * J_omegas_outer + 3 * J_omegai_outer + 4 * J_0_outer)
    
    return R1_outer, R2_outer

if __name__ == "__main__":
    B = 14.1
    r = 2.5e-10
    A = 1e6

   
    R1_SBM, R2_SBM = SBM(r, B)
    R1_curie, R2_curie = curie(B, r)
    R1_con, R2_con = contactrelax(A, B)
    R1_outer, R2_outer = OuterSphere(B)
    
    print(R1_SBM, R2_SBM)
    print(R1_curie, R2_curie)
    print(R1_con, R2_con)
    print(R1_outer, R2_outer)