#!/usr/bin/env python3


def water_viscosity(T):
    """
    Returns the viscosity of water in centipoise at temperature T (Kelvin)
    Equation from CRC Handbook of Chemistry and Physics, 96th Edition, 2015
    
    Parameters:
    -----------
    T : float
        Temperature in Kelvin
    
    Returns:
    --------
    eta_water : float
        Viscosity of water in centipoise
    """
    
    eta_water = 1.7753 - 0.0565 * (T - 273) + 1.0751e-3 * (T - 273)**2 - 9.222e-6 * (T - 273)**3
    return eta_water

def recompute_tau(tau, T_new, T_old=298.15):
    """
    Recompute the correlation time tau at a new temperature T_new, given the old temperature T_old and the old correlation time tau. Assumes that the correlation time is proportional to the viscosity of water at the given temperature, and inversely proportional to temperature, according to the Stokes-Einstein relation.
    
    Parameters:
    -----------
    tau : float
        Correlation time at the old temperature T_old, in seconds
    T_new : float
        New temperature in Kelvin
    T_old : float, optional
        Old temperature in Kelvin (default is 298.15 K)
    
    Returns:
    --------
    tau_new : float
        Recomputed correlation time at the new temperature T_new, in seconds
    """
    
    eta_old = water_viscosity(T_old)
    eta_new = water_viscosity(T_new)
    tau_new = tau * (eta_new / eta_old) * (T_old / T_new)
    return tau_new
