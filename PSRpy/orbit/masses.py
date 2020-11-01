from ..const import T_sun
import numpy as np

def mass_function(pb, x):
    """
    Computes Keplerian mass function, given projected size and orbital period.

    Inputs:
        - pb = orbital period [days]
        - x = projected semimajor axis [lt-s]

    Output:
        - mass function [solar mass]
    """

    nb = 2 * np.pi / pb / 86400
    return nb**2 * x**3 / T_sun

def mass_companion(pb, x, mp, sini, mc=0.5, tolerance=1e-12):
    """
    Computes the companion mass from the Keplerian mass function. This function 
    uses a Newton-Raphson method since the equation is transcendental.
    """

    # first, compute the mass function.
    mf = mass_function(pb, x)

    # use a Newton-Raphson method for determining the companion mass 
    # from the mass function, for an arbitrary value of sini.
    mc_current = mc
    mc_before = mc

    for ii in range(100):
        g = (mc_current * sini)**3 / (mp + mc_current)**2 - mf
        dgdmc = mc_current**2 * sini**3 * (mc_current + 3 * mp) / (mp + mc_current)**3
        mc_current -= (g / dgdmc)

        if (np.fabs(mc_current - mc_before) < tolerance):
            break

        mc_before = mc_current

    return mc_current

def mass_pulsar(pb, x, mc, sini):
    """
    Computes the companion mass from the Keplerian mass function. This function 
    uses a Newton-Raphson method since the equation is transcendental.
    """

    mf = mass_function(pb, x)
    return np.sqrt((mc * sini)**3 / mf) - mc

def mass_total(pb, ecc, omdot, omdot_error=None):
    """
    Computes the total mass of the pulsar-binary system from the Keplerian elements and 
    the observed periastron advance as predicted by general relativity.
    """

    """
    Calculate the total mass of the binary system, given a measurement of OMDOT.
    """

    pb_in = pb * 86400
    omdot_in = omdot * np.pi / 180 / 365.25 / 86400
    total_mass = (omdot_in / 3 * (pb_in / 2 / np.pi)**(5./3.) * (1 - ecc**2))**(1.5) / T_sun

    if (omdot_error is not None):
        omdot_error_in = omdot_error * np.pi / 180 / 365.25 / 86400
        total_mass_error = ((pb_in / 2 / np.pi)**(5./3.) / 3 * (1 - ecc**2))**(1.5) * \
                           (1.5 * np.sqrt(omdot_in)) * omdot_error_in / T_sun
        return (total_mass, total_mass_error)

    else:
        return total_mass
