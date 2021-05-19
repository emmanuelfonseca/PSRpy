from ..const import T_sun
import numpy as np

def mass_function(period: float, axis_semimajor_projected: float):
    """
    Computes Keplerian mass function, given projected size and orbital period.

    Parameters
    ----------
    period : float 
        period of the orbit, in units of days

    axis_semimajor_projected : float 
        component of Keplerian semimajor axis, projected onto the line of sight, 
        in units of light-seconds

    Returns
    -------
    mass_function : float 
        Keplerian mass function, in units of solar masses
    """

    frequency = 2 * np.pi / period / 86400
    mass_function = frequency**2 * axis_semimajor_projected**3 / T_sun.value

    return mass_function

def mass_companion(period: float, axis_semimajor_projected: float, mass_pulsar: float, 
    inclination: float, initial_guess: float = 0.5, nr_attempts: int = 100, 
    nr_tolerance: float = 1e-12):
    """
    Computes the companion mass from the Keplerian mass function, assuming values of 
    the pulsar mass and system inclination. This function uses a Newton-Raphson method as 
    the equation, rewritten using the mass function, is transcendental in companion mass.

    Parameters
    ----------
    period : float 
        period of the orbit, in units of days

    axis_semimajor_projected : float 
        component of Keplerian semimajor axis, projected onto the line of sight, 
        in units of light-seconds

    mass_pulsar : float
        mass of pulsar, in units of solar masses

    inclination : float 
        geometric inclination of orbital system, in units of degrees

    initial_guess : float, optional 
        initial guess of companion mass for NR method.        

    nr_attempts : int, optional
        number of iterations for attempting the Newton-Raphson method 

    nr_tolerance : float, optional
        tolerance used for evaluating accuracy of the Newton-Raphson method

    Returns
    -------
    mass_function : float 
        Keplerian mass function, in units of solar masses
    """

    # first, compute the mass function.
    mf = mass_function(period, axis_semimajor_projected)
    si = np.sin(inclination * np.pi / 180)

    # define variables needed for NR method.
    mass_companion = 0.
    mc_mid = initial_guess
    mc0 = initial_guess

    # loop over desired number of attempts to perform NR calculation.
    for idx in range(nr_attempts):
        mass_total = mass_pulsar + mc_mid
        func = (mc_mid * si)**3 / mass_total**2 - mf
        func_deriv = mc_mid**2 * si**3 * (mc_mid + 3 * mass_pulsar) / mass_total**3
        mc_mid -= func / func_deriv

        # if the mass value has sub-threshold accuracy, break out.
        if np.fabs(mc_mid - mc0) < nr_tolerance:
            mass_companion = mc_mid
            break

        # otherwise, update prior-step array and redo calculation.
        else:
            mc0 = mc_mid

    return mass_companion

def mass_pulsar(period: float, axis_semimajor_projected: float, mass_companion: float, 
    inclination: float):
    """
    Computes the pulsar mass from the Keplerian mass function, using radio-timing data 
    and assuming values of the companion mass and system inclination.

    Parameters
    ----------
    period : float 
        period of the orbit, in units of days

    axis_semimajor_projected : float 
        component of Keplerian semimajor axis, projected onto the line of sight, 
        in units of light-seconds

    mass_companion : float
        mass of pulsar, in units of solar masses

    Returns
    -------
    mass_pulsar : float 
        Keplerian mass function, in units of solar masses
    """

    # first, compute the mass function and sine of inclination angle.
    mf = mass_function(period, axis_semimajor_projected)
    si = np.sin(inclination * np.pi / 180)

    # now use mass function to compute pulsar mass.
    mass_pulsar = np.sqrt((mass_companion * si)**3 / mf) - mass_companion

    return mass_pulsar

def mass_total_from_omdot(period: float, eccentricity: float, omdot: float, omdot_error: float = None):
    """
    Computes the total mass of the pulsar-binary system from the Keplerian elements and 
    the observed periastron advance as predicted by general relativity.

    Parameters
    ----------
    period : float 
        period of the orbit, in units of days

    eccentricity : float
        Keplerian eccentricity for closed orbits (i.e., 0 <= eccentricity < 1)

    omdot : float
        rate of change in argument of periastron, in units of degrees per year

    omdot_error : float, optional
        uncertainty in omdot, in units of degrees per year

    Returns
    -------
    mass_total_list : array_like
        a two-element list containing the total mass and its corresponding uncertainty; 
        the uncertainty is set to None if uncertainty in omdot is not provided
    """

    # compute total mass from equation for relativistic advance of periastron.
    pb_in = period * 86400
    omdot_in = omdot * np.pi / 180 / 365.25 / 86400
    term_ecc = 1 - eccentricity**2
    mass_total = (omdot_in / 3 * (pb_in / 2 / np.pi)**(5./3.) * term_ecc)**(1.5) / T_sun.value

    # if an uncertainty is provided, then compute propagated uncertainty in total mass.
    mass_total_error = None

    if omdot_error is not None:
        omdot_error_in = omdot_error * np.pi / 180 / 365.25 / 86400
        mass_total_error = ((pb_in / 2 / np.pi)**(5./3.) / 3 * term_ecc)**(1.5) * \
                           (1.5 * np.sqrt(omdot_in)) * omdot_error_in / T_sun.value

    # load data to return as a two-element list.
    mass_total_list = np.array([mass_total, mass_total_error])
    
    return mass_total_list

