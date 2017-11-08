"""
Dispersive elements used to separate light into its constituent wavelengths
(or frequencies, or energies) for recording a spectrum.
"""
import astropy.units as u
from astropy.units import (UnitConversionError, UnitsError, quantity_input,
                           Quantity)
from ...constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi)
import numpy as np


def bragg_angle(d_lattice, wavelength, order=1):
    r"""Calculate the thermal deBroglie wavelength for electrons.

    Parameters
    ----------
    d_lattice: Quantity
        Interlattice spacing of the Bragg diffraction crystal.
        
    wavelength: Quantity
        Wavelength of photon being diffracted.
    
    order: Quantity, int
        Diffraction order. Defaults to first order diffraction

    Returns
    -------
    theta: Quantity
        The Bragg angle at which the crystal diffracts the photon.

    Raises
    ------
    TypeError
        If argument is not a Quantity

    UnitConversionError
        If argument is in incorrect units

    ValueError
        If argument contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    The Bragg angle is the angle at which a photon of a given wavelength
    will be diffracted by a crystal with given interlattice spacing

    .. math::
    n\lambda = 2 d \sin(\theta)

    See also
    --------
    

    Example
    -------
    >>> from astropy import units as u
    >>> wavelength = 4.188655 * u.angstrom
    >>> d_lattice = 0.6708 / 2 * u.nm
    >>> bragg_angle(d_lattice, wavelength)
    <Quantity 0.6743974698099433 rad>
    

    """
    theta = np.arcsin(order * wavelength / (2 * d_lattice))
    return theta.to(u.rad)