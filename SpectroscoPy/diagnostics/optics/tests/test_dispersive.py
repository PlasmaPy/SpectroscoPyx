"""Tests for functions that uses Distribution functions."""

import numpy as np
import pytest
from astropy import units as u

from ....constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi)
from ..dispersive import (bragg_angle,
                          )

# test class for bragg_angle function:
class Test_bragg_angle(object):
    def setup_method(self):
        """initializing parameters for tests """
        self.d_lattice = 0.6708 / 2 * u.nm # HOPG spacing
        self.wavelength = 4.188655 * u.angstrom
    def test_typical_angle(self):
        """
        Checks whether function returns expected value for a typical
        crystal and photon configuration.
        """
        angleTrue = 0.6743974698099433 * u.rad
        angle = bragg_angle(d_lattice=self.d_lattice,
                            wavelength=self.wavelength,
                            order=1)
        errStr = (f"Bragg angle should be {angleTrue} "
                  f"and not {angle}.")
        assert np.isclose(angle.value, 
                          angleTrue.value,
                          rtol=1e-8,
                          atol=0.0), errStr