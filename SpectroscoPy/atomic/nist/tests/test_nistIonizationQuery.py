#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 11:49:30 2017

Unit tests for plasma_functions.py module.

@author: Pawel M. Kozlowski
"""

# importing Python modules
import pytest
import numpy as np
from astropy import units as u

# importing custom modules
from .. import nistIonizationQuery as niq


class Test_nistIonizationQuery(object):
    """Testing NIST ionization query functions"""
    def setup_method(self):
        """ Setting up for the test """
        ## Set element, charge state, and temperature
        self.element = 'Na'
        self.chargeState = 3
    def test_queryEnergy(self):
        """Test ionization energy query"""
        energy, uncert, dataType = niq.nistIonEnergyQuery(self.element,
                                                          self.chargeState)
        queryVal = np.array([energy.value,
                             uncert.value])
        queryExpect = np.array([9.8936e1,
                                1.2e-2])
        testTrue = np.allclose(queryVal, 
                               queryExpect,
                               rtol=1e-3,
                               atol=0)
        testType = dataType == "Interpolated/extrapolated"
        testTrue = testTrue and testType
        errStr = (f"nistIonEnergyQuery() method is wrong. Got {queryVal}, "
                  f"should be {queryExpect}.")
        assert testTrue, errStr