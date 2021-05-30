#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 15:18:45 2017

Queries NIST database on ionization energies.
Mainly used for generating Saha-Boltzmann plots.

TODO:
- Add tests for nistIonEnergiesQuery()
- Move status of connection into a log file
-improve this query function by making it a class and allowing for the use
    of methods to get different column values
- Use code for generalized query of all energies for all charge states
    as this will be a useful speedup
- implement an object that can be returned from this query and methods
  that take slices through the table

@author: Pawel M. Kozlowski
"""

# importing Python modules
import requests as req
import parse as parse
import re
import numpy as np
import astropy.units as u
import mendeleev as lev


def cleanhtml(raw_html):
    """
    Cleaning HTML tags using regular expressions
    """
    cleanr = re.compile('<.*?>')
    cleantext = re.sub(cleanr, '', raw_html)
    return cleantext


# URL to nist ionization energy database
nistIonizationURL1 = 'https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html'
nistIonizationURL2 = 'https://physics.nist.gov/cgi-bin/ASD/ie.pl'

# Checking website status
reqObj = req.post(nistIonizationURL1)
print("STATUS:")
print(reqObj.status_code)
if reqObj.status_code == req.codes.ok:
    print("OK!")
else:
    print("CONNECTION FAILED!")
reqObj.raise_for_status()

# for testing purposes only
# print("HEADERS:")
# print(reqObj.headers['content-type'])
# print("ENCODING:")
# print(reqObj.encoding)
# print("TEXT:")
# print(reqObj.text)

#%% NIST levels


def nistIonEnergiesQuery(element):
    """
    Given an element, and a charge state this
    script queries the NIST ionization energy level database and retrieves 
    the value of ionization energy and its uncertainty. The uncertainty
    is given in units of the last decimal place of the value.
    Elements (string)
    chargeState [dimless integer]
    Ionization energy [eV]
    Uncertainty [eV]
    """
    # checking inputs
    if not isinstance(element, str):
        raise ValueError("Element must be a string!")
    elementCharge = element
    # Starting requests session
    with req.session() as s:
        # values to be filled into website form
        payload = {'spectra': elementCharge,  # element and charge state
                   'submit': 'Retrieve Data',
                   'units': '1',  # selects energy units eV
                   'format': '1',  # ASCII output
                   'order': '0',  # order by Z or sequence
                   'at_num_out': 'on',  # output atomic number
                   'sp_name_out': 'on',  # output spectrum name
                   'ion_charge_out': 'on',  # output ion charge
                   'el_name_out': 'on',  # output element name
                   'seq_out': 'on',  # output isoelectronic sequence
                   'shells_out': 'on',  # output ground-state electronic shells
                   'conf_out': 'on',  # outpuse ground-state configuration
                   'level_out': 'on',  # output ground-state level
                   'ion_conf_out': 'on',  # output ionization configuration
                   'e_out': '0',  # selects ionization energy not binding
                   'unc_out': 'on',  # output uncertainty
                   'biblio': 'on'  # output bibliographic references
                   }
        # website's repsonse to query
        response = s.post(nistIonizationURL2, data=payload)
#        # for testing purposes only
#        print("URL:")
#        print(response.url)
#        print("STATUS:")
#        print(response.status_code)
#        print("TEXT:")
#        print(response.text)

    # Extract ionization energy value string from HTML/ASCII
    # Getting entire ASCII table
#    searchStringTable = "<pre>{}</pre>"
        searchStringTable = "<pre>{}<p"
    tableStr = parse.search(searchStringTable, response.text)
#    print(tableStr)
    # splitting table into lines (rows)
    tableLines = tableStr[0].splitlines()
    # cutting header and footer of table
    header = 4
    footer = -1
    tableLinesShort = tableLines[header:footer]

    searchStringRow = "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}</a>"
    # could get rid of this footer if we had a dictionary or similar of
    # element info on possible charge states (equivalently number of electrons
    # in the neutral atom).

    elementObj = lev.element(element)
    chargeStatesNum = elementObj.atomic_number
#    print(f"number of charge states is {chargeStatesNum}")
    energiesArr = np.zeros((chargeStatesNum))
    uncertaintiesArr = np.zeros((chargeStatesNum))
    energyTypeArr = np.empty((chargeStatesNum), dtype=object)
#    print(energiesArr)
    for chargeState in np.arange(chargeStatesNum):
        #        print(f"####Charge state is {chargeState}####")
        #        print(tableLinesShort[chargeState])
        parsedRow = parse.search(searchStringRow, tableLinesShort[chargeState])
        (atNum,     # atomic number
         spName,     # species name
         ionCharge,     # ion charge
         elName,     # element name
         isoElSeq,    # isoelectronic sequence
         gndShells,     # ground shells
         shellRedun,     # redundant ending for ground shells
         gndLevel,     # ground level
         ionLevel,     # ionized level
         ionEnergyHTML,     # ionization energy (eV)
         uncertainty,  # uncertainty in the ionization energy (eV)
         refs    # references
         ) = [entry.strip() for entry in parsedRow]
    #    boom = [entry.strip() for entry in parsedRow]

        # further parsing to get at ionization energy and uncertainty
        ionEnergyStrip = cleanhtml(ionEnergyHTML)

        # test if value is measured, interpolated/extrapolated, or
        # if it is purely theoretical (nothing, square brackets, or parens)
    #    print(ionEnergyStrip)
        if ionEnergyStrip[0].isdigit():
            energyType = 'Measured'
            ionEnergyStr = ionEnergyStrip
        elif ionEnergyStrip[0] == '[':
            energyType = 'Interpolated/extrapolated'
            ionEnergyStr = ionEnergyStrip[1:-1]
        elif ionEnergyStrip[0] == '(':
            energyType = 'Theoretical'
            ionEnergyStr = ionEnergyStrip[1:-1]
        else:
            raise Exception('Something went wrong. String not parsed!')

    #    print(ionEnergyStr)
    #    print(uncertainty)
        # converting extracted strings to floats and adding units
        ionEnergy = float(ionEnergyStr)
        ionEnergyUnc = float(uncertainty)

        # record results into array
        energiesArr[chargeState] = ionEnergy
        uncertaintiesArr[chargeState] = ionEnergyUnc
        energyTypeArr[chargeState] = energyType

    # adding units
#    print(energiesArr)
    energiesArr = energiesArr * u.eV
    uncertaintiesArr = uncertaintiesArr * u.eV

    return energiesArr, uncertaintiesArr, energyTypeArr


def nistIonEnergyQuery(element, chargeState):
    """
    Given an element, and a charge state this
    script queries the NIST ionization energy level database and retrieves 
    the value of ionization energy and its uncertainty. The uncertainty
    is given in units of the last decimal place of the value.
    Elements (string)
    chargeState [dimless integer]
    Ionization energy [eV]
    Uncertainty [eV]
    """
    # checking inputs
    if not isinstance(element, str):
        raise ValueError("Element must be a string!")
    if not isinstance(chargeState, (int, np.int64)):
        raise ValueError("Charge state must be an integer!")
    # test when charge state is greater than atomic number
    elementObj = lev.element(element)
    atomicNum = elementObj.atomic_number
    if chargeState > (atomicNum - 1):
        raise Exception("Charge state cannot be greater than element's "
                        "atomic number!")
    # converting integer charge state to Roman numeral
#    chargeStateStr = intToRoman(chargeState)
    # formatting input strings
#    elementCharge = element + ' ' + chargeStateStr
    elementCharge = element
    # Starting requests session
    with req.session() as s:
        # values to be filled into website form
        payload = {'spectra': elementCharge,  # element and charge state
                   'submit': 'Retrieve Data',
                   'units': '1',  # selects energy units eV
                   'format': '1',  # ASCII output
                   'order': '0',  # order by Z or sequence
                   'at_num_out': 'on',  # output atomic number
                   'sp_name_out': 'on',  # output spectrum name
                   'ion_charge_out': 'on',  # output ion charge
                   'el_name_out': 'on',  # output element name
                   'seq_out': 'on',  # output isoelectronic sequence
                   'shells_out': 'on',  # output ground-state electronic shells
                   'conf_out': 'on',  # outpuse ground-state configuration
                   'level_out': 'on',  # output ground-state level
                   'ion_conf_out': 'on',  # output ionization configuration
                   'e_out': '0',  # selects ionization energy not binding
                   'unc_out': 'on',  # output uncertainty
                   'biblio': 'on'  # output bibliographic references
                   }
        # website's repsonse to query
        response = s.post(nistIonizationURL2, data=payload)
#        # for testing purposes only
#        print("URL:")
#        print(response.url)
#        print("STATUS:")
#        print(response.status_code)
#        print("TEXT:")
#        print(response.text)

    # Extract ionization energy value string from HTML/ASCII
    # Getting entire ASCII table
#    searchStringTable = "<pre>{}</pre>"
    searchStringTable = "<pre>{}<p"
    tableStr = parse.search(searchStringTable, response.text)
#    print(tableStr)
    # splitting table into lines (rows)
    tableLines = tableStr[0].splitlines()
    # cutting header and footer of table
    header = 4
    footer = -1
    tableLinesShort = tableLines[header:footer]

    # Here is where we can add a temporary storage for the queried
    # table so that the query doesn't have to be done multiple times
    # It may be worth splitting the function here into multiple functions.
    # One that grabs the table and another that parses it.

    # Each row has 10 columns of data
    searchStringRow = "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}</a>"
    parsedRow = parse.search(searchStringRow, tableLinesShort[chargeState - 0])
    (atNum,     # atomic number
     spName,     # species name
     ionCharge,     # ion charge
     elName,     # element name
     isoElSeq,    # isoelectronic sequence
     gndShells,     # ground shells
     shellRedun,     # redundant ending for ground shells
     gndLevel,     # ground level
     ionLevel,     # ionized level
     ionEnergyHTML,     # ionization energy (eV)
     uncertainty,  # uncertainty in the ionization energy (eV)
     refs    # references
     ) = [entry.strip() for entry in parsedRow]
#    boom = [entry.strip() for entry in parsedRow]

    # further parsing to get at ionization energy and uncertainty
    ionEnergyStrip = cleanhtml(ionEnergyHTML)

    # test if value is measured, interpolated/extrapolated, or
    # if it is purely theoretical (nothing, square brackets, or parens)
#    print(ionEnergyStrip)
    if ionEnergyStrip[0].isdigit():
        energyType = 'Measured'
        ionEnergyStr = ionEnergyStrip
    elif ionEnergyStrip[0] == '[':
        energyType = 'Interpolated/extrapolated'
        ionEnergyStr = ionEnergyStrip[1:-1]
    elif ionEnergyStrip[0] == '(':
        energyType = 'Theoretical'
        ionEnergyStr = ionEnergyStrip[1:-1]
    else:
        raise Exception('Something went wrong. String not parsed!')

#    print(ionEnergyStr)
#    print(uncertainty)
    # converting extracted strings to floats and adding units
    ionEnergy = float(ionEnergyStr) * u.eV
    ionEnergyUnc = float(uncertainty) * u.eV
#    # Return ionization energy and uncertainty
#    return ionEnergy, ionEnergyUnc
    return ionEnergy, ionEnergyUnc, energyType


#%% Calling function
# Set element, charge state, and temperature
#element = 'Na'
##element = 'He'
#chargeState = 4
#stuff = nistIonEnergyQuery(element, chargeState)
#print(f"Ionization energy: {stuff[0]}")
#
#energies, uncertainties, energyTypes = nistIonEnergiesQuery(element)
#print(f"Ionization energy: {energies}")
#
