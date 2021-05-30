#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 20:13:55 2018

Query tools for Henke filters site.

@author: Pawel M. Kozlowski
"""

# python modules
import numpy as np
import requests as req
import parse as parse
import re
import astropy.units as u
import matplotlib.pyplot as plt



#%% simple utilities

# URL to nist ionization energy database
henkeURL1 = 'http://henke.lbl.gov/optical_constants/filter2.html'
henkeURL2 = 'http://henke.lbl.gov/cgi-bin/filter.pl'
henkeURL3 = 'http://henke.lbl.gov/'

# Checking website status
#reqObj = req.post(henkeURL2)
#print("STATUS:")
#print(reqObj.status_code)
#if reqObj.status_code == req.codes.ok:
#    print("OK!")
#else:
#    print("CONNECTION FAILED!")
#reqObj.raise_for_status()
    
# for testing purposes only
#print("HEADERS:")
#print(reqObj.headers['content-type'])
#print("ENCODING:")
#print(reqObj.encoding)
#print("TEXT:")
#print(reqObj.text)


#%% query tool(s)
# list of default materials in Henke filters database
# material name : (chemical formula, density in gm/cm^3)
defaultMaterials = {'polyimide' : ('C22H10N2O5', 1.43),
                    'boron nitride' : ('BN', 2.25),
                    'silicon nitride' : ('Si3N4', 3.44),
                    'polypropylene' : ('C3H6', 0.90),
                    'PMMA' : ('C5H8O2', 1.19),
                    'polycarbonate' : ('C16H14O3', 1.2),
                    'mylar' : ('C10H8O4', 1.4),
                    'Teflon' : ('C2F4', 2.2),
                    'Parylene-C' : ('C8H7Cl', 1.29),
                    'Parylene-N' : ('C8H8', 1.11)}



def fetchLink(text):
    """
    Parses text resulting from Henke filter query to retrieve temporary
    hyperlink to data table.
    """
    searchString = "<a HREF={}>"
    tableStr = parse.search(searchString, text)
    urlEnding = tableStr[0][2:-1]
    fullURL = henkeURL3 + urlEnding
    return fullURL


def henkeParse(text):
    """
    Parses transmission data table from Henke filters page and returns
    an array of values.
    """
    # split along new lines
    splitText = text.splitlines()
    
    # parsing first header line, which is of the form
    # ' Cu Density=8.96 Thickness=0.2 microns'
    header1 = splitText[0].strip().split(' ')
    chemical = header1[0]
    density = float(header1[1].split('=')[1])
    thickness = float(header1[2].split('=')[1])
    
    
    # parsing second header line, which is of the form
    # ' Photon Energy (eV), Transmission'
    header2 = splitText[1].strip().split(',')
    spectralUnit = header2[0]
    
    
    # getting array of transmission data strings
    transmissionStr = splitText[2:]
    # stripping whitespace on either side, splitting spectral and
    # transmission values
    transmissionStr2 = [transText.strip().split() for transText in transmissionStr]
#    print(transmissionStr2)
    # wrap into numpy array
    transmissionArr = np.array(transmissionStr2)
    # splitting spectral and transmission values into separate arrays
    spectralStr = transmissionArr[:,0]
    transmissionStr3 = transmissionArr[:,1]
    # forcing transmission data to floats
    spectralArr = [float(spect) for spect in spectralStr]
    transmissionArr = [float(trans) for trans in transmissionStr3]
    # adding units to transmission data
    if spectralUnit == "Photon Energy (eV)":
        spectralArr = spectralArr * u.eV
    elif spectralUnit == "Wavelength (nm)":
        spectralArr = spectralArr * u.nm
    else:
        raise Exception("Something changed with Henke site. Contact "
                        "the developer to raise this bug as an issue.")
    transmissionArr = transmissionArr * u.dimensionless_unscaled
    
    return (spectralArr,
            transmissionArr,
            chemical,
            density,
            thickness,
            spectralUnit)


@u.quantity_input(density=u.g * u.cm ** -3,
                  thickness=u.um)
def queryHenke(formula="Si3N4",
               density=-1 * u.g * u.cm ** -3,
               thickness=0.2 * u.um,
               spectralUnit="Energy",
               minRange=10,
               maxRange=1000,
               npts=100,
               plotScaling="Linear"):
    """
    Queries Henke filters page.
    density defaults to -1, which uses tabulated values.
    thickness defaults to 0.2 um.
    
    """
    # checking inputs
    if not isinstance(formula, str):
        raise ValueError("formula must be a string!")
    # check if plot scaling is valid
    if not plotScaling in ["Linear", "Log", "LogLin", "LinLog"]:
        raise Exception(f"Invalid plot scaling {plotScaling} given.")
    if minRange >= maxRange:
        raise Exception("minRange must be less than maxRange!")
    # ensuring spectral range is within limits of Henke tool
    if spectralUnit == "Energy":
        if minRange < 10 or maxRange > 30000:
            raise Exception("Photon energies must be between 10 and 30,000 eV")
    elif spectralUnit == "Wave":
        if minRange < 0.041 or maxRange > 124:
            raise Exception("Wavelengths must be between 0.041 and 124 nm")
    else:
        raise Exception(f"SpectralUnit must be either Energy or Wave, and "
                        f"not {spectralUnit}")
        
    
    # Starting requests session
    with req.session() as s:
        # values to be filled into website form
        payload = {'Material': 'Enter Formula', # element and charge state
                   'Formula': formula,
                   'Density' : (density.to(u.g * u.cm ** -3)).value,
                   'Thickness' : (thickness.to(u.um)).value,
                   'Scan' : spectralUnit,
                   'Min' : minRange,
                   'Max' : maxRange,
                   'Npts' : npts,
                   'Plot' : plotScaling,
                   'Output' : 'Plot'
                   }
        # website's repsonse to query
        response1 = s.post(henkeURL2, data=payload)    
        # extracting hyperlink to the transmission data table    
        dataLink = fetchLink(response1.text)
        # querying transmission data table
        response2 = s.get(dataLink)
    # parsing transmission data table
    data = henkeParse(response2.text)
    return data


def henkeFilterPlot(data):
    """
    Takes data from queryHenke() and plots it.
    """
    plt.plot(data[0], data[1])
    plt.xlabel(data[5])
    plt.ylabel('Transmission')
    plt.title(data[2])
    plt.show()
    return

#%% testing query
#data = queryHenke(formula="C",
#                  thickness=1*u.um,
#                  minRange=0.05,
#                  maxRange=100,
#                  npts=100,
#                  spectralUnit="Wave")
data = queryHenke(formula="C",
                  thickness=1*u.um,
                  minRange=10,
                  maxRange=1000,
                  npts=100,
                  spectralUnit="Energy")

henkeFilterPlot(data)

