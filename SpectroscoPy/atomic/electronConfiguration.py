#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 02:49:32 2017

Class for describing electron configuration for an ion. This will be
used with the term symbol to give a full description of the
energy state of an ion. Two such energy states then define a transition.
There will be a converter function which converts between this class
object description and spectroscopic notation as found in databases
such as NIST.


TODO:
    -add square bracket functionality NIST
    - __bool__ behavior for configurations, terms, energies, etc
    -include L-S coupling
    -include functions/methods without L-S coupling assumption
    -use subshellConfiguration() in electronConfiguration()
    -use electronConfiguration() and termSymbol() for energyLevel()
    -use 2 energyLevel() objects for transition() object
    -Create system for managing/saving/converting list of transition() objs
    -add tests for allowed/forbidden transitions

@author: Pawel M. Kozlowski
"""

# importing python packages
import numpy as np
import pandas as pd
# importing custom packages

# defining superscript numbers
superNums = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
# defining subscript numbers
subNums = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

# ion (element and ionization stage)
# this defines the total number of electrons


# maximum electrons in a shell
def shellElectrons(n):
    """
    Given principal quantum number, n, return maximum number of electrons
    contained in shell n.
    """
    return 2 * n**2


# magnetic quantum number (defines specific orbital within subshell)
"""
ml range from -l to +l by integer steps
"""

# spin projection quantum number (describes spin of electron in given orbital)
"""
ms range from -s to +s by integer steps
for an electron this mean -1/2 or +1/2, therefore only 2 electrons per
orbital due to Pauli exclusion principle
"""


# multiplicity (singlet, doublet, triplet, etc)
"""
multiplicity = 2(S+1)
where S = total spin angular momentum, determined by the total spin of
unpaired electrons. Electrons which are already paired in an orbital are
of opposite spin and therefore cancel out.
multiplicity = 1  singlet
multiplicity = 2  doublet
multiplicity = 3  triplet
multiplicity = 4  quartet
multiplicity = 5  quintet

The number of levels is further determined by the azimuthal quantum number
as well. When S <= L the multiplicity equals the number of spin orientations,
but when S > L there are only 2L+1 orientations of total angular momentum
possible (where total angular momentum is J = L + S). 

example:
    S = 3/2 (3 unpaired electrons)
    2S +1 = 4 (quartet)
    For an S state, L = 0
    S can range from -S to + S in integer steps, but L is fixed, and J must
    be positive. So only possible value of J is J = 0 + 3/2, and there is
    only one level even though multiplicity is 4.
    
Note that electrons will first fill empty orbitals, so it should be easy to
calculate the number of unpaired electrons, and therefore the  multiplicity.


One way to get multiplicity is to take total # of possible electrons in 
a subshell and divide it by 2 (which then equals 2l+1).
With nElec being the total number of electrons in a subshell...
if nElec <= 2l+1:
    # no electrons are paired
    S = nElec * 1/2
elif nElec > 2l+1:
    # some electrons will be paired, but how many?
    paired = nElec - 2l+1
    # S is then nElec, less the paired electrons, times 1/2
    S = (nElec - 2*paired) * 1/2
    # reducing, which is just the maximum electrons in the subshell, less
    # the total number of present electrons. This works because the empty
    # spots tell us there are unpaird electrons, since over half the
    # subshell is filled!
    S = (2(2l+1) - nElec) * 1/2
"""


# fine structure splitting is due to spin-orbit coupling of electron

# hyperfine structure splitting (would need to include magnetic dipole
# moment of the nucleus)


#%%
# term, level, state, parity

# x-ray spectroscopic notation for principal quantum numbers (shells)
xrayNotation = ['K',
                'L',
                'M',
                'N',
                'O',
                'P',
                'Q',
                'R',
                'S',
                'T',
                'U',
                'V',
                'W',
                'X',
                'Y',
                'Z',
                ]
principals = np.arange(0, len(xrayNotation)) + 1
# constructing dataframe for converting x-ray notation symbols to
# principal quantum numbers
xrayShellsDf = pd.DataFrame(data=principals, index=xrayNotation, columns=['n'])
# example looking up number using symbol, in this case looking up
# shell 'M'
print(xrayShellsDf['n']['M'])


symbols = ['s', # sharp
           'p', # principal
           'd', # diffuse
           'f', # fundamental/fine
           'g',
           'h',
           'i',
           'k',
           'l',
           'm',
           'n',
           'o',
           'q',
           'r',
           't',
           'u',
           'v',
           'w',
           'x',
           'y',
           'z',
           ]
azis = np.arange(0, len(symbols))

# constructing dataframe for converting spdf symbols to azimuthal quantum
# number and back again
subshellsDf = pd.DataFrame(data=azis, index=symbols, columns=['l'])
subshellsDf['sym'] = symbols

# example looking up number using symbol, in this case looking up
# subshell 'f'
print(subshellsDf['l']['f'])
# example looking up symbol using number, in this case looking up
# symbol for subshell l=3
print(subshellsDf['sym'][3])

#%% subshell configurations
class SubshellConfiguration():
    """
    Object for describing atomic subshells containing electrons.
    """

    def __init__(self, principalNum, aziNum, electronsNum):
        """
        Populating a specific subshell, described by principal quantum
        number and azimuthal quantum number, with electronsNum number of
        electrons.
        This will check whether the number of electrons is allowed.
        """
        # check that given numbers are within allowable range (positive)
        if not principalNum >= 1:
            raise ValueError("Principal quantum number must be n >= 1.")
        if not aziNum >= 0:
            raise ValueError("Azimuthal quantum number must be l >= 0.")
        # assigning principal and azimuthal quantum numbers to object
        self.principal = principalNum
        self.aziNum = aziNum
        # check if azimuthal quantum number exists for principal quantum number
        self.existsSubshell()
        # check if number of electrons is allowed
        maxElectrons = self.subshellElectrons()
        if electronsNum > maxElectrons:
            raise ValueError(f"Given number of electrons is {electronsNum}, "
                             f"which is more than maximum, {maxElectrons}.")
        # assign electrons
        self.electrons = electronsNum
    def n(self):
        """
        return principal quantum number.
        """
        return self.principal

    def l(self):
        """
        return azimuthal quantum number.
        """
        return self.aziNum

    def elec(self):
        """
        return number of electrons in subshell.
        """
        return self.electrons
    def __str__(self):
        """
        Displays subshell configuration in standard spectroscopic notation.
        """
        # convert principal and azimuthal quantum numbers into symbols
        principalSym = str(self.principal)
        aziSym = subshellsDf['sym'][self.aziNum]
        # convert number of electrons to superscript
        electronsSym = (str(self.electrons)).translate(superNums)
        # concatenating
        subshellSym = principalSym + aziSym + electronsSym
        return subshellSym
    def __eq__(self, other):
        """
        Test if 2 subshell configurations are identical.
        """
        # check if other is a SubshellConfigration() object
        if not isinstance(other, SubshellConfiguration):
            raise ValueError(f"{other} is not a SubshellConfiguration().")
        principalEq = self.n() == other.n()
        aziEq = self.l() == other.l()
        electronsEq = self.elec() == other.elec()
        totalEq = principalEq and aziEq and electronsEq
        return totalEq
    def existsSubshell(self, principal=np.nan, aziNum=np.nan):
        """
        Checks if azimuthal quantum number (subshell) exists within the
        principal quantum number.
        """
        # if arguments aren't passed then use parameters given in object
        if np.isnan(principal):
            principal = self.principal
        if np.isnan(aziNum):
            aziNum = self.aziNum
        # azimuthal number must be 0 <= l <= n-1
        if not aziNum <= (principal - 1):
            raise ValueError(f"Azimuthal number {aziNum} does not exist "
                             f"for principal quantum number {principal}.")
        return

    def subshellElectrons(self, aziNum=None):
        """
        Given azimuthal quantum number, return maximum number of electrons
        in the subshell.
        Symbol for azimuthal quantum number is lowercase, cursive L.
        """
        # if parameters aren't passed then use parameters given in object
        if aziNum == None:
            aziNum = self.aziNum
        return 2 * (2 * aziNum + 1)
    def str2obj(configStr):
        """
        Method for generating a subshell configuration from spectroscopic
        notation string.
        """
        # split string into components
        
        # test that string is valid
        
        # generate subshell configuration object
        subshellConfig = 
        return subshellConfig


# testing subshell class
testSubshell = SubshellConfiguration(principalNum=2,
                                     aziNum=1,
                                     electronsNum=5)
print(testSubshell)


#%% full electron configurations

# One way to describe the electron configuration is as a collection of
# subshell configurations. Undefined subshells are considered empty.
# Could add function to subshell configuration to define and populate
# the object by using spectroscopic notation.

class ElectronConfiguration():
    """
    Object for describing the electron configuration of an atom/ion.
    """

    def __init__(self, subshells):
        """
        Initializes electron configuration using a list of 
        SubshellConfiguration() objects. It is assumed that
        any subshells which are not given are empty.
        """
        # test that subshells is a list or numpy array
        
        # convert subshells into a numpy array if it is a list
        
        
        # test that all objects in list are SubshellConfiguration()
        # objects
        
        # test that there are no conflicting subshells, i.e., same
        # principal and azimuthal quantum numbers
        
        # drop any empty subshells from the list
        
        # organize the list of subshells according to principal
        # quantum number and azimuthal number
        subshellsOrdered = 
        
        # if an empty list is passed then an empty configuration object 
        # is created. This represents a fully ionized particle
        self.subshells = subshellsOrdered
    def __eq__(self, other):
        """
        Check if 2 electron configurations are equal. This is done
        by checking if lists of subshells are equal.
        """
        # test if argument is an ElectronConfiguration()
        if not isinstance(other, ElectronConfiguration):
            raise ValueError(f"{other} is not an ElectronConfiguration().")
        # test equivalence of subshells
        eq = (other.subshells == self.subshells).all()
        return eq
    def __bool__(self):
        """
        Method for determining if class value is True or False.
        ElectronConfigration() is False when it is empty, and True when
        there is at least one subshell.
        """
        # since subshells is a list which describes the configuration, when
        # the list is empty then the configuration is False.
        return bool(self.subshells)
    def __str__(self):
        """
        An electron configuration consists of an ordered list of subshell
        configurations.
        !!!May also include quantities in square brackets, but still
        need to figure out what those mean!!!
        """
    def str2obj(self):
        """
        Initialize an electron configuration using spectroscopic notation
        string.
        """
        # break string into substrings, with each one representing
        # a subshell configuration.
        
        # test that string is valid
        
        # combine subshell configurations to form full electron configuration.
        
        return configObj
    def subshell(self, principalNum, aziNum):
        """
        Return SubhshellConfiguration() object for a particular subshell
        corresponding to principal and azimuthal quantum numbers.
        """
        return
    def empty():
        """
        Generate an empty configuration representing a fully ionized
        particle.
        """
        return ElectronConfiguration([])
    


#%% term symbol description

class TermSymbol():
    """
    """

    def __init__(self):
        """
        """

    def __str__(self):
        """
        """
        # spin state (e.g. doublet, triplet). This becomes the superscript

        # azimuthal quantum number, L. This becomes the letter symbol

        # total angular momentum quantum number, J. This becomes the
        # subscript.


#%% full energy level description, this includes an electron configuration
# and a term symbol

class EnergyLevel():
    """
    """
    def __init__(self, electronConfiguration, termSymbol):
        """
        """
        # initialize electron configuration
        self.electronConfiguration = electronConfiguration
        # initialize term symbol
        self.termSymbol = termSymbol

    def __str__(self):
        """
        String representation of enery level. This is just the string
        representations of the electron configuration and term symbol.
        """
        # fetching electron configuration and term symbol strings.
        configStr = self.electronConfiguration.__str__()
        termStr = self.termSymbol.__str__()
        # concatenating
        energyStr = configStr + ' ' + termStr
        return energyStr
    def configuration(self):
        """
        Returns full electron configuration.
        """
        return self.electronConfiguration
    def term(self):
        """
        Returns term symbol.
        """
        return self.termSymbol
