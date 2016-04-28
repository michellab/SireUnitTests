
from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *

import os,re,sys
import shutil

combining_rules = "arithmetic"
temperature = 25 * celsius
pressure = 1 * atm
coulomb_cutoff = 1000 * angstrom
coulomb_feather = 999.5 * angstrom
lj_cutoff = 1000 * angstrom
lj_feather = 999.5 * angstrom
#############################################################

def test_nrg(verbose=False):

    top_file = "../io/ose.top"
    crd_file = "../io/ose.crd"

    amber = Amber()
    molecules, space = amber.readCrdTop(crd_file, top_file)
    # Overload, we want to calc the energy in a non periodic box for comparison with Sander
    space = Cartesian()

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    system = System()

    solute = MoleculeGroup("solute", moleculeList[0])

    # Add these groups to the System
    system.add(solute)

    # Now solute bond, angle, dihedral energy
    solute_intraff = InternalFF("solute_intraff")
    solute_intraff.add(solute)

    # Now solute intramolecular CLJ energy
    solute_intraclj = IntraCLJFF("solute_intraclj")
    solute_intraclj.add(solute)

    # Here is the list of all forcefields
    forcefields = [ solute_intraff, solute_intraclj ]
 
    # Add these forcefields to the system
    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty( "space", space )
    system.setProperty( "switchingFunction", 
                        HarmonicSwitchingFunction(coulomb_cutoff, coulomb_feather,
                                                  lj_cutoff, lj_feather) ) 
    system.setProperty( "combiningRules", VariantProperty(combining_rules) )

    total_nrg = solute_intraclj.components().total() + solute_intraff.components().total()

    e_total = system.totalComponent()
    system.setComponent( e_total, total_nrg )

    if verbose:
        print("\nTotal energy ")
        print(system.energy())

        print("Components energies ")
        for component in list(system.energyComponents().keys()):
            print(component, system.energyComponents().value(component) * kcal_per_mol)

        print("The AMBER11/sander energies for this system are ") 
        print("""
""")

if __name__ == "__main__":
    test_nrg(True)

