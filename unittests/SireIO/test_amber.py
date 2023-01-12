try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

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

import os, re, sys
import shutil

from nose.tools import assert_almost_equal

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
    forcefields = [solute_intraff, solute_intraclj]

    # Add these forcefields to the system
    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty(
        "switchingFunction",
        HarmonicSwitchingFunction(
            coulomb_cutoff, coulomb_feather, lj_cutoff, lj_feather
        ),
    )
    system.setProperty("combiningRules", wrap(combining_rules))

    total_nrg = (
        solute_intraclj.components().total()
        + solute_intraff.components().total()
    )

    e_total = system.totalComponent()
    system.setComponent(e_total, total_nrg)

    if verbose:
        print("\nTotal energy ")
        print(system.energy())

        print("Components energies ")
        for component in list(system.energyComponents().keys()):
            print(
                component,
                system.energyComponents().value(component) * kcal_per_mol,
            )

        print("The AMBER14/sander energies for this system are ")
        print(
            """
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1      -3.4880E+01     1.3388E+01     6.2845E+01     C2         20

 BOND    =        5.1844  ANGLE   =       10.0783  DIHED      =       20.1271
 VDWAALS =       -1.7278  EEL     =     -256.9757  HBOND      =        0.0000
 1-4 VDW =       10.6377  1-4 EEL =      177.7958  RESTRAINT  =        0.0000
"""
        )

        diff = abs(-34.880 - system.energy().value())
        print("Difference = %s" % diff)

    e_bond = system.energy(solute_intraff.components().bond()).value()
    e_ang = system.energy(solute_intraff.components().angle()).value()
    e_dih = (
        system.energy(solute_intraff.components().dihedral()).value()
        + system.energy(solute_intraff.components().improper()).value()
    )

    assert_almost_equal(e_bond, 5.1844, 2)
    assert_almost_equal(e_ang, 10.0783, 2)
    assert_almost_equal(e_dih, 20.1271, 2)

    e_coul = system.energy(solute_intraclj.components().coulomb()).value()
    e_lj = system.energy(solute_intraclj.components().lj()).value()

    assert_almost_equal(e_coul, -256.9757 + 177.7958, 2)
    assert_almost_equal(e_lj, -1.7278 + 10.6377, 2)

    assert_almost_equal(system.energy().value(), -34.880, 2)


if __name__ == "__main__":
    test_nrg(True)
