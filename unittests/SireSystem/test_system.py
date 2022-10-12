from sire.legacy.System import *
from sire.legacy.FF import *
from sire.legacy.MM import *
from sire.legacy.IO import *
from sire.legacy.Mol import *
from sire.legacy.Maths import *
from sire.legacy.CAS import *
from sire.legacy.Vol import *
from sire.legacy.Units import *
from sire.legacy.Qt import *

from nose.tools import assert_almost_equal

def test_system(verbose=False):

    cljff = InterCLJFF()

    mincoords = Vector(-18.3854, -18.66855, -18.4445)
    maxcoords = Vector( 18.3854,  18.66855,  18.4445)

    vol = PeriodicBox(mincoords, maxcoords)

    cljff.setSpace(vol)

    if verbose:
        print("Reading pdb...")

    mols = PDB().read("../io/water.pdb")
                                                
    if verbose:
        print(mols)
        print("Read in %d molecules!" % mols.nMolecules())

    i = 0

    if verbose:
        print("Getting the first molecule...")

    mol = mols.moleculeAt(0).molecule()

    if verbose:
        print("About to edit...")

    mol = mol.edit().atom( AtomName("O00") ) \
                        .setProperty("LJ", LJParameter(3.15363*angstrom,  \
                                                   0.1550*kcal_per_mol)).molecule() \
                    .atom( AtomName("H01") ) \
                        .setProperty("charge", 0.520 * mod_electron).molecule() \
                    .atom( AtomName("H02") ) \
                        .setProperty("charge", 0.520 * mod_electron).molecule() \
                    .atom( AtomName("M03") ) \
                        .setProperty("charge", -1.04 * mod_electron).molecule() \
             .commit()

    if verbose:
        print("Edited molecule...")

    charges = mol.property("charge")
    ljs = mol.property("LJ")

    cljff.add(mol)

    if verbose:
        print("Setting charges...")

    for i in range(1, mols.nMolecules()):
        mol = mols.moleculeAt(i).molecule()

        mol = mol.edit().setProperty("charge", charges) \
                        .setProperty("LJ", ljs) \
                 .commit()

        cljff.add(mol)

    if verbose:
        print("Creating system...")

    system = System()

    system.add(cljff)

    nrg = system.energy()

    if verbose:
        print("System energy = %s" % (nrg))

    copy_system = System(system)

    nrg2 = copy_system.energy()

    if verbose:
        print("Copy energy: %s (should be %s)" % (nrg2,nrg))

    assert_almost_equal(nrg.value(), nrg2.value(), 5)

    copy_system.mustNowRecalculateFromScratch()
    nrg3 = copy_system.energy()

    if verbose:
        print("Copy energy: %s (should be %s)" % (nrg2,nrg))

    assert_almost_equal(nrg.value(), nrg3.value(), 5)

if __name__ == "__main__":
    test_system(True)
