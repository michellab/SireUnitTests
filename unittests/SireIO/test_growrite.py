
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.FF import *

from nose.tools import assert_equal, assert_almost_equal

def _getEnergies(s):
    intraclj = IntraFF("intraclj")
    intraclj.add(s.molecules())

    intraff = InternalFF("intraff")
    intraff.setUse14Calculation(True)
    intraff.add(s.molecules())

    interff = InterFF("interff")
    interff.add(s.molecules())

    ffs = ForceFields()
    ffs.add(intraclj)
    ffs.add(intraff)
    ffs.add(interff)

    return ffs.energies()

def _printCompareEnergies(oldnrgs, newnrgs):
    keys = list(oldnrgs.keys())
    keys.sort()

    for key in keys:
        print("%s: %s  %s" % (key, oldnrgs[key],newnrgs[key]))

def _assert_almost_equal(oldnrgs, newnrgs):
    oldkeys = list(oldnrgs.keys())
    newkeys = list(newnrgs.keys())

    oldkeys.sort()
    newkeys.sort()

    assert_equal( oldkeys, newkeys )

    for key in oldkeys:
        assert_almost_equal( oldnrgs[key], newnrgs[key], 5 )

def test_growrite(verbose=False):

    try:
        s = MoleculeParser.read
    except:
        return

    if verbose:
        print("Reading...")

    s = MoleculeParser.read("../io/urea.top", "../io/urea.gro", 
                            {"GROMACS_PATH":"../io/gromacs"})

    # calculate the initial internal energy
    oldnrgs = _getEnergies(s)

    if verbose:
        print("Writing...")

    filenames = MoleculeParser.write(s, "test")

    if verbose:
        print("Saved the system to file(s): %s" % filenames)

    # read this back in and check the energies
    s = MoleculeParser.read(filenames)    

    newnrgs = _getEnergies(s)

    if verbose:
        _printCompareEnergies(oldnrgs,newnrgs)

    _assert_almost_equal(oldnrgs, newnrgs)

if __name__ == "__main__":
    test_growrite(True)

