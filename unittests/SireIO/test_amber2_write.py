
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

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
        assert_almost_equal( oldnrgs[key], newnrgs[key], 3 )

count = 0

def _test_write(files, verbose=False):

    global count
    count += 1

    try:
        s = MoleculeParser.read
    except:
        return

    if verbose:
        print("Reading... %s | %s" % (files[0],files[1]))

    s = MoleculeParser.read(files)

    # calculate the initial internal energy
    oldnrgs = _getEnergies(s)

    if verbose:
        print("Writing...")

    filenames = MoleculeParser.write(s, f"amber2_{count}_test")

    if verbose:
        print("Saved the system to file(s): %s" % filenames)

    # read this back in and check the energies
    s = MoleculeParser.read(filenames)    

    newnrgs = _getEnergies(s)

    if verbose:
        _printCompareEnergies(oldnrgs,newnrgs)

    _assert_almost_equal(oldnrgs, newnrgs)

def test_write(verbose=False):
    _test_write( ["../io/ose.top", "../io/ose.crd"], verbose )
    _test_write( ["../io/thrombin.top", "../io/thrombin.rst7"], verbose )
    _test_write( ["../io/proteinbox.top", "../io/proteinbox.crd"], verbose )

if __name__ == "__main__":
    test_write(True)

