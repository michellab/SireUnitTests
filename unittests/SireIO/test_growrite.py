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
        print("%s: %s  %s" % (key, oldnrgs[key], newnrgs[key]))


def _assert_almost_equal(oldnrgs, newnrgs):
    oldkeys = list(oldnrgs.keys())
    newkeys = list(newnrgs.keys())

    oldkeys.sort()
    newkeys.sort()

    assert_equal(oldkeys, newkeys)

    newsum = 0
    oldsum = 0

    for key in oldkeys:
        if str(key).find("dihedral") == -1 and str(key).find("improper") == -1:
            assert_almost_equal(oldnrgs[key], newnrgs[key], 3)
        else:
            newsum += newnrgs[key]
            oldsum += oldnrgs[key]

    # this is the sum of improper and dihedral energy
    assert_almost_equal(oldsum, newsum, 3)


def _test_ambergro(files, verbose=False):

    if verbose:
        print("Reading Amber file...")

    s = MoleculeParser.read(files)

    if verbose:
        print("Calculating energies...")

    oldnrgs = _getEnergies(s)

    if verbose:
        print("Converting to gromacs format... topology...")

    g = GroTop(s)
    g.writeToFile("test.top")

    if verbose:
        print("...coordinates/velocities...")

    g87 = Gro87(s)
    g87.writeToFile("test.gro")

    if verbose:
        print("Re-reading from the gromacs file...")

    s = MoleculeParser.read("test.top", "test.gro")

    if verbose:
        print("Calculating energies...")

    newnrgs = _getEnergies(s)

    if verbose:
        _printCompareEnergies(oldnrgs, newnrgs)

    _assert_almost_equal(oldnrgs, newnrgs)

    if verbose:
        print("Writing back to amber...")

    # write back to amber and test
    a = AmberPrm(s)
    a.writeToFile("test.prm")

    if verbose:
        print("Reading back...")

    s = MoleculeParser.read("test.prm", "test.gro")

    if verbose:
        print("Getting energies...")

    newnrgs = _getEnergies(s)

    if verbose:
        _printCompareEnergies(oldnrgs, newnrgs)

    _assert_almost_equal(oldnrgs, newnrgs)


def test_ambergro(verbose=False):
    if verbose:
        print("\nTesting ose.top/ose.crd")

    _test_ambergro(["../io/ose.top", "../io/ose.crd"], verbose)

    # if verbose:
    #    print("\nTesting ethanol.grotop/ethanol.gro")

    # _test_ambergro(["../io/ethanol.grotop","../io/ethanol.gro"], verbose)


def test_growrite(verbose=False):

    if verbose:
        print("Reading...")

    s = MoleculeParser.read(
        "../io/urea.top", "../io/urea.gro", {"GROMACS_PATH": "../io/gromacs"}
    )

    # calculate the initial internal energy
    oldnrgs = _getEnergies(s)

    if verbose:
        print("Writing...")

    filenames = MoleculeParser.write(s, "test")

    if verbose:
        print("Saved the system to file(s): %s" % filenames)

    # read this back in and check the energies
    s = MoleculeParser.read(filenames)

    newnrgs = _getEnergies(s)

    if verbose:
        _printCompareEnergies(oldnrgs, newnrgs)

    _assert_almost_equal(oldnrgs, newnrgs)


if __name__ == "__main__":
    test_ambergro(True)
    test_growrite(True)
