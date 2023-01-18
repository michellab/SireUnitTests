try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.MM import *
from Sire.FF import *

import os

from nose.tools import assert_equal, assert_almost_equal


def _assert_vec_equal(vec0, vec1):
    assert_almost_equal(vec0.x(), vec1.x(), 5)
    assert_almost_equal(vec0.y(), vec1.y(), 5)
    assert_almost_equal(vec0.z(), vec1.z(), 5)


def _assert_vecs_equal(vals0, vals1):
    assert_equal(len(vals0), len(vals1))

    for i in range(0, len(vals0)):
        _assert_vec_equal(vals0[i], vals1[i])


def _getEnergies(s):
    intraclj = IntraFF("intraclj")
    intraclj.add(s.molecules())

    intraff = InternalFF("intraff")
    intraff.setUse14Calculation(True)
    intraff.add(s.molecules())

    ffs = ForceFields()
    ffs.add(intraclj)
    ffs.add(intraff)

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

    for key in oldkeys:
        assert_almost_equal(oldnrgs[key], newnrgs[key], 5)


def test_amberrst(verbose=False):

    if verbose:
        print("Reading...")

    s = MoleculeParser.read("../io/NA16.top", "../io/NA16.rst")

    if verbose:
        print("Converting...")

    newfile = AmberRst(s)

    if verbose:
        print(newfile)
        print(newfile.creatorApplication())

    if verbose:
        print("Loading original...")

    oldfile = AmberRst("../io/NA16.rst")

    if verbose:
        print("Comparing...")

    assert_equal(newfile.nFrames(), oldfile.nFrames())
    _assert_vecs_equal(newfile.coordinates(), oldfile.coordinates())
    _assert_vecs_equal(newfile.velocities(), oldfile.velocities())
    _assert_vecs_equal(newfile.forces(), oldfile.forces())
    _assert_vec_equal(newfile.boxDimensions(), oldfile.boxDimensions())
    _assert_vec_equal(newfile.boxAngles(), oldfile.boxAngles())
    assert_equal(newfile.createdFromRestart(), oldfile.createdFromRestart())
    assert_almost_equal(newfile.time(), oldfile.time())

    # save the file and reload
    if verbose:
        print("Writing to a temporary file...")

    newfile.writeToFile("file_amberrst_test.rst")

    if verbose:
        print("Reloading from the temporary file...")

    new2file = AmberRst("file_amberrst_test.rst")

    print(new2file.warnings())

    if verbose:
        print("Comparing the data...")

    assert_equal(new2file.nFrames(), oldfile.nFrames())
    _assert_vecs_equal(new2file.coordinates(), oldfile.coordinates())
    _assert_vecs_equal(new2file.velocities(), oldfile.velocities())
    _assert_vecs_equal(new2file.forces(), oldfile.forces())
    _assert_vec_equal(new2file.boxDimensions(), oldfile.boxDimensions())
    _assert_vec_equal(new2file.boxAngles(), oldfile.boxAngles())
    assert_equal(new2file.createdFromRestart(), oldfile.createdFromRestart())
    assert_almost_equal(new2file.time(), oldfile.time())

    if verbose:
        print("Writing back to a full file...")

    for t in ("file_amberrst_test.rst", "file_amberrst_test.prm7"):
        try:
            os.remove(t)
        except:
            pass

    if verbose:
        print("Writing...")

    MoleculeParser.write(s, "file2_amberrst_test")

    if verbose:
        print("Re-reading the data from the written file...")

    s2 = MoleculeParser.read("file2_amberrst_test.prm7", "file2_amberrst_test.rst")

    if verbose:
        print("Comparing energies...")

    oldnrgs = _getEnergies(s)
    newnrgs = _getEnergies(s2)

    if verbose:
        _printCompareEnergies(oldnrgs, newnrgs)

    _assert_almost_equal(oldnrgs, newnrgs)


if __name__ == "__main__":
    test_amberrst(True)
