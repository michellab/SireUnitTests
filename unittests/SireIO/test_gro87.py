try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *

from nose.tools import assert_equal, assert_almost_equal

# check that we have Gro87 support in this version of Sire
has_gro87 = True

try:
    p = Gro87()
except:
    # No Gro87 support
    has_gro87 = False


def _assert_strings_equal(x, y):
    assert_equal(len(x), len(y))

    for i in range(0, len(x)):
        assert_equal(x[i], y[i])


def _assert_equal(x, y, diff):
    if abs(x - y) > diff:
        assert_equal(x, y)


def _assert_almost_equal(v0, v1):
    assert_equal(len(v0), len(v1))

    is_vec = False

    try:
        v0[0].x()
        is_vec = True
    except:
        pass

    if is_vec:
        for i in range(0, len(v0)):
            _assert_equal(v0[i].x(), v1[i].x(), 0.01)
            _assert_equal(v0[i].y(), v1[i].y(), 0.01)
            _assert_equal(v0[i].z(), v1[i].z(), 0.01)
    else:
        for i in range(0, len(v0)):
            _assert_equal(v0[i], v1[i], 0.01)


def test_gro_rst(verbose=False):
    if not has_gro87:
        return

    # read in the system using a different format
    if verbose:
        print("Reading in rst/top file...")

    s = MoleculeParser.read("../io/NA16.rst", "../io/NA16.top")

    # Convert the data to Gro87
    if verbose:
        print("Converting...")

    g = Gro87(s)

    # Ensure that the numbers of atoms etc. are all correct
    if verbose:
        print("Converting to RST...")

    a = AmberRst(s)

    if verbose:
        print("Comparing...")

    assert_equal(g.nAtoms(), a.nAtoms())
    assert_equal(g.hasVelocities(), a.hasVelocities())

    _assert_almost_equal(g.coordinates(), a.coordinates())

    avels = a.velocities()
    # convert avels from angstrom / (20.455 ps) to nanometer / ps
    for i in range(0, len(avels)):
        avels[i] = (0.1 / 20.455) * avels[i]

    _assert_almost_equal(g.velocities(), avels)

    if verbose:
        print("Writing...")

    g.writeToFile("test.gro")

    if verbose:
        print("Re-reading...")

    g2 = Gro87("test.gro")

    if verbose:
        print("Comparing...")

    assert_equal(g.nAtoms(), g2.nAtoms())
    assert_equal(g.hasVelocities(), g2.hasVelocities())

    _assert_almost_equal(g.coordinates(), g2.coordinates())
    _assert_almost_equal(g.velocities(), g2.velocities())
    _assert_almost_equal(g.atomNumbers(), g2.atomNumbers())
    _assert_strings_equal(g.atomNames(), g2.atomNames())
    _assert_almost_equal(g.residueNumbers(), g2.residueNumbers())
    _assert_strings_equal(g.residueNames(), g2.residueNames())

    assert_equal(g.boxV1(), g2.boxV1())
    assert_equal(g.boxV2(), g2.boxV2())
    assert_equal(g.boxV3(), g2.boxV3())


def test_gro87(verbose=False):
    if not has_gro87:
        return

    if verbose:
        print("Reading the Gro87 file using the Gro87 parser...")

    grofile = "../io/water.gro"
    # grofile = "../io/urea.gro"

    gro87 = Gro87(grofile)

    if verbose:
        print(gro87)


if __name__ == "__main__":
    test_gro87(True)
    test_gro_rst(True)
