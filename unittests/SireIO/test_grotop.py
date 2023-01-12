try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *
from Sire.Base import *

from nose.tools import assert_equal, assert_almost_equal
from nose.tools import assert_true, assert_false

# check that we have GroTop support in this version of Sire
has_grotop = True

# Set the GROMACS_PATH to "../io/gromacs" as this is the
# directory that contains all of the gromacs include files
gromacs_path = StringProperty("../io/gromacs")

try:
    p = GroTop()
except:
    # No GroTop support
    has_grotop = False


def test_grotop(verbose=False):
    if not has_grotop:
        return

    # read in the system using a different format
    if verbose:
        print("Reading in gromacs top file...")

    files = ["../io/urea.top", "../io/urea.gro"]

    s = MoleculeParser.read(files, {"GROMACS_PATH": gromacs_path})

    # check parameters...
    for i in range(0, s.nMolecules()):
        m = s[MolIdx(i)]

        print(m.propertyKeys())

        for atom in m.atoms():
            print(atom, atom.property("charge"), atom.property("LJ"))
            print(atom.property("atomtype"), atom.property("charge_group"))

    p = PDB2(s)

    p.writeToFile("test.pdb")

    m = Mol2(s)

    m.writeToFile("test.mol2")

    g = Gro87(s)

    g.writeToFile("test.gro")

    prm = AmberPrm(s)

    prm.writeToFile("test.prm")


def test_grosys(verbose=False):

    try:
        s = GroSystem()
    except:
        # no GroSystem support
        return

    assert_true(s.isNull())
    assert_true(s.isEmpty())
    assert_equal(s.nMolecules(), 0)

    s.setName("Test System")

    assert_equal(s.name(), "Test System")

    assert_false(s.isNull())
    assert_true(s.isEmpty())

    s.add("water")

    assert_false(s.isEmpty())
    assert_false(s.isNull())
    assert_equal(s.nMolecules(), 1)
    assert_equal(len(s.uniqueTypes()), 1)
    assert_equal(s.uniqueTypes()[0], "water")

    assert_equal(s[0], "water")

    s.add("water", 100)

    assert_equal(s.nMolecules(), 101)
    assert_equal(len(s.uniqueTypes()), 1)
    assert_equal(s.uniqueTypes()[0], "water")

    for i in range(0, s.nMolecules()):
        assert_equal(s[i], "water")

    s.add("sodium", 30)

    assert_equal(s.nMolecules(), 131)
    assert_equal(len(s.uniqueTypes()), 2)
    assert_equal(s.uniqueTypes()[0], "water")
    assert_equal(s.uniqueTypes()[1], "sodium")

    for i in range(0, s.nMolecules()):
        if i <= 100:
            assert_equal(s[i], "water")
        else:
            assert_equal(s[i], "sodium")


def test_pairs(verbose=False):

    # Load cyclohexane molecule and create System object.
    files = ["../io/cyclohexane.gro", "../io/cyclohexane.top"]
    s = MoleculeParser.read(files, {"GROMACS_PATH": gromacs_path})

    # Convert to a GroTop object and write to file.
    top = GroTop(s)
    top.writeToFile("cyclohexane.top")

    # Load validation toplogy file.
    lines0 = []
    with open("../io/cyclohexane_check.top", "r") as file:
        for line in enumerate(file):
            lines0.append(line)

    # Load written topology file.
    lines1 = []
    with open("cyclohexane.top", "r") as file:
        for line in enumerate(file):
            lines1.append(line)

    # Make sure files are the same size.
    assert len(lines0) == len(lines1)

    # Validate each record.
    for x in range(0, len(lines0)):
        # Ignore timestamp.
        if x != 1:
            assert lines0[x] == lines1[x]


if __name__ == "__main__":
    test_grotop(True)
    test_grosys(True)
    test_pairs(True)
