
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

    #g = GroTop("../io/urea.top", {"GROMACS_PATH":gromacs_path})
    g = GroTop("../io/dppc.top", {"GROMACS_PATH":gromacs_path})

    print(g.warnings())

    print(g)

    print(g.groSystem())

    #print(g.atomTypes())
    #print(g.bondPotentials())
    #print(g.anglePotentials())
    #print(g.dihedralPotentials())

    for mol in g.moleculeTypes():
        print("MOL: %s" % mol)

        for atom in mol.atoms():
            print("ATOM: %s" % atom)

        #for bond, param in mol.bonds().items():
        #    print("BOND: %s | %s" % (bond,param))

        #for angle, param in mol.angles().items():
        #    print("ANGLE: %s | %s" % (angle,param))

        #for dihedral, param in mol.dihedrals().items():
        #    print("DIHEDRAL: %s | %s" % (dihedral,param))

def test_grosys(verbose=False):

    try:
        s = GroSystem()
    except:
        # no GroSystem support
        return

    assert_true( s.isNull() )
    assert_true( s.isEmpty() )
    assert_equal( s.nMolecules(), 0 )

    s.setName("Test System")

    assert_equal( s.name(), "Test System" )

    assert_false( s.isNull() )
    assert_true( s.isEmpty() )

    s.add("water")

    assert_false( s.isEmpty() )
    assert_false( s.isNull() )
    assert_equal( s.nMolecules(), 1 )
    assert_equal( len(s.uniqueTypes()), 1 )
    assert_equal( s.uniqueTypes()[0], "water" )

    assert_equal( s[0], "water" )

    s.add("water", 100)

    assert_equal( s.nMolecules(), 101 )
    assert_equal( len(s.uniqueTypes()), 1 )
    assert_equal( s.uniqueTypes()[0], "water" )

    for i in range(0, s.nMolecules()):
        assert_equal( s[i], "water" )

    s.add("sodium", 30)

    assert_equal( s.nMolecules(), 131 )
    assert_equal( len(s.uniqueTypes()), 2 )
    assert_equal( s.uniqueTypes()[0], "water" )
    assert_equal( s.uniqueTypes()[1], "sodium" )

    for i in range(0, s.nMolecules()):
        if i <= 100:
            assert_equal( s[i], "water" )
        else:
            assert_equal( s[i], "sodium" )

if __name__ == "__main__":
    test_grotop(True)
    test_grosys(True)
