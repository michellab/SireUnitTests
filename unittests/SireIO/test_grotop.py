
from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *
from Sire.Base import *

from nose.tools import assert_equal, assert_almost_equal

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

    g = GroTop("../io/urea.top", {"GROMACS_PATH":gromacs_path})
    #g = GroTop("../io/dppc.top", {"GROMACS_PATH":gromacs_path})

    print(g)
    #print(g.atomTypes())
    #print(g.bondPotentials())
    #print(g.anglePotentials())
    #print(g.dihedralPotentials())

    for mol in g.moleculeTypes():
        print("MOL: %s" % mol)

        for atom in mol.atoms():
            print("ATOM: %s" % atom)

        for bond, param in mol.bonds().items():
            print("BOND: %s | %s" % (bond,param))

        for angle, param in mol.angles().items():
            print("ANGLE: %s | %s" % (angle,param))

        for dihedral, param in mol.dihedrals().items():
            print("DIHEDRAL: %s | %s" % (dihedral,param))

if __name__ == "__main__":
    test_grotop(True)
