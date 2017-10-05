from Sire.Base import *
from Sire.IO import *
from Sire.Mol import *

from glob import glob
from nose.tools import assert_equal, assert_almost_equal

# check that we have Mol2 support in this version of Sire
has_mol2 = True

try:
    p = Mol2()
except:
    # No Mol2 support
    has_mol2 = False

# General test of ability to read Mol2 files.
# All Mol2 files in the "../io/" directory are parsed.
def test_read(verbose=False):
    if not has_mol2:
        return

    # Glob all of the Mol2 files in the example file directory.
    mol2files = glob('../io/*mol2')

    for file in mol2files:
        # Parse the file into a Mol2 object.
        # Errors should be thrown if the record data in a file
        # doesn't match the Mol2 format.
        p = Mol2(file)

# Specific atom coordinate data validation test for file "../io/complex.mol2".
def test_atom_coords(verbose=False):
    if not has_mol2:
        return

    # Test atoms.
    atoms = ["N", "CA", "C", "O", "CB"]

    # Test coordinates.
    coords = [[ -2.9880,  -2.0590,  -2.6220],
              [ -3.8400,  -2.0910,  -7.4260],
              [ -6.4250,  -3.9190, -10.9580],
              [ -6.1980,  -6.0090, -14.2910],
              [ -9.8700,  -6.5500, -15.2480]]

    # Parse the Mol2 file.
    p = Mol2("../io/complex.mol2")

    # Create a molecular system.
    s = p.toSystem()

    # Get the first molecule.
    m = s[MolIdx(0)]

    # Loop over all of the atoms.
    for i in range(0, len(atoms)):

        # Extract the atom from the residue "i + 1".
        a = m.atom(AtomName(atoms[i]) + ResNum(i+1))

        # Extract the atom coordinates.
        c = a.property("coordinates")

        # Validate parsed coordinates against known values.
        assert_almost_equal( c[0], coords[i][0] )
        assert_almost_equal( c[1], coords[i][1] )
        assert_equal( c[2], coords[i][2] )

# Residue and chain validation test for file "../io/complex.mol2".
def test_residues(verbose=False):
    if not has_mol2:
        return

    # Parse the Mol2 file.
    p = Mol2("../io/complex.mol2")

    # Create a molecular system.
    s = p.toSystem()

    # Get the two molecules.
    m1 = s[MolIdx(0)]
    m2 = s[MolIdx(1)]

    # Get the chains from the molecules.
    c1 = m1.chains()
    c2 = m2.chains()

    # Check the number of chains in each molecule.
    assert_equal( len(c1), 3 )
    assert_equal( len(c2), 1 )

    # Check the number of residues in each chain of the first molecule.
    assert_equal( len(c1[0].residues()), 118 )
    assert_equal( len(c1[1].residues()), 114 )
    assert_equal( len(c1[2].residues()), 118 )

    # Check the number of residues in the single chain of the second molecule.
    assert_equal( len(c2[0].residues()), 1 )

    # Check some specific residue names in the first chain from the first molecule.
    assert_equal( c1[0].residues()[0].name().toString(), "ResName('PRO1')" )
    assert_equal( c1[1].residues()[1].name().toString(), "ResName('MET2')" )
    assert_equal( c1[1].residues()[2].name().toString(), "ResName('PHE3')" )

if __name__ == "__main__":
    test_read(True)
    test_atom_coords(True)
    test_residues(True)
