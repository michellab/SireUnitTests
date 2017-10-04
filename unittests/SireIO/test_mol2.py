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

# Specific atom coordinate data validation test for file "../io/complex_6.mol2".
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
    p = Mol2("../io/complex_6.mol2")

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

if __name__ == "__main__":
    test_read(True)
    test_atom_coords(True)
