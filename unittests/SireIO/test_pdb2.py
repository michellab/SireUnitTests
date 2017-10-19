from Sire.Base import *
from Sire.IO import *
from Sire.Mol import *

from glob import glob
from nose.tools import assert_equal, assert_almost_equal

# check that we have PDB2 support in this version of Sire
has_pdb2 = True

try:
    p = PDB2()
except:
    # No PDB2 support
    has_pdb2 = False

# General test of ability to read PDB files.
# All PDB files in the "../io/" directory are parsed and data
# is validated for all files containing a MASTER record.
def test_read(verbose=False):
    if not has_pdb2:
        return

    # Glob all of the PDB files in the example file directory.
    pdbfiles = glob('../io/*pdb')

    for file in pdbfiles:
        # Parse the file into a PDB2 object.
        # Errors should be thrown if the record data in a file
        # doesn't match the PDB format.
        p = PDB2(file)

        # If there is a master record for this file, then validate
        # the parsed data against it.
        if p.hasMaster():
            # Extract the master record.
            m = p.getMaster()

            # Extract the title record.
            t = p.getTitle();

            # Work out the number of coordinate transformation records.
            # A complete object counts as 3 records, i.e. 1 for each dimension.
            num_transform = 3 * ( p.hasTransOrig() + p.hasTransScale() + p.hasTransMatrix() )

            # Validate data.
            assert_equal( t.nRemarks(), m.nRemarks() )
            assert_equal( p.nAtoms(), m.nAtoms() )
            assert_equal( p.nHelices(), m.nHelices() )
            assert_equal( p.nSheets(), m.nSheets() )
            assert_equal( p.nTers(), m.nTers() )
            assert_equal( num_transform, m.nTransforms() )

# Specific atom coordinate data validation test for file "../io/ntrc.pdb".
def test_atom_coords(verbose=False):
    if not has_pdb2:
        return

    # Test atoms.
    atoms = ["CA", "CB", "N", "O", "HB"]

    # Test coordinates.
    coords = [[-13.721,  -3.484, 14.690],
              [-10.695,  -0.294, 14.759],
              [ -8.536,  -2.557, 13.277],
              [ -7.037,  -1.615,  9.350],
              [ -5.045,   2.118,  8.812]]

    # Parse the PDB file.
    p = PDB2("../io/ntrc.pdb")

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

# Test that all files the parser can convert a Sire molecule
# back into the correct data format, ready to be written to file.
def test_write(verbose=False):
    if not has_pdb2:
        return

    # Glob all of the PDB files in the example file directory.
    pdbfiles = glob('../io/*pdb')

    for file in pdbfiles:
        print(file)
        # Parse the file into a PDB object.
        # Errors should be thrown if the record data in a file
        # doesn't match the Mol2 format.
        p = PDB2(file)

        # Construct a Sire molecular system.
        s = p.toSystem()

        # Now re-parse the molecular system.
        p = PDB2(s)

if __name__ == "__main__":
    test_read(True)
    test_write(True)
    test_atom_coords(True)
