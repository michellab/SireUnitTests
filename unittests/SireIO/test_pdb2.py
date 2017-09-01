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

def test_pdb2(verbose=False):
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

if __name__ == "__main__":
    test_pdb2(True)
