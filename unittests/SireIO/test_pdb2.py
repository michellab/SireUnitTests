from Sire.IO import *
from Sire.Mol import *

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

    # These are useful test files since they contain a large assortment of PDB
    # record types, including a MASTER record for data validation.
    pdbfiles = [ "../io/ntrc.pdb", "../io/1P38.pdb", "../io/dioxin.pdb" ]

    for file in pdbfiles:
        # Parse the file into a PDB2 object.
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
