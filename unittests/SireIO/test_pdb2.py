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

    # This is a useful test file since it contains a large
    # assortment of PDB records.
    pdbfile = "../io/ntrc.pdb"

    # Parse the file into a PDB2 object.
    p = PDB2(pdbfile)

    # First check that we're parsing the various records correctly.
    assert_equal( p.nTitles(), 119 )
    assert_equal( p.nAtoms(), 1910 )
    assert_equal( p.nHelices(), 4 )
    assert_equal( p.nSheets(), 3 )
    assert( p.hasCrystal() )
    assert( p.hasTransOrig() )
    assert( p.hasTransScale() )

    # Now validate data for specifc objects...

if __name__ == "__main__":
    test_pdb2(True)
