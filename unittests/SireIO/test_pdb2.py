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

    pdbfile = "../io/ntrc.pdb"

    p = PDB2(pdbfile)

    # Testing that we're reading data correctly.
    # Currently reading ATOM, HELIX, and SHEET records ok.
    assert_equal( p.num_atom, 1910 )
    assert_equal( p.num_helix, 4 )
    assert_equal( p.num_sheet, 3 )

if __name__ == "__main__":
    test_pdb2(True)
