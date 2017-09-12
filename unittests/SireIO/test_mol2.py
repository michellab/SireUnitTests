from Sire.IO import *
from Sire.Mol import *

from nose.tools import assert_equal, assert_almost_equal

# check that we have Mol2 support in this version of Sire
has_mol2 = True

try:
    p = Mol2()
except:
    # No Mol2 support
    has_mol2 = False

def test_mol2(verbose=False):
    if not has_mol2:
        return

    mol2file = "../io/dioxin.mol2"

    m = Mol2(mol2file);

if __name__ == "__main__":
    test_mol2(True)
