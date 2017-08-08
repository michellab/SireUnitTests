
from Sire.IO import *
from Sire.Mol import *

from nose.tools import assert_equal, assert_almost_equal

# check that we have Gro87 support in this version of Sire
has_gro87 = True

try:
    p = Gro87()
except:
    # No Gro87 support
    has_gro87 = False

def test_gro87(verbose=False):
    if not has_gro87:
        return

    grofile = "../io/water.gro"
    grofile = "../io/urea.gro"

    if verbose:
        print("Reading the Gro87 file using the Gro87 parser...")

    gro87 = Gro87(grofile)

    if verbose:
        print(gro87)

if __name__ == "__main__":
    test_gro87(True)
