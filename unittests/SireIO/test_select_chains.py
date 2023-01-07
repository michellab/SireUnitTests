try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

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

# Regression test for bug in AtomSelection::selectedChains()
def test_chains(verbose=False):
    if not has_mol2:
        return

    # Parse a Mol2 file containing multiple chains.
    p = Mol2("../io/complex.mol2")

    # Create a Sire molecular system.
    s = p.toSystem()

    # Extract the first molecule.
    m = s[MolIdx(0)]

    # Extract the residues from the molecule.
    r = m.residues()

    # Assert that the residues belong to the correct chains.
    assert_equal(r[0].chain().name().value(), "A")
    assert_equal(r[-1].chain().name().value(), "C")


if __name__ == "__main__":
    test_chains(True)
