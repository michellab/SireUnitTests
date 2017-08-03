
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

    if verbose:
        print("Reading the file using the old PDB parser...")

    oldmols = PDB().read(pdbfile)

    if verbose:
        print("Reading the file using the new PDB parser...")

    newmols = MoleculeParser.read(pdbfile)

    # make sure the same number of molecules have been written
    assert_equal( oldmols.nMolecules(), newmols.nMolecules() )

    # now compare the protein loaded from the file by both parsers
    oldprot = oldmols[MolWithResID("ALA")]
    newprot = newmols[MolWithResID("ALA")]

    # write some tests that compare these two proteins
    assert_equal( oldprot.nAtoms(), newprot.nAtoms() )


if __name__ == "__main__":
    test_pdb2(True)
