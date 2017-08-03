
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

    pdbfile = "../io/dioxin.pdb"
    mol2file = "../io/dioxin.mol2"

    if verbose:
        print("Reading the PDB file using the PDB parser...")

    pdbmols = MoleculeParser.read(pdbfile)

    if verbose:
        print("Reading the Mol2 file using the Mol2 parser...")

    mol2mols = MoleculeParser.read(mol2file)

    # make sure the same number of molecules have been written
    assert_equal( pdbmols.nMolecules(), mol2mols.nMolecules() )

    # now compare the first molecule loaded from the file by both parsers
    pdbmol = oldmols[MolIdx(0)]
    mol2mol = newmols[MolIdx(0)]

    # write some tests that compare these two proteins
    assert_equal( pdbmol.nAtoms(), mol2mol.nAtoms() )


if __name__ == "__main__":
    test_mol2(True)
