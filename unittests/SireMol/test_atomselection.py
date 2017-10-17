
from Sire.Mol import *
from Sire.IO import *
from Sire.ID import *

import sys

from nose.tools import assert_equal, assert_true, assert_false

mol = PDB().readMolecule("../io/p38.pdb")

def test_select(verbose=False):
    s = mol.selection()

    if verbose:
        print("Test select all")

    s = s.selectAll()
    assert_true( s.selectedAll() )
    assert_false( s.selectedNone() )
    assert_equal( s.nSelected(), mol.nAtoms() )
    assert_true( s.selectedAllCutGroups() )
    assert_true( s.selectedAllAtoms() )
    assert_true( s.selectedAllResidues() )
    assert_true( s.selectedAllChains() )
    assert_true( s.selectedAllSegments() )

    if verbose:
        print("Test select none")

    s = s.deselectAll()
    assert_false( s.selectedAll() )
    assert_true( s.selectedNone() )
    assert_equal( s.nSelected(), 0 )
    assert_false( s.selectedAllCutGroups() )
    assert_false( s.selectedAllAtoms() )
    assert_false( s.selectedAllResidues() )
    assert_false( s.selectedAllChains() )
    assert_false( s.selectedAllSegments() )

    if verbose:
        print("Test select one atom")

    for i in range(0, mol.nAtoms()):
        if verbose and i % 100 == 0:
            print(".", end="")
            sys.stdout.flush()

        s = s.selectNone()

        s = s.select( AtomIdx(i) )

        assert_false( s.selectedAll() )
        assert_false( s.selectedNone() )
        assert_equal( s.nSelected(), 1 )
        assert_false( s.selectedAllCutGroups() )
        assert_false( s.selectedAllAtoms() )
        assert_false( s.selectedAllResidues() )
        assert_true( s.selected(AtomIdx(i)) )

    if verbose:
        print("...done")
        print("Test select one residue")

    for i in range(0,mol.nResidues()):
        if verbose and i % 10 == 0:
            print(".", end="")
            sys.stdout.flush()

        res = mol.residue( ResIdx(i) )

        s = s.selectNone()
        s = s.select( ResIdx(i) )

        assert_false( s.selectedAll() )
        assert_false( s.selectedNone() )
        assert_equal( s.nSelected(), res.nAtoms() )
        assert_false( s.selectedAllCutGroups() )
        assert_false( s.selectedAllAtoms() )
        assert_false( s.selectedAllResidues() )
        assert_true( s.selected(ResIdx(i)) )
        assert_true( s.selectedAll(ResIdx(i)) )

        for atom in res.atoms():
            assert_true( s.selected(atom.index()) )

    if verbose:
        print("...done")

if __name__ == "__main__":
    test_select(True)

