try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.Mol import *
from Sire.IO import *
from Sire.ID import *
from Sire.Maths import *

import sys

from nose.tools import assert_equal, assert_true, assert_false

mol = PDB().readMolecule("../io/p38.pdb")

rangen = RanGenerator()


def test_select(verbose=False):
    s = mol.selection()

    if verbose:
        print("Test select all")

    s = s.selectAll()
    assert_true(s.selectedAll())
    assert_false(s.selectedNone())
    assert_equal(s.nSelected(), mol.nAtoms())
    assert_true(s.selectedAllCutGroups())
    assert_true(s.selectedAllAtoms())
    assert_true(s.selectedAllResidues())
    assert_true(s.selectedAllChains())
    assert_true(s.selectedAllSegments())

    for i in range(0, mol.nAtoms()):
        assert_true(s.selected(AtomIdx(i)))

    if verbose:
        print("Test select none")

    s = s.deselectAll()
    assert_false(s.selectedAll())
    assert_true(s.selectedNone())
    assert_equal(s.nSelected(), 0)
    assert_false(s.selectedAllCutGroups())
    assert_false(s.selectedAllAtoms())
    assert_false(s.selectedAllResidues())
    assert_false(s.selectedAllChains())
    assert_false(s.selectedAllSegments())

    for i in range(0, mol.nAtoms()):
        assert_false(s.selected(AtomIdx(i)))

    if verbose:
        print("Test select one atom")

    for i in range(0, mol.nAtoms()):
        if verbose and i % 100 == 0:
            print(".", end="")
            sys.stdout.flush()

        s = s.selectNone()

        s = s.select(AtomIdx(i))

        assert_false(s.selectedAll())
        assert_false(s.selectedNone())
        assert_equal(s.nSelected(), 1)
        assert_false(s.selectedAllCutGroups())
        assert_false(s.selectedAllAtoms())
        assert_false(s.selectedAllResidues())
        assert_true(s.selected(AtomIdx(i)))

    if verbose:
        print("...done")
        print("Test select one residue")

    for i in range(0, mol.nResidues()):
        if verbose and i % 10 == 0:
            print(".", end="")
            sys.stdout.flush()

        res = mol.residue(ResIdx(i))

        s = s.selectNone()
        s = s.select(ResIdx(i))

        assert_false(s.selectedAll())
        assert_false(s.selectedNone())
        assert_equal(s.nSelected(), res.nAtoms())
        assert_false(s.selectedAllCutGroups())
        assert_false(s.selectedAllAtoms())
        assert_false(s.selectedAllResidues())
        assert_true(s.selected(ResIdx(i)))
        assert_true(s.selectedAll(ResIdx(i)))

        for atom in res.atoms():
            assert_true(s.selected(atom.index()))

    if verbose:
        print("...done")


def test_multi_select(verbose=False):
    s = mol.selection()

    if verbose:
        print("\nTest multi select...", end="")
        sys.stdout.flush()

    s = s.selectNone()

    # generate random sets of atoms and make
    # sure that the selection is correct
    for i in range(0, 100):
        if verbose and i % 10 == 0:
            print(".", end="")
            sys.stdout.flush()

        s = s.selectNone()
        assert_true(s.selectedNone())

        natms = rangen.randInt(2, 50)

        selected = []
        cgs = []

        for j in range(0, natms):
            idx = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

            while idx in selected:
                assert_true(s.selected(idx))
                idx = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

            cgidx = mol.info().parentCutGroup(idx)

            assert_false(s.selected(idx))

            if cgidx in cgs:
                assert_true(s.selected(cgidx))
            else:
                assert_false(s.selected(cgidx))
                cgs.append(cgidx)

            selected.append(idx)
            s = s.select(idx)
            assert_equal(s.nSelected(), j + 1)
            assert_true(s.selected(cgidx))
            assert_equal(s.nSelectedCutGroups(), len(cgs))

            t = s.selectedAtoms()
            assert_equal(len(selected), len(t))
            for atom in t:
                assert_true(atom in selected)

            t = s.selectedCutGroups()
            assert_equal(len(cgs), len(t))
            for cg in t:
                assert_true(cg in cgs)

        for j in range(0, natms):
            s = s.deselect(selected[j])
            assert_false(s.selected(selected[j]))
            assert_equal(s.nSelected(), natms - j - 1)

        assert_true(s.selectedNone())

    if verbose:
        print("...done")


def test_intersect(verbose=False):
    if verbose:
        print("Testing intersections/masking...", end="")
        sys.stdout.flush()

    s0 = mol.selection()
    s0 = s0.selectNone()
    s1 = mol.selection()
    s1 = s1.selectNone()

    assert_true(s0.selectedNone())
    assert_true(s1.selectedNone())

    # create two selections and make sure that
    # intersect, subtract and unite all work
    for i in range(0, 100):
        if verbose and i % 10 == 0:
            print(".", end="")
            sys.stdout.flush()

        s0 = s0.selectNone()
        assert_true(s0.selectedNone())
        s1 = s1.selectNone()
        assert_true(s1.selectedNone())

        natms = rangen.randInt(1, 50)

        for j in range(0, natms):
            idx = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

            while s0.selected(idx):
                idx = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

            s0.select(idx)
            assert_true(s0.selected(idx))
            assert_equal(s0.nSelected(), j + 1)

        natms = rangen.randInt(1, 50)

        for j in range(0, natms):
            idx = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

            while s1.selected(idx):
                idx = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

            s1 = s1.select(idx)
            assert_true(s1.selected(idx))
            assert_equal(s1.nSelected(), j + 1)

        natms = rangen.randInt(1, 50)

        # make sure that there are some overlapping atoms
        for j in range(0, natms):
            idx = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))
            s0 = s0.select(idx)
            s1 = s1.select(idx)

        # test intersect/mask
        t = AtomSelection(s0)
        t = t.mask(s1)

        for atom in t.selectedAtoms():
            assert_true(s0.selected(atom) and s1.selected(atom))

        for atom in s0.selectedAtoms():
            if s1.selected(atom):
                assert_true(t.selected(atom))

        for cg in t.selectedCutGroups():
            assert_true(s0.selected(cg) and s1.selected(cg))

        # test subtraction
        t = AtomSelection(s0)
        t = t.subtract(s1)

        for atom in t.selectedAtoms():
            assert_false(s1.selected(atom))

        for atom in s1.selectedAtoms():
            assert_false(t.selected(atom))

        t = AtomSelection(s1)
        t = t.subtract(s0)

        for atom in t.selectedAtoms():
            assert_false(s0.selected(atom))

        for atom in s0.selectedAtoms():
            assert_false(t.selected(atom))

        # test unite
        t = AtomSelection(s0)
        t = t.unite(s1)

        for atom in t.selectedAtoms():
            assert_true(s0.selected(atom) or s1.selected(atom))

        for atom in s0.selectedAtoms():
            assert_true(t.selected(atom))

        for atom in s1.selectedAtoms():
            assert_true(t.selected(atom))

        for cg in t.selectedCutGroups():
            assert_true(s0.selected(cg) or s1.selected(cg))

        for cg in s0.selectedCutGroups():
            assert_true(t.selected(cg))

        for cg in s1.selectedCutGroups():
            assert_true(t.selected(cg))

    if verbose:
        print("...done")


def test_deselect(verbose=False):
    if verbose:
        print("Testing deselecting from full set...", end="")
        sys.stdout.flush()

    s = mol.selection()
    s = s.selectAll()

    assert_true(s.selectedAll())

    # test deselecting one atom
    for i in range(0, 40):
        if verbose:
            print(".", end="")
            sys.stdout.flush()

        s = s.selectAll()
        assert_true(s.selectedAll())
        assert_equal(s.nSelected(), mol.nAtoms())

        idx = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

        s = s.deselect(idx)

        assert_false(s.selected(idx))

        assert_equal(s.nSelected(), mol.nAtoms() - 1)

        for j in range(0, 40):
            idx2 = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

            while idx2 == idx:
                idx2 = AtomIdx(rangen.randInt(0, mol.nAtoms() - 1))

            s = s.deselect(idx2)

            assert_false(s.selected(idx))
            assert_false(s.selected(idx2))
            assert_equal(s.nSelected(), mol.nAtoms() - 2)

            s = s.select(idx2)

            assert_false(s.selected(idx))
            assert_true(s.selected(idx2))
            assert_equal(s.nSelected(), mol.nAtoms() - 1)

    for i in range(0, 40):
        s = s.selectAll()
        assert_true(s.selectedAll())
        assert_equal(s.nSelected(), mol.nAtoms())

        idx = CGIdx(rangen.randInt(0, mol.nCutGroups() - 1))

        s = s.deselect(idx)

        assert_false(s.selected(idx))
        assert_equal(s.nSelected(), mol.nAtoms() - mol.info().nAtoms(idx))

        for j in range(0, 40):
            idx2 = CGIdx(rangen.randInt(0, mol.nCutGroups() - 1))

            while idx2 == idx:
                idx2 = CGIdx(rangen.randInt(0, mol.nCutGroups() - 1))

            s = s.deselect(idx2)

            assert_false(s.selected(idx))
            assert_false(s.selected(idx2))
            assert_equal(
                s.nSelected(),
                mol.nAtoms()
                - mol.info().nAtoms(idx)
                - mol.info().nAtoms(idx2),
            )

            s = s.select(idx2)

            assert_false(s.selected(idx))
            assert_true(s.selected(idx2))
            assert_equal(s.nSelected(), mol.nAtoms() - mol.info().nAtoms(idx))

    if verbose:
        print("...done")


def test_invert(verbose=False):
    if verbose:
        print("Testing inversion...", end="")
        sys.stdout.flush()

    s0 = mol.selection()
    s0 = s0.selectAll()
    assert_true(s0.selectedAll())

    for i in range(0, 50):
        if verbose and i % 5 == 0:
            print(".", end="")
            sys.stdout.flush()

        s1 = mol.selection()
        s1 = s1.deselectAll()
        s0 = s0.selectAll()

        # generate a random set to subtract
        for j in range(0, rangen.randInt(1, 300)):
            s1 = s1.select(AtomIdx(rangen.randInt(0, mol.nAtoms() - 1)))

        s0 = s0.subtract(s1)

        assert_equal(s0.nSelected(), mol.nAtoms() - s1.nSelected())

        for atom in s0.selectedAtoms():
            assert_false(s1.selected(atom))

        s0 = s0.invert()

        assert_equal(s0.nSelected(), s1.nSelected())

        for atom in s0.selectedAtoms():
            assert_true(s1.selected(atom))

        s0 = s0.invert()

        assert_equal(s0.nSelected(), mol.nAtoms() - s1.nSelected())

        for atom in s0.selectedAtoms():
            assert_false(s1.selected(atom))

    if verbose:
        print("...done!")


if __name__ == "__main__":
    test_select(True)
    test_multi_select(True)
    test_intersect(True)
    test_deselect(True)
    test_invert(True)
