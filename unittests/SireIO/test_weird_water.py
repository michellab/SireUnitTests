try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *

import Sire.Stream

import os, re, sys
import shutil

from nose.tools import assert_equal, assert_almost_equal

combining_rules = "arithmetic"
temperature = 25 * celsius
pressure = 1 * atm
coulomb_cutoff = 1000 * angstrom
coulomb_feather = 999.5 * angstrom
lj_cutoff = 1000 * angstrom
lj_feather = 999.5 * angstrom
#############################################################


def _assert_same_function(f1, f2, x, start, end):
    delta = (end - start) / 25

    for i in range(0, 26):
        val = {x: (start + i * delta)}
        v1 = f1.evaluate(val)
        v2 = f2.evaluate(val)

        if abs(v2 - v1) > 0.001:
            assert_equal(f1, f2)


def _pvt_compare_molecules(mol1, mol2, verbose):
    # compare the residues
    assert_equal(mol1.nResidues(), mol2.nResidues())

    for i in range(0, mol1.nResidues()):
        res1 = mol1.residue(ResIdx(i))
        res2 = mol2.residue(ResIdx(i))

        assert_equal(res1.name(), res2.name())
        assert_equal(res1.number(), res2.number())

    # compare the atoms
    assert_equal(mol1.nAtoms(), mol2.nAtoms())

    for i in range(0, mol1.nAtoms()):
        atom1 = mol1.atom(AtomIdx(i))
        atom2 = mol2.atom(AtomIdx(i))

        assert_equal(atom1.name(), atom2.name())
        assert_equal(atom1.number(), atom2.number())

        for p in [
            "coordinates",
            "charge",
            "LJ",
            "mass",
            "element",
            "ambertype",
        ]:
            p1 = atom1.property(p)
            p2 = atom2.property(p)
            assert_equal(p1, p2)

    if verbose:
        print("Compared nAtoms = %s : all equal" % mol1.nAtoms())

    # compare the connectivity
    try:
        conn1 = mol1.property("connectivity")
        conn2 = mol2.property("connectivity")
        have_connectivity = True
    except:
        have_connectivity = False

    if have_connectivity:
        assert_equal(conn1.nConnections(), conn2.nConnections())

        for i in range(0, mol1.nAtoms()):
            idx = AtomIdx(i)

            assert_equal(conn1.nConnections(idx), conn2.nConnections(idx))

            bonded1 = conn1.connectionsTo(idx)
            bonded2 = conn2.connectionsTo(idx)

            assert_equal(len(bonded1), len(bonded2))

            for atom in bonded1:
                assert atom in bonded2

        if verbose:
            print("Connectivity is equal")

    # compare the bond, angle, dihedral and improper functions
    try:
        bonds1 = mol1.property("bond")
        bonds2 = mol2.property("bond")
        have_bonds = True
    except:
        have_bonds = False

    if have_bonds:
        assert_equal(bonds1.nFunctions(), bonds2.nFunctions())

        for func in bonds1.potentials():
            p1 = func.function()
            p2 = bonds2.potential(func.atom0(), func.atom1())
            assert_equal(p1, p2)

        if verbose:
            print("Compared nBonds = %s : all equal" % bonds1.nFunctions())

    try:
        angles1 = mol1.property("angle")
        angles2 = mol2.property("angle")
        have_angles = True
    except:
        have_angles = False

    if have_angles:
        assert_equal(angles1.nFunctions(), angles2.nFunctions())

        for func in angles1.potentials():
            p1 = func.function()
            p2 = angles2.potential(func.atom0(), func.atom1(), func.atom2())
            assert_equal(p1, p2)

        if verbose:
            print("Compared nAngles = %s : all equal" % angles1.nFunctions())

    try:
        dihedrals1 = mol1.property("dihedral")
        dihedrals2 = mol2.property("dihedral")
        have_dihs = True
    except:
        have_dihs = False

    if have_dihs:
        assert_equal(dihedrals1.nFunctions(), dihedrals2.nFunctions())

        for func in dihedrals1.potentials():
            p1 = func.function()
            p2 = dihedrals2.potential(
                func.atom0(), func.atom1(), func.atom2(), func.atom3()
            )

            _assert_same_function(p1, p2, Symbol("phi"), 0, 3.141)

        if verbose:
            print(
                "Compared nDihedrals = %s : all equal"
                % dihedrals1.nFunctions()
            )

    try:
        impropers1 = mol1.property("improper")
        impropers2 = mol2.property("improper")
        have_imps = True
    except:
        have_imps = False

    if have_imps:
        assert_equal(impropers1.nFunctions(), impropers2.nFunctions())

        for func in impropers1.potentials():
            p1 = func.function()
            p2 = impropers2.potential(
                func.atom0(), func.atom1(), func.atom2(), func.atom3()
            )

            _assert_same_function(p1, p2, Symbol("phi"), 0, 3.141)

        if verbose:
            print(
                "Compared nImpropers = %s : all equal"
                % impropers1.nFunctions()
            )

    if have_connectivity:
        # now calculate an intramolecular energy
        if verbose:
            print("Compare IntraCLJ energy...")

        ff1 = IntraFF()
        ff2 = IntraFF()

        ff1.add(mol1)
        ff2.add(mol2)

        nrg1 = ff1.energy().value()
        nrg2 = ff2.energy().value()

        if verbose:
            print("%s versus %s" % (nrg1, nrg2))

        assert_almost_equal(nrg1, nrg2)

        if verbose:
            print("Compare Internal bond/angle/dihedral energy...")

        ff1 = InternalFF()
        ff2 = InternalFF()

        ff1.add(mol1)
        ff2.add(mol2)

        nrg1 = ff1.energy().value()
        nrg2 = ff2.energy().value()

        if verbose:
            print("%s versus %s" % (nrg1, nrg2))

        assert_almost_equal(nrg1, nrg2)


def test_weird_water(verbose=False):
    try:
        # check if we have this
        r = AmberRst()
    except:
        return

    rst_file = "../io/proteinbox.crd"
    top_file = "../io/proteinbox.top"

    # Load the box of molecules using the old and new parser
    if verbose:
        print("Loading box of molecules using old parser...")

    (mols1, space1) = Amber().readCrdTop(rst_file, top_file)

    if verbose:
        print("Loading box of molecules using new parser...")

    system2 = MoleculeParser.read(rst_file, top_file)

    if verbose:
        print("Comparing all molecules in the boxes...")

    mols2 = system2[MGIdx(0)]
    space2 = system2.property("space")

    assert_equal(space1, space2)
    assert_equal(mols1.nMolecules(), mols2.nMolecules())

    # this water causes problems
    water1 = mols1[MolIdx(333)].molecule()
    water2 = mols2[MolIdx(333)].molecule()

    try:
        _pvt_compare_molecules(water1, water2, verbose)
    except:
        Sire.Stream.save([water1, water2], "broken.s3")
        raise


if __name__ == "__main__":
    test_weird_water(True)
