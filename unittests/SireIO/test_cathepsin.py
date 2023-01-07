try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.System import *
from Sire.MM import *
from Sire.CAS import *

import glob
import os

from nose.tools import assert_equal

d = "../io/cathepsin"

prms = glob.glob("%s/*.parm7" % d)

inputs = {}

for prm in prms:
    root = prm[0:-6]
    rst = glob.glob("%s*.rst7" % root)[0]

    inputs[prm] = rst


def _get_first_molecules(s):
    m = MoleculeGroup("all")
    intraff = InternalFF("intraff")
    intraclj = IntraCLJFF("intraclj")

    for i in range(0, s.nMolecules()):
        mol = s[MolIdx(i)]

        # only add non-water molecules
        if mol.nAtoms() > 5:
            m.add(mol)
            intraff.add(mol)
            intraclj.add(mol)

    s = System()
    s.add(m)
    s.add(intraff)
    s.add(intraclj)

    return s


def _combine_dih_imp(nrgs):
    dih = Symbol("E_{intraff}^{dihedral}")
    imp = Symbol("E_{intraff}^{improper}")
    nrgs[dih] += nrgs[imp]
    nrgs[imp] = 0


def _get_energies(s):
    nrgs = s.energies()
    n = {}

    for key in nrgs.keys():
        n[key] = nrgs[key]

    _combine_dih_imp(n)
    return n


def _print_energies(nrgs1, nrgs2):
    keys = list(nrgs1.keys())
    keys.sort()

    for key in keys:
        print("%s  %s  %s" % (key, nrgs1[key], nrgs2[key]))


def _compare_energies(nrgs1, nrgs2):
    keys = list(nrgs1.keys())
    keys.sort()

    for key in keys:
        diff = abs(nrgs1[key] - nrgs2[key])
        if diff > 0.1:
            print(key)
            assert_equal(nrgs1[key], nrgs2[key])


def _test_input(prm, crd, verbose=False, slow_tests=False):
    root = prm.split("/")[-1][0:-6]

    if verbose:
        print("\n************\nTesting %s (%s | %s)" % (root, prm, crd))
        print("\nTesting straight read and write")

    r = AmberRst7(crd)

    a = AmberPrm(prm)

    d = os.path.dirname("test-%s" % root)
    if not os.path.isdir(d) and d != "":
        os.mkdir(d)

    a.writeToFile("test-%s.prm" % root)

    a2 = AmberPrm("test-%s.prm" % root)

    if verbose:
        print(a, a2)

    assert_equal(a.nAtoms(), a2.nAtoms())

    full_s = a.toSystem(r)

    s = _get_first_molecules(full_s)
    nrgs = _get_energies(s)

    r = AmberRst7(s)
    a = AmberPrm(s)

    if verbose:
        print(a)

    if slow_tests:
        if verbose:
            print("\nTesting triangle conversion to grotop")

        g = GroTop(s)
        g.writeToFile("test-%s.grotop" % root)

        if verbose:
            print(g)

        g = GroTop("test-%s.grotop" % root)
        s2 = _get_first_molecules(g.toSystem(r))

        nrgs2 = _get_energies(s2)

        if verbose:
            _print_energies(nrgs, nrgs2)

        _compare_energies(nrgs, nrgs2)

    if verbose:
        print("\nComparing backwards/forwards conversion")

    a2 = AmberPrm(s)

    if verbose:
        print(a2)

    a2.writeToFile("test-%s.prm7" % root)
    a2 = AmberPrm("test-%s.prm7" % root)

    s2 = _get_first_molecules(a2.toSystem(r))

    nrgs2 = _get_energies(s2)

    if verbose:
        _print_energies(nrgs, nrgs2)

    _compare_energies(nrgs, nrgs2)

    if verbose:
        print("Comparing writing to coordinate files")

    r = AmberRst(full_s)
    r.writeToFile("test-%s.rst" % root)

    s2 = MoleculeParser.read(prm, "test-%s.rst" % root)

    s2 = _get_first_molecules(s2)

    nrgs2 = _get_energies(s2)

    if verbose:
        _print_energies(nrgs, nrgs2)

    _compare_energies(nrgs, nrgs2)


def test_cathepsin(verbose=False, slow_tests=False):
    keys = list(inputs.keys())
    keys.sort()

    for prm in keys:
        _test_input(prm, inputs[prm], verbose, slow_tests)


if __name__ == "__main__":
    test_cathepsin(True, slow_tests=False)
