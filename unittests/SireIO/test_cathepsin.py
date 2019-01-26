
from Sire.IO import * 
from Sire.Mol import *
from Sire.System import *
from Sire.MM import *
from Sire.CAS import *

import glob

from nose.tools import assert_equal

d = "../io/cathepsin"

prms = glob.glob("%s/*.parm7" % d)

inputs = {}

for prm in prms:
    root = prm[0:-6]
    rst = glob.glob("%s*.rst7" % root)[0]

    inputs[prm] = rst

def _get_first_molecule(s):
    m = MoleculeGroup("all")
    mol = s[MolIdx(0)]

    intraff = InternalFF("intraff")
    intraclj = IntraCLJFF("intraclj")

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
        if diff > 0.001:
            print(key)
            assert_equal(nrgs1[key], nrgs2[key])

def _test_input(prm, crd, verbose=False):
    root = prm.split("/")[-1][0:-6]

    if verbose:
        print("\n************\nTesting %s (%s | %s)" % (root, prm, crd))
        print("\nTesting straight read and write")

    r = AmberRst7(crd)

    a = AmberPrm(prm)

    a.writeToFile("test-%s.prm" % root)

    a2 = AmberPrm("test-%s.prm" % root)

    if verbose:
        print(a, a2)

    assert_equal(a.nAtoms(), a2.nAtoms())

    if verbose:
        print("\nTesting triangle conversion to grotop")

    s = _get_first_molecule(a.toSystem(r))
    r = AmberRst7(s)

    g = GroTop(s)
    g.writeToFile("test-%s.grotop" % root)
    print(g)

    g = GroTop("test-%s.grotop" % root)
    s2 = _get_first_molecule(g.toSystem(r))

    nrgs = _get_energies(s)
    nrgs2 = _get_energies(s2)

    if verbose:
        _print_energies(nrgs, nrgs2)

    _compare_energies(nrgs, nrgs2)

    if verbose:
        print("\nComparing backwards/forwards conversion")

    a2 = AmberPrm(s)
    a2.writeToFile("test-%s.prm7" % root)
    a2 = AmberPrm("test-%s.prm7" % root)

    s2 = _get_first_molecule(a2.toSystem(r))

    nrgs2 = _get_energies(s2)

    if verbose:
        _print_energies(nrgs, nrgs2)

    _compare_energies(nrgs, nrgs2)
    


def test_cathepsin(verbose=False):
    for prm in inputs.keys():
        _test_input(prm, inputs[prm], verbose)

if __name__ == "__main__":
    test_cathepsin(True)

