
from Sire.IO import * 
from Sire.Mol import *
from Sire.System import *
from Sire.MM import *
from Sire.CAS import *

from nose.tools import assert_equal

def _get_first_molecules(s):
    m = MoleculeGroup("all")
    intraff = InternalFF("intraff")
    intraclj = IntraCLJFF("intraclj")

    for i in range(0,s.nMolecules()):
        mol = s[MolIdx(i)]

        # only add non-water molecules
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

prm = "../io/ssbonds/solvated-cyx-ssbond.parm7"
rst = "../io/ssbonds/solvated-cyx-ssbond.rst7"

noss_prm = "../io/ssbonds/solvated-cyx-NOssbond.parm7"
noss_rst = "../io/ssbonds/solvated-cyx-NOssbond.rst7"

def test_read_write(verbose=False):
    if verbose:
        print("Testing read/write for %s/%s" % (prm, rst))

    s = MoleculeParser.read(prm, rst)

    if verbose:
        print("Testing writing...")

    a = AmberPrm(s)
    r = AmberRst(s)

    a.writeToFile("test_ssbonds.prm7")
    r.writeToFile("test_ssbonds.rst")

    s2 = MoleculeParser.read("test_ssbonds.prm7", "test_ssbonds.rst")

    if verbose:
        print("Comparing...")
        print(s, s2)

    assert_equal(s.nMolecules(), s2.nMolecules())

    s = _get_first_molecules(s)
    s2 = _get_first_molecules(s2)

    if verbose:
        print("Calculating energies...")

    nrgs = _get_energies(s)
    nrgs2 = _get_energies(s2)

    if verbose:
        _print_energies(nrgs, nrgs2)

    _compare_energies(nrgs, nrgs2)    
    

def test_ssbonds(verbose=False):
    if verbose:
        print("Reading %s" % prm)

    a = AmberPrm(prm)

    if verbose:
        print("Converting to system")

    s = a.toSystem()

if __name__ == "__main__":
    test_ssbonds(True)
    test_read_write(True)


