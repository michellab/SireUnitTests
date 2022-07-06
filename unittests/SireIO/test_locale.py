
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.Base import *

import Sire.Config

import subprocess
import shlex
import sys
import os

from nose.tools import assert_equal

sire_python = getBinDir() + "/sire_python"

gromacs_path = StringProperty("../io/gromacs")

def _test_broken_pdb():
    try:
        verbose = os.getenv("VERBOSE_TEST")
    except:
        verbose = False

    if verbose:
        print("Loading the original PDB...")

    s = MoleculeParser.read("../io/aladip.pdb")

    m = s[MolIdx(0)]

    if verbose:
        print("Writing the PDB using locale '%s'" % os.getenv("LC_ALL"))

    # write this PDB using the environment locale
    PDB().write(m, "test_broken_pdb.pdb")

    if verbose:
        print("Re-reading the PDB...")

    # read this back to see if we are using commas or numbers
    s = MoleculeParser.read("test_broken_pdb.pdb")

    m2 = s[MolIdx(0)]

    assert_equal( m.nAtoms(), m2.nAtoms() )

    for i in range(0, m.nAtoms()):
        if verbose:
            print("%s vs %s" % (m.atoms()[i].property("coordinates"),
                                m2.atoms()[i].property("coordinates")))

        assert_equal( m.atoms()[i].property("coordinates"),
                      m2.atoms()[i].property("coordinates") )

    # now try the same with PDB2
    if verbose:
        print("Checking PDB2. Writing to file...")

    p = PDB2(s)
    p.writeToFile("test_broken_pdb.pdb")

    if verbose:
        print("Re-reading from the file...")

    s = MoleculeParser.read("test_broken_pdb.pdb")

    m2 = s[MolIdx(0)]

    assert_equal( m.nAtoms(), m2.nAtoms() )

    for i in range(0, m.nAtoms()):
        if verbose:
            print("%s vs %s" % (m.atoms()[i].property("coordinates"),
                                m2.atoms()[i].property("coordinates")))

        assert_equal( m.atoms()[i].property("coordinates"),
                      m2.atoms()[i].property("coordinates") )

def _test_broken_rst7():
    try:
        verbose = os.getenv("VERBOSE_TEST")
    except:
        verbose = False

    if verbose:
        print("Loading the original RST7...")

    s = MoleculeParser.read("../io/ala.top", "../io/ala.crd")

    if verbose:
        print("Writing out again...")

    a = AmberRst7(s)

    a.writeToFile("test.rst7")

    if verbose:
        print("Re-reading...")

    s2 = MoleculeParser.read("../io/ala.top", "test.rst7")

    if verbose:
        print("Validating...")

    for i in range(0,s.nMolecules()):
        m1 = s[MolIdx(i)]
        m2 = s[MolIdx(i)]

        for j in range(0,m1.nAtoms()):
            v1 = m1.atoms()[j].property("coordinates")
            v2 = m2.atoms()[j].property("coordinates")

            assert_equal(v1, v2)

    if verbose:
        print("All ok :-)")

def _test_broken_gro():
    try:
        verbose = os.getenv("VERBOSE_TEST")
    except:
        verbose = False

    if verbose:
        print("Loading the original GRO...")

    s = MoleculeParser.read("../io/urea.top", "../io/urea.gro",
                            {"GROMACS_PATH":gromacs_path})

    if verbose:
        print("Writing to a new file...")

    g = Gro87(s)
    g.writeToFile("test.gro")

    if verbose:
        print("Re-reading...")

    s = MoleculeParser.read("../io/urea.top", "test.gro",
                            {"GROMACS_PATH":gromacs_path})

    if verbose:
        print("Validating...")

    for i in range(0,s.nMolecules()):
        m1 = s[MolIdx(i)]
        m2 = s[MolIdx(i)]

        for j in range(0,m1.nAtoms()):
            v1 = m1.atoms()[j].property("coordinates")
            v2 = m2.atoms()[j].property("coordinates")

            assert_equal(v1, v2)

    if verbose:
        print("All ok :-)")

def _test_broken_function(verbose, function):
    cmd = "%s %s %s" % (sire_python, "test_locale.py", function)

    if verbose:
        print(cmd)

    env = dict(os.environ, LC_ALL="it_IT.UTF-8")

    if verbose:
        env["VERBOSE_TEST"] = "1"

    p = subprocess.Popen(shlex.split(cmd),
                         env=env,
                         stdout=subprocess.PIPE)
    p.wait()

    if verbose:
        for line in p.stdout.readlines():
            print(str(line))

    if p.returncode != 0:
        assert(False)

def test_broken_pdb(verbose=False):
    _test_broken_function(verbose, "--test_broken_pdb")

def test_broken_rst7(verbose=False):
    _test_broken_function(verbose, "--test_broken_rst7")

def test_broken_gro(verbose=False):
    _test_broken_function(verbose, "--test_broken_gro")

funcs = {}
funcs["--test_broken_pdb"] = _test_broken_pdb
funcs["--test_broken_rst7"] = _test_broken_rst7
funcs["--test_broken_gro"] = _test_broken_gro

if __name__ == "__main__":
    if len(sys.argv) > 1:
        funcs[sys.argv[1]]()
    else:
        test_broken_pdb(True)
        test_broken_rst7(True)
        test_broken_gro(True)
    #print("DISABLING LOCALE TEST AS INFINITE LOOP")
