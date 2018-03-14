
from Sire.IO import *
from Sire.Mol import *
import Sire.Config

import subprocess
import shlex
import sys
import os

from nose.tools import assert_equal

sire_python = "%s/python" % Sire.Config.binary_directory

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

def test_broken_pdb(verbose=False):
    cmd = "%s %s --test_broken_pdb" % (sire_python, sys.argv[0])

    if verbose:
        print(cmd)

    env = dict(os.environ, LC_ALL="fr_FR")

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

funcs = {}
funcs["--test_broken_pdb"] = _test_broken_pdb

if __name__ == "__main__":
    if len(sys.argv) > 1:
        funcs[sys.argv[1]]()
    else:
        test_broken_pdb(True)
