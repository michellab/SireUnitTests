
from Sire.IO import * 

import glob

from nose.tools import assert_equal

d = "../io/cathepsin"

prms = glob.glob("%s/*.parm7" % d)

inputs = {}

for prm in prms:
    root = prm[0:-6]
    rst = glob.glob("%s*.rst7" % root)[0]

    inputs[prm] = rst

def _test_input(prm, crd, verbose=False):
    root = prm.split("/")[-1][0:-6]

    if verbose:
        print("Testing %s (%s | %s)" % (root, prm, crd))
        print("\nTesting straight read and write")

    a = AmberPrm(prm)

    a.writeToFile("test-%s.prm" % root)

    a2 = AmberPrm("test-%s.prm" % root)

    if verbose:
        print(a, a2)

    assert_equal(a.nAtoms(), a2.nAtoms())

    if verbose:
        print("\nTesting to and from system")

    s = a.toSystem()    
    s2 = a2.toSystem()

    if verbose:
        print(s, s2)

    a2 = AmberPrm(s)

    if verbose:
        print(a, a2)

    assert_equal(a.nAtoms(), a2.nAtoms())

    s2 = a2.toSystem()

    if verbose:
        print(s, s2)

    if verbose:
        print("\nTesting from MoleculeParser")

    s = MoleculeParser.read(prm, crd)

    a2 = AmberPrm(s)

    if verbose:
        print(a, a2)

    assert_equal(a.nAtoms(), a2.nAtoms())

    if verbose:
        print("\nTesting to and from a file")

    a.writeToFile("test-%s-2.prm7" % root)

    a2 = AmberPrm("test-%s-2.prm7" % root)

    assert_equal(a.nAtoms(), a2.nAtoms())

    s2 = MoleculeParser.read("test-%s-2.prm7" % root, crd)

    if verbose:
        print(s, s2)

    


def test_cathepsin(verbose=False):
    for prm in inputs.keys():
        _test_input(prm, inputs[prm], verbose)

if __name__ == "__main__":
    test_cathepsin(True)

