from Sire.IO import *
from Sire.Mol import *
from Sire.System import *

from nose.tools import assert_equal

def test_mapping(verbose=False):
    if verbose:
        print("Loading the molecules...")

    ose = MoleculeParser.read("../io/ose.crd", "../io/ose.top") \
                        .search("resname /ose/i")[0]

    zan = MoleculeParser.read("../io/zan.rst", "../io/zan.top") \
                        .search("resname /zan/i")[0]

    if verbose:
        print("Performing the mapping...")

    # force H1 in OSE to equal H18 in zan
    user_match = { AtomName("H1") : AtomName("H18"),
                   AtomName("C9") : AtomName("O4") }

    mapping = ose.evaluate().findMCS(zan, AtomIDMatcher(user_match), True)

    if verbose:
        keys = list(mapping.keys())

        for key in keys:
            atom0 = ose.atom(key)
            atom1 = zan.atom(mapping[key])

            print("%s => %s" % (atom0,atom1))

    zan = zan.move().align(ose, AtomResultMatcher(mapping,True))

    m = MoleculeGroup("all")
    m.add(ose)
    m.add(zan)

    s = System()
    s.add(m)

    p = PDB2(s)
    p.writeToFile("test_mapping.pdb")

def test_planar(verbose=False):
    """A test for mapping with near planar molecule."""

    if verbose:
        print("Loading the molecules...")

    s0 = MoleculeParser.read("../io/CHEMBL151480.mol2")
    s1 = MoleculeParser.read("../io/CHEMBL95097.mol2")

    m0 = s0.molecule(MolIdx(0))
    m1 = s1.molecule(MolIdx(0))

    if verbose:
        print("Performing the mapping...")

    mapping = m0.evaluate().findMCS(m1, AtomIDMatcher(), True)

    if verbose:
        keys = list(mapping.keys())

        for key in keys:
            atom0 = m0.atom(key)
            atom1 = m1.atom(mapping[key])

            print("%s => %s" % (atom0,atom1))

    m1 = m1.move().align(m0, AtomResultMatcher(mapping,True))

    m = MoleculeGroup("all")
    m.add(m0)
    m.add(m1)

    s = System()
    s.add(m)

    p = PDB2(s)
    p.writeToFile("test_planar_mapping.pdb")

if __name__ == "__main__":
    #test_mapping(True)
    test_planar(True)

