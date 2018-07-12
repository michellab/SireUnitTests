
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
    p.writeToFile("test.pdb")

if __name__ == "__main__":
    test_mapping(True)

