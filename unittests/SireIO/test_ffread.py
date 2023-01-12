try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *

from nose.tools import assert_equal


def test_ffread(verbose=False):
    try:
        gaff = getForceField("amber::gaff")
    except:
        return

    a = AmberPrm("../io/ala.top", {"forcefield": gaff})

    if verbose:
        print(a.forcefield())

    assert_equal(a.forcefield(), gaff)

    try:
        # this should raise an exception
        a = AmberPrm("../io/ala.top", {"forcefield": Atom()})
        assert_equal(True, False)
    except:
        pass

    a = AmberPrm("../io/ala.top", {"forcefield": "test"})

    assert_equal(a.forcefield(), getForceField("amber::ff"))

    if verbose:
        print("Testing that the molecules load with the right forcefield...")

    s = MoleculeParser.read(
        "../io/ala.top", "../io/ala.crd", {"forcefield": gaff}
    )

    for i in range(0, s.nMolecules()):
        ff = s[MolIdx(i)].property("forcefield")

        if verbose:
            print(s[MolIdx(i)])
            print(ff)

        assert_equal(ff, gaff)


if __name__ == "__main__":
    test_ffread(True)
