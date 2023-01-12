try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Base import *

from nose.tools import assert_equal


def test_broken_prm(verbose=False):
    if verbose:
        print("Reading the file...")

    a = AmberPrm("../io/hsp90.top")

    if verbose:
        print("Creating the system in serial...")

    s1 = a.toSystem({"parallel": wrap(False)})

    if verbose:
        print("Creating the system in parallel...")

    s2 = a.toSystem()

    if verbose:
        print("Creating the system with coordinates...")

    s = MoleculeParser.read("../io/hsp90.top", "../io/hsp90.crd")

    if verbose:
        print(s)

    p = PDB2(s)
    p.writeToFile("test.pdb")


def test_broken_zan(verbose=False):

    if verbose:
        print("Reading the file...")

    a = AmberPrm("../io/ZWT.top")

    if verbose:
        print("Creating the system...")

    s = a.toSystem()


if __name__ == "__main__":
    test_broken_prm(True)
    test_broken_zan(True)
