try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *

from nose.tools import assert_equal


def test_props(verbose=False):

    if verbose:
        print("Loading molecules...")

    # load the molecules
    s1 = MoleculeParser.read(
        "../io/ala.top",
        "../io/ala.crd",
        {"coordinates": "fish", "charge": "puppy"},
    )

    s2 = MoleculeParser.read("../io/ala.top", "../io/ala.crd")

    if verbose:
        # print the coordinates and charge of atoms of the first
        # molecule - at the "fish"/"puppy" properties
        for atom in s1[MolIdx(0)].atoms():
            print(
                "%s : coords = %s, charge = %s"
                % (atom, atom.property("fish"), atom.property("puppy"))
            )

    for i in range(0, s1.nMolecules()):
        atoms1 = s1[MolIdx(i)].atoms()
        atoms2 = s2[MolIdx(i)].atoms()

        for j in range(0, atoms1.count()):
            assert_equal(
                atoms1[j].property("fish"), atoms2[j].property("coordinates")
            )

            assert_equal(
                atoms1[j].property("puppy"), atoms2[j].property("charge")
            )


if __name__ == "__main__":
    test_props(True)
