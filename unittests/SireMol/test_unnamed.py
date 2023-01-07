try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.Maths import *
from Sire.Mol import *

from nose.tools import assert_equal


def test_unnamed(verbose=False):
    na = Molecule()
    na = na.edit().add(CGName("1")).add(AtomName("Na")).molecule().commit()

    na = (
        na.edit()
        .atom(AtomIdx(0))
        .renumber(AtomNum(5))
        .setProperty("coordinates", Vector(1, 2, 3))
        .molecule()
        .add(ResName("ION"))
        .renumber(ResNum(42))
        .molecule()
        .atom(AtomIdx(0))
        .reparent(ResName("ION"))
        .molecule()
        .commit()
    )

    if verbose:
        print(na)

        for atom in na.atoms():
            print(atom, atom.property("coordinates"))

        for residue in na.residues():
            print(residue)

    assert_equal(na.nAtoms(), 1)
    assert na.number().value() != 0

    assert_equal(na.atom(AtomIdx(0)).property("coordinates"), Vector(1, 2, 3))
    assert_equal(na.atom(AtomIdx(0)).name().value(), "Na")
    assert_equal(na.atom(AtomIdx(0)).number().value(), 5)
    assert_equal(na.atom(AtomIdx(0)).residue(), na.residue(ResIdx(0)))
    assert_equal(na.residue(ResIdx(0)).name().value(), "ION")
    assert_equal(na.residue(ResIdx(0)).number().value(), 42)


if __name__ == "__main__":
    test_unnamed(True)
