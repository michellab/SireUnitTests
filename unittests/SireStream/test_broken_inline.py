try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

import Sire.Stream

from Sire.Mol import *
from Sire.Base import *
from Sire.Maths import *
from Sire.Units import *

## This test is needed as this causes a crash when using Sire
## on Linux or Windows that is compiled with GCC >= 5 using -O3
## (-O2 seems to be ok)


def test_should_crash(verbose=False):
    mol = Molecule("test")

    mol = (
        mol.edit()
        .add(CGName("1"))
        .add(AtomName("CA"))
        .cutGroup()
        .add(AtomName("CB"))
        .cutGroup()
        .add(AtomName("CC"))
        .cutGroup()
        .add(AtomName("CD"))
        .cutGroup()
        .add(AtomName("CE"))
        .cutGroup()
        .add(AtomName("CF"))
        .molecule()
        .commit()
    )

    mol = (
        mol.atom(AtomIdx(0))
        .edit()
        .setProperty("test", 5.0 * mod_electron)
        .molecule()
        .commit()
    )

    Sire.Stream.save(mol, "testmol.s3")

    if verbose:
        print("Simple test")

    test = mol.property("test")

    if verbose:
        print(test.__class__)

    t = test.toVector()

    if verbose:
        print(t)

    if verbose:
        print("Stream test")

    mol = Sire.Stream.load("testmol.s3")
    test = mol.property("test")

    if verbose:
        print(test.__class__)

    t = test.toVector()

    if verbose:
        print(t)

    mol = Sire.Stream.load("../io/saved_molecule.s3").molecule()

    mol = (
        mol.atom(AtomIdx(0))
        .edit()
        .setProperty("test", 5.0)
        .setProperty("chg3", 5 * mod_electron)
        .molecule()
        .commit()
    )

    if verbose:
        print("test property")

    test = mol.property("test")

    if verbose:
        print(test.__class__)
        print("chg property")

    chg = mol.property("charge")

    if verbose:
        print(chg.__class__)
        print("chg2 property")

    chg2 = Sire.Stream.load("../io/saved_charges.s3")

    if verbose:
        print(chg2.__class__)
        print("chg3 property")

    chg3 = mol.property("chg3")

    if verbose:
        print(chg3.__class__)
        print("lj property")

    lj = mol.property("LJ")

    if verbose:
        print(lj.__class__)

    coords = mol.property("coordinates")

    if verbose:
        print("test.toVector()")

    t = test.toVector()

    if verbose:
        print(t)

    if verbose:
        print("chg(all).toVector()")

    t3 = chg3.toVector()
    t2 = chg2.toVector()
    t1 = chg.toVector()

    if verbose:
        print(t1)
        print(t2)
        print(t3)
        print("lj.toVector()")

    t = lj.toVector()

    if verbose:
        print(t)
        print(lj)

        print("coords.toVector()")

    t = coords.toVector()

    if verbose:
        print(t)
        print(coords)

        print("\nTHIS WILL CRASH IF BROKEN COMPILATION")

    s = str(test)

    if verbose:
        print(s)
        print("\nTHESE WILL CRASH IF BROKEN COMPILATION")

    s1 = str(chg3)
    s2 = str(chg)

    if verbose:
        print(s1)
        print(s2)


if __name__ == "__main__":
    test_should_crash(True)
