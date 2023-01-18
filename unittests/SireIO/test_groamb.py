try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.System import *
from Sire.MM import *
from Sire.Mol import *

from nose.tools import assert_equal, assert_almost_equal


def _addForceFields(s):
    cljff = InterFF("cljff")
    cljff.add(s.molecules())
    intraclj = IntraFF("intraclj")
    intraclj.add(s.molecules())
    intraff = InternalFF("intraff")
    intraff.enable14Calculation()
    intraff.add(s.molecules())

    s = System(s)
    s.add(cljff)
    s.add(intraff)
    s.add(intraclj)

    return s


def test_groamb(verbose=False):
    if verbose:
        print("Reading gromacs file...")

    # load up the gromacs file
    s = MoleculeParser.read(
        "../io/urea.top", "../io/urea.gro", {"GROMACS_PATH": "../io/gromacs"}
    )

    # assert that all of the molecules have an amber-style forcefield
    for i in range(0, s.nMolecules()):
        ff = s[MolIdx(i)].property("forcefield")

        if verbose:
            print("Forcefield %s ==\n%s" % (s[MolIdx(i)], ff))

        assert_equal(ff.isAmberStyle(), True)

    if verbose:
        print("Saving to amber files...")

    # save to amber, and then reload
    a = AmberPrm(s)
    c = AmberRst(s)

    a.writeToFile("test_urea.top")
    c.writeToFile("test_urea.rst")

    if verbose:
        print("Reading from amber files...")

    s2 = MoleculeParser.read("test_urea.top", "test_urea.rst")

    if verbose:
        print("Adding forcefields...")

    s = _addForceFields(s)
    s2 = _addForceFields(s2)

    if verbose:
        print("Calculating energies...")

    nrgs = s.energies()
    nrgs2 = s2.energies()

    keys = list(nrgs.keys())
    keys.sort()

    if verbose:
        print("Energies:")
        for key in keys:
            print("%s:  %s  versus  %s" % (key, nrgs[key], nrgs2[key]))

    for key in keys:
        assert_almost_equal(nrgs[key], nrgs2[key], 2)


def test_ryckaert_bellemans(verbose=False):
    if verbose:
        print("\nTesting Ryckaert-Bellemans dihedral conversion...")
        print("Reading gromacs file...")

    # load up the gromacs file
    s = MoleculeParser.read(
        "../io/cage_quin1.top",
        "../io/cage_quin1.gro",
        {"GROMACS_PATH": "../io/cage_quin1"},
    )

    if verbose:
        print(f"Read {s.nAtoms()} atoms...")
        print("Saving to amber files...")

    # save to amber, and then reload
    a = AmberPrm(s)
    c = AmberRst(s)

    a.writeToFile("test_cage.top")
    c.writeToFile("test_cage.rst")

    if verbose:
        print(f"Written {a.nAtoms()} / {c.nAtoms()} atoms")
        print("Reading from amber files...")

    s2 = MoleculeParser.read("test_cage.top", "test_cage.rst")

    if verbose:
        print("Adding forcefields...")

    s = _addForceFields(s)
    s2 = _addForceFields(s2)

    if verbose:
        print("Calculating energies...")

    nrgs = s.energies()
    nrgs2 = s2.energies()

    keys = list(nrgs.keys())
    keys.sort()

    if verbose:
        print("Energies:")
        for key in keys:
            print("%s:  %s  versus  %s" % (key, nrgs[key], nrgs2[key]))

    for key in keys:
        assert_almost_equal(nrgs[key], nrgs2[key], 2)


if __name__ == "__main__":
    test_groamb(True)
    test_ryckaert_bellemans(True)
