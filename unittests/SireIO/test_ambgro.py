try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.System import *
from Sire.MM import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.Units import *

from nose.tools import assert_equal, assert_almost_equal


def _assert_almost_equal(a, b, tol=0.5):
    if abs(a - b) > tol:
        assert_almost_equal(a, b, 2)


def _addForceFields(s):
    cljff = InterFF("cljff")
    cljff.add(s.molecules())
    intraclj = IntraFF("intraclj")
    intraff = InternalFF("intraff")
    intraff.enable14Calculation()

    # only calculate intramolecular energy of non-solvent molecules
    for molnum in s.molNums():
        mol = s[molnum]
        if mol.nAtoms() > 3:
            intraclj.add(mol)
            intraff.add(mol)

    s = System(s)
    s.add(cljff)
    s.add(intraff)
    s.add(intraclj)

    return s


def test_small(verbose=False):
    if verbose:
        print("Reading small amber files...")

    files = [
        ("../io/ose.top", "../io/ose.crd"),
        ("../io/ala.top", "../io/ala.crd"),
    ]

    for file in files:
        if verbose:
            print("\nReading %s/%s" % file)

        s = MoleculeParser.read(file[0], file[1])

        if verbose:
            print("Converting to gromacs...")

        gtop = GroTop(s)
        g87 = Gro87(s)

        gtop.writeToFile("test.grotop")
        g87.writeToFile("test.gro")

        if verbose:
            print("Re-reading from gromacs...")

        s2 = MoleculeParser.read("test.gro", "test.grotop")

        if verbose:
            print("Adding forcefields...")

        s = _addForceFields(s)
        s2 = _addForceFields(s2)

        if verbose:
            print("Calculating energies...")

        nrgs = s.energies()
        # merge the dihedral and improper energies...
        dihsym = Symbol("E_{intraff}^{dihedral}")
        impsym = Symbol("E_{intraff}^{improper}")
        nrgs.set(dihsym, nrgs[dihsym] + nrgs[impsym])
        nrgs.set(impsym, 0.0)

        nrgs2 = s2.energies()

        keys = list(nrgs.keys())
        keys.sort()

        if verbose:
            print("Energies:")
            for key in keys:
                print("%s:  %s  %s" % (key, nrgs[key], nrgs2[key]))

        for key in keys:
            _assert_almost_equal(nrgs[key], nrgs2[key])


def test_ambgro(verbose=False):
    if verbose:
        print("Reading amber files...")

    # load up the gromacs file
    s = MoleculeParser.read("../io/cchept.parm7", "../io/cchept.rst7")

    if verbose:
        print("Saving to gromacs files...")

    # save to amber, and then reload
    gtop = GroTop(s)
    g87 = Gro87(s)

    gtop.writeToFile("test.grotop")
    g87.writeToFile("test.g87")

    if verbose:
        print("Reading from gromacs files...")

    s2 = MoleculeParser.read("test.grotop", "test.g87")

    if verbose:
        print("Reading from the parmed gromacs files...")

    s3 = MoleculeParser.read("../io/cchept.grotop", "../io/cchept.gro")

    if verbose:
        print("Writing back to amber...")

    a = AmberPrm(s2)
    r = AmberRst7(s2)

    a.writeToFile("test.prm7")
    r.writeToFile("test.rst7")

    if verbose:
        print("Reading back from amber...")

    s4 = MoleculeParser.read("test.prm7", "test.rst7")

    if verbose:
        print("Writing parmed to amber...")

    a = AmberPrm(s3)
    r = AmberRst7(s3)

    a.writeToFile("test.prm7")
    r.writeToFile("test.rst7")

    if verbose:
        print("Reading parmed back from amber...")

    s5 = MoleculeParser.read("test.prm7", "test.rst7")

    if verbose:
        print("Adding forcefields...")

    s = _addForceFields(s)
    s2 = _addForceFields(s2)
    s3 = _addForceFields(s3)
    s4 = _addForceFields(s4)
    s5 = _addForceFields(s5)

    if verbose:
        print("Calculating energies...")

    nrgs = s.energies()

    # merge the dihedral and improper energies...
    dihsym = Symbol("E_{intraff}^{dihedral}")
    impsym = Symbol("E_{intraff}^{improper}")
    nrgs.set(dihsym, nrgs[dihsym] + nrgs[impsym])
    nrgs.set(impsym, 0.0)

    nrgs2 = s2.energies()
    nrgs3 = s3.energies()
    nrgs4 = s4.energies()
    nrgs5 = s5.energies()

    nrgs4.set(dihsym, nrgs4[dihsym] + nrgs4[impsym])
    nrgs4.set(impsym, 0.0)
    nrgs5.set(dihsym, nrgs5[dihsym] + nrgs5[impsym])
    nrgs5.set(impsym, 0.0)

    keys = list(nrgs.keys())
    keys.sort()

    if verbose:
        print("Energies:")
        for key in keys:
            print(
                "%s:  %s  %s  %s  %s  %s"
                % (
                    key,
                    nrgs[key],
                    nrgs2[key],
                    nrgs3[key],
                    nrgs4[key],
                    nrgs5[key],
                )
            )

    for key in keys:
        _assert_almost_equal(nrgs[key], nrgs2[key])
        _assert_almost_equal(nrgs2[key], nrgs4[key])
        _assert_almost_equal(nrgs3[key], nrgs5[key])
        # _assert_almost_equal( nrgs[key], nrgs3[key] )


if __name__ == "__main__":
    test_small(True)
    test_ambgro(True)
