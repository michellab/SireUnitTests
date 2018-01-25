
from Sire.IO import *
from Sire.System import *
from Sire.MM import *

from nose.tools import assert_equal, assert_almost_equal

def _addForceFields(s):
    cljff = InterFF("cljff")
    cljff.add(s.molecules())
    intraclj = IntraFF("intraclj")
    intraclj.add(s.molecules())
    intraff = InternalFF("intraff")
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
    s = MoleculeParser.read("../io/urea.top", "../io/urea.gro",
                            {"GROMACS_PATH":"../io/gromacs"})

    if verbose:
        print("Saving to amber files...")

    # save to amber, and then reload
    a = AmberPrm(s)
    c = AmberRst(s)

    a.writeToFile("test.top")
    c.writeToFile("test.rst")

    if verbose:
        print("Reading from amber files...")

    s2 = MoleculeParser.read("test.top", "test.rst")

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
        assert_almost_equal( nrgs[key], nrgs2[key], 2 )

if __name__ == "__main__":
    test_groamb(True)

