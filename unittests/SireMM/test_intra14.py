
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.System import *
from Sire.MM import *
from Sire.Units import *

from nose.tools import assert_almost_equal

s = MoleculeParser.read("../io/NA16.top", "../io/NA16.rst")

coul_cutoff = 7.5 * angstrom
lj_cutoff = 7.5 * angstrom

def _make_old(s):
    s2 = System(s)

    intraclj = IntraCLJFF("intraclj")
    intraclj.setShiftElectrostatics(True)
    intraclj.setSwitchingFunction( HarmonicSwitchingFunction(
                                      coul_cutoff,coul_cutoff,
                                      lj_cutoff,lj_cutoff) )
    intraclj.setSpace(s.property("space"))

    intraclj.add(s.molecules())

    intraff = InternalFF("intraff")
    intraff.add(s.molecules())

    s2.add(intraff)
    s2.add(intraclj)

    return s2

def _make_new(s):
    s2 = System(s)

    cljfunc = CLJIntraShiftFunction()
    cljfunc.setCoulombCutoff(coul_cutoff)
    cljfunc.setLJCutoff(lj_cutoff)
    cljfunc.setSpace(s.property("space"))

    intraclj = IntraFF("intraclj")
    intraclj.setCLJFunction(cljfunc)
    intraclj.add(s.molecules())

    intraff = InternalFF("intraff")
    intraff.enable14Calculation()
    intraff.add(s.molecules())

    s2.add(intraclj)
    s2.add(intraff)

    return s2

def _print_energies(s1, s2):
    nrgs1 = s1.energies()
    nrgs2 = s2.energies()

    keys1 = list(nrgs1.keys())
    keys1.sort()
    keys2 = list(nrgs2.keys())
    keys2.sort()

    for i in range(0,len(keys1)):
        print("%s : %s versus %s" % (keys1[i], nrgs1[keys1[i]], nrgs2[keys2[i]]))

def test_intra14(verbose=True):
    if verbose:
        print("Making the system with 1-4 in IntraCLJFF...")

    s1 = _make_old(s)

    if verbose:
        print("Making the system with 1-4 in InternalFF...")

    s2 = _make_new(s)

    if verbose:
        print("Calculating energies...")

    nrg1 = s1.energy()
    nrg2 = s2.energy()

    if verbose:
        print("%s versus %s" % (nrg1, nrg2))
        _print_energies(s1,s2)

    assert_almost_equal(nrg1.value(), nrg2.value(), 2)

if __name__ == "__main__":
    test_intra14(True)
