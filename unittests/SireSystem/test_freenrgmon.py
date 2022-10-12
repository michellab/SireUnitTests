from sire.legacy.System import *
from sire.legacy.IO import *
from sire.legacy.Mol import *
from sire.legacy.MM import *
from sire.legacy.CAS import *
from sire.legacy.Maths import *
from sire.legacy.Base import *

from nose.tools import assert_almost_equal

s = MoleculeParser.read("../io/waterbox.crd", "../io/waterbox.top")

mols = s.molecules()

molnums = mols.molNums()
molnums.sort()

water0 = mols[molnums[10]].molecule()
water1 = mols[molnums[11]].molecule()
water2 = mols[molnums[12]].molecule()

# translate water0 to the center
water0 = water0.move().translate( -water0.evaluate().center() ).commit()

# translate water1 and water2 to (0,0,2)
water1 = water1.move().translate( Vector(0,0,2) - water1.evaluate().center() ).commit()
water2 = water2.move().translate( Vector(0,0,2) - water2.evaluate().center() ).commit()

def _pvt_calculateEnergy(lamval, verbose):
    cljff01 = InterGroupCLJFF("cljff01")
    cljff02 = InterGroupCLJFF("cljff02")

    cljff01.add( water0, MGIdx(0) )
    cljff01.add( water1, MGIdx(1) )

    cljff02.add( water0, MGIdx(0) )
    cljff02.add( water2, MGIdx(1) )

    soft_cljff01 = InterGroupSoftCLJFF("soft_cljff01")
    soft_cljff02 = InterGroupSoftCLJFF("soft_cljff02")

    soft_cljff01.add( water0, MGIdx(0) )
    soft_cljff01.add( water1, MGIdx(1) )

    soft_cljff02.add( water0, MGIdx(0) )
    soft_cljff02.add( water2, MGIdx(1) )

    ref = MoleculeGroup("ref", water0)
    group_a = MoleculeGroup("group_a", water1)
    group_b = MoleculeGroup("group_b", water2)

    dlam = 0.001
    lamval_f = lamval + dlam

    lam = Symbol("lambda")
    lam_f = Symbol("lambda_f")

    soft_nrg = (1-lam) * soft_cljff01.components().total(0) + lam * soft_cljff02.components().total(0)
    soft_nrg_f = (1-lam_f) * soft_cljff01.components().total(1) + lam_f * soft_cljff02.components().total(1)

    de_soft = soft_nrg_f - soft_nrg

    nrg = ((1-lam) * cljff01.components().total()) + (lam * cljff02.components().total())
    nrg_f = ((1-lam_f) * cljff01.components().total()) + (lam_f * cljff02.components().total())

    de = nrg_f - nrg

    soft_cljff01.setProperty("alpha0", wrap(lamval))
    soft_cljff02.setProperty("alpha0", wrap(1-lamval))
    soft_cljff01.setProperty("alpha1", wrap(lamval_f))
    soft_cljff02.setProperty("alpha1", wrap(1-lamval_f))

    soft_cljff01.setProperty("coulombPower", wrap(0))
    soft_cljff01.setProperty("shiftDelta", wrap(1.1))
    soft_cljff02.setProperty("coulombPower", wrap(0))
    soft_cljff02.setProperty("shiftDelta", wrap(1.1))

    sys = System()
    sys.add(cljff01)
    sys.add(cljff02)
    sys.add(soft_cljff01)
    sys.add(soft_cljff02)
    sys.add(ref)
    sys.add(group_a)
    sys.add(group_b)
    sys.setComponent(lam, lamval)
    sys.setComponent(lam_f, lamval_f)
    sys.setComponent(sys.totalComponent(), nrg)
    sys.setComponent(Symbol("E_{total_f}"), nrg_f)
    sys.setComponent(Symbol("dE"), de)
    sys.setComponent(Symbol("E_soft_{total}"), soft_nrg)
    sys.setComponent(Symbol("E_soft_{total_f}"), soft_nrg_f)
    sys.setComponent(Symbol("dE_soft"), de_soft)

    nrgmon = FreeEnergyMonitor(ref, group_a, group_b)

    soft_nrgmon = FreeEnergyMonitor(ref, group_a, group_b)
    soft_nrgmon.setCoulombPower(0)
    soft_nrgmon.setShiftDelta(1.1)

    sys.add( "nrgmon", nrgmon )
    sys.add( "soft_nrgmon", soft_nrgmon)

    sys.collectStats()

    nrgmon = sys[ MonitorName("nrgmon") ]
    dg = nrgmon.freeEnergies()[0].average()

    soft_nrgmon = sys[ MonitorName("soft_nrgmon") ]
    soft_dg = soft_nrgmon.freeEnergies()[0].average()

    sys_dg = sys.energy(Symbol("dE")).value()
    sys_soft_dg = sys.energy(Symbol("dE_soft")).value()

    if verbose:
        print("%s : %s versus %s (should be equal)" % (lamval,dg,sys_dg))
        print("%s : %s versus %s (should be equal)" % (lamval,soft_dg,sys_soft_dg))

    assert_almost_equal(dg, sys_dg, 3)
    assert_almost_equal(soft_dg, sys_soft_dg, 3)

def test_1(verbose=False):
    _pvt_calculateEnergy(0.0,verbose)

def test_2(verbose=False):
    _pvt_calculateEnergy(0.1,verbose)

def test_3(verbose=False):
    _pvt_calculateEnergy(0.5,verbose)

def test_4(verbose=False):
    _pvt_calculateEnergy(0.9,verbose)

def test_5(verbose=False):
    _pvt_calculateEnergy(0.999,verbose)

if __name__ == "__main__":
    test_1(True)
    test_2(True)
    test_3(True)
    test_4(True)
    test_5(True)

