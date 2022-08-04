
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass
  
from Sire.MM import *
from Sire.FF import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Vol import *
from Sire.Base import *
from Sire.Units import *
from Sire.Qt import *

from nose.tools import assert_almost_equal, assert_equal

(mols, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

slow_ff1 = InterCLJFF("ff1")
slow_ff2 = InterCoulombFF("ff2")
slow_ff3 = InterLJFF("ff3")

slow_ff1.add(mols)
slow_ff2.add(mols)
slow_ff3.add(mols)

fast_ff1 = InterFF("ff1")
fast_ff1.setCLJFunction( CLJShiftFunction(100*angstrom, 100*angstrom) )
fast_ff1.setUseParallelCalculation(False)

fast_ff2 = InterFF("ff2")
fast_ff2.setCLJFunction( CLJShiftFunction(100*angstrom, 0*angstrom) )
fast_ff2.setUseParallelCalculation(False)

fast_ff3 = InterFF("ff3")
fast_ff3.setCLJFunction( CLJShiftFunction(0*angstrom, 100*angstrom) )
fast_ff3.setUseParallelCalculation(False)

fast_ff1.add(mols)
fast_ff2.add(mols)
fast_ff3.add(mols)

def _pvt_test(ff1, ff2, ff3, par=False, verbose=False):

    ffields = ForceFields()
    assert_equal( ffields.nForceFields(), 0 )
    assert_equal( ffields.nMolecules(), 0 )

    if par:
        ff1.setUseParallelCalculation(True)
        ff2.setUseParallelCalculation(True)
        ff3.setUseParallelCalculation(True)

    ffields.add(ff1)
    assert_equal( ffields.nForceFields(), 1 )
    assert_equal( ffields.nMolecules(), ff1.nMolecules() )

    ffields.add(ff2)
    assert_equal( ffields.nForceFields(), 2 )
    assert_equal( ffields.nMolecules(), ff1.nMolecules() )

    ffields.add(ff3)
    assert_equal( ffields.nForceFields(), 3 )
    assert_equal( ffields.nMolecules(), ff1.nMolecules() )

    ffields.mustNowRecalculateFromScratch()
    assert( ffields.isDirty() )

    if verbose:
        print("Calculating total energy...")

    nrg = ffields.energy()

    if verbose:
        print("...equals %s" % nrg)
        print("Calculating individual energies...")

    t = QElapsedTimer()
    t.start()
    nrg1 = ff1.energy()
    nrg2 = ff2.energy()
    nrg3 = ff3.energy()
    ns = t.nsecsElapsed()

    if verbose:
        print("Serial evaluation took %f ms" % (0.000001*ns))

    total_nrg = nrg1+nrg2+nrg3

    if verbose:
        print("...equals %s" % total_nrg)

    assert_almost_equal( nrg.value(), total_nrg.value(), 2 )

    ffields.setComponent( ffields.totalComponent(),
                          ff1.components().total() + ff2.components().total() +
                          ff3.components().total() )

    ffields.mustNowRecalculateFromScratch()
    assert(ffields.isDirty())

    if verbose:
        print("Calculating total energy...")

    nrg = ffields.energy()

    assert_almost_equal( nrg.value(), total_nrg.value(), 2 )

    lam = Symbol("lambda")

    ffields.setComponent( ffields.totalComponent(),
                          lam * (ff1.components().total() + ff2.components().total()) +
                          ff3.components().total() )

    for i in range(0,11):
        l = 0.1*i

        ffields.mustNowRecalculateFromScratch()
        assert(ffields.isDirty())

        if verbose:
            print("Calculating energy for lambda = %f" % l)

        ffields.setComponent(lam, l)

        nrg = ffields.energy()
        total_nrg = l * (nrg1 + nrg2) + nrg3

        assert_almost_equal( nrg.value(), total_nrg.value(), 2 )

    part = Symbol("E_part")

    ffields.setComponent( part, lam * (ff1.components().total() + ff2.components().total()) )
    ffields.setComponent( ffields.totalComponent(),
                          part + ff3.components().total() )

    for i in range(0,11):
        l = 0.1*i
    
        ffields.mustNowRecalculateFromScratch()
        assert(ffields.isDirty())

        if verbose:
            print("Calculating energy for lambda = %f" % l)

        ffields.setComponent(lam, l)

        t.start()
        nrg = ffields.energy()
        ns = t.nsecsElapsed()
        total_nrg = l * (nrg1 + nrg2) + nrg3

        if verbose:
            print("Evaluation took %f ms" % (0.000001*ns))

        assert_almost_equal( nrg.value(), total_nrg.value(), 2 )    

    if par:
        ff1.setUseParallelCalculation(False)
        ff2.setUseParallelCalculation(False)                                          
        ff3.setUseParallelCalculation(False)

def test_slow_ffs(verbose=False):
    if verbose:
        print("\nTesting parallel running of (slow) serial forcefields")

    _pvt_test(slow_ff1,slow_ff2,slow_ff3,False,verbose)

def test_parallel_ffs(verbose=False):
    if verbose:
        print("\nTesting parallel running of parallel forcefields")

    _pvt_test(fast_ff1,fast_ff2,fast_ff3,True,verbose)

def test_fast_ffs(verbose=False):
    if verbose:
        print("\nTesting parallel running of fast (serial) forcefields")

    _pvt_test(fast_ff1,fast_ff2,fast_ff3,False,verbose)

if __name__ == "__main__":
    test_fast_ffs(True)
    test_parallel_ffs(True)
    test_slow_ffs(True)

