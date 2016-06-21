  
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

ff1 = InterCLJFF("ff1")
ff2 = InterCoulombFF("ff2")
ff3 = InterLJFF("ff3")

ff1.add(mols)
ff2.add(mols)
ff3.add(mols)

def test_ffs_total(verbose=False):

    ffields = ForceFields()
    assert_equal( ffields.nForceFields(), 0 )
    assert_equal( ffields.nMolecules(), 0 )

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

    nrg1 = ff1.energy()
    nrg2 = ff2.energy()
    nrg3 = ff3.energy()

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

if __name__ == "__main__":
    test_ffs_total(True)
