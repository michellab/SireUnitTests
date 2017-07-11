
from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.Maths import *

from nose.tools import assert_equal, assert_almost_equal

rangen = RanGenerator()

s = MoleculeParser.read("../io/proteinbox.top","../io/proteinbox.crd")

mols = s.molecules()
#mols.add( s[ MolWithResID(ResName("ALA")) ] )

def _assert_almost_equal(oldpots, newpots):

    assert_equal( len(oldpots), len(newpots) )

    has_atom2 = True
    has_atom3 = True

    try:
        oldpots[0].atom2()
    except:
        has_atom2 = False

    try:
        oldpots[0].atom3()
    except:
        has_atom3 = False

    for i in range(0,len(oldpots)):
        op = oldpots[i]
        np = newpots[i]

        if op != np:
            assert_equal(op.atom0(), np.atom0())
            assert_equal(op.atom1(), np.atom1())
            if has_atom2:
                assert_equal(op.atom2(), np.atom2())
            if has_atom3:
                assert_equal(op.atom3(), np.atom3())

            of = op.function()
            nf = np.function()

            if of != nf:
                s = of.symbols()[0]

                for i in range(0,50):
                    vals = {s:rangen.rand(-5.0,5.0)}

                    assert_almost_equal( of.evaluate(vals), nf.evaluate(vals), 5 )

def _assert_equal(oldparams, newparams):
    if oldparams == newparams:
        return

    assert_equal( oldparams.charges(), newparams.charges() )
    assert_equal( oldparams.masses(), newparams.masses() )
    assert_equal( oldparams.elements(), newparams.elements() )
    assert_equal( oldparams.ljs(), newparams.ljs() )
    assert_equal( oldparams.amberTypes(), newparams.amberTypes() )
    assert_equal( oldparams.gbRadii(), newparams.gbRadii() )
    assert_equal( oldparams.gbScreening(), newparams.gbScreening() )
    assert_equal( oldparams.treeChains(), newparams.treeChains() )
    assert_equal( oldparams.radiusSet(), newparams.radiusSet() )
    assert_equal( oldparams.excludedAtoms(), newparams.excludedAtoms() )
    assert_equal( oldparams.connectivity(), newparams.connectivity() )
    _assert_almost_equal( oldparams.bondFunctions().potentials(), 
                          newparams.bondFunctions().potentials() )
    _assert_almost_equal( oldparams.angleFunctions().potentials(), 
                          newparams.angleFunctions().potentials() )
    _assert_almost_equal( oldparams.dihedralFunctions().potentials(), 
                          newparams.dihedralFunctions().potentials() )
    _assert_almost_equal( oldparams.improperFunctions().potentials(), 
                          newparams.improperFunctions().potentials() )
    assert_equal( oldparams.nb14s(), newparams.nb14s() )
    assert_equal( oldparams.cljScaleFactors(), newparams.cljScaleFactors() )

def test_params(verbose=False):

    for molnum in mols.molNums():
        mol = mols[molnum].molecule()
        oldparams = mol.property("parameters")

        if verbose:
            print("Old parameters = %s" % oldparams)

        newparams = AmberParams(mol)

        if verbose:
            print("New parameters = %s" % newparams)

        _assert_equal(oldparams, newparams)

if __name__ == "__main__":
    test_params(True)


