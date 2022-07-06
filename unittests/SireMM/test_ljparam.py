
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.MM import *
from Sire.Units import *
from Sire.Maths import *

try:
    from sire.legacy.MM import LJParameter
except Exception:
    pass

from nose.tools import assert_almost_equal

import math

def test_params(verbose=False):

    sigma = 3.451 * angstrom
    epsilon = 0.105 * kcal_per_mol

    lj = LJParameter(sigma, epsilon)

    assert_almost_equal( lj.sigma().value(), sigma.value(), 3 )
    assert_almost_equal( lj.epsilon().value(), epsilon.value(), 3 )

    A = lj.A()
    B = lj.B()

    if verbose:
        print("LJ sig/eps %s : %s equals A/B %s : %s" % \
             (sigma.value(),epsilon.value(),A,B))

    lj2 = LJParameter.fromAAndB(A, B)
 
    sigma2 = lj2.sigma()
    epsilon2 = lj2.epsilon()

    if verbose:
        print("%s versus %s" % (lj, lj2))

    assert_almost_equal( sigma.value(), sigma2.value(), 3 )
    assert_almost_equal( epsilon.value(), epsilon2.value(), 3 )

def _pvt_test_equal(lj, rule):
    com = lj.combine(lj, rule)

    assert_almost_equal( lj.sigma().value(), com.sigma().value(), 3 )
    assert_almost_equal( lj.epsilon().value(), com.epsilon().value(), 3 )

def test_combining(verbose=False):

    # ensure that LJParameter has the combine function
    try:
        a = LJParameter.combine
    except:
        if verbose:
            print("Skipping test as no support for LJParameter.combine")
        return

    rangen = RanGenerator()

    for i in range(0,100):
        if verbose:
            print("Test %d of 100..." % (i+1))

        lj0 = LJParameter( rangen.rand(0.0, 5.0)*angstrom,
                           rangen.rand(0.0,1.5)*kcal_per_mol )

        lj1 = LJParameter( rangen.rand(0.0, 5.0)*angstrom,
                           rangen.rand(0.0,1.5)*kcal_per_mol )    

        _pvt_test_equal(lj0, LJParameter.ARITHMETIC)
        _pvt_test_equal(lj0, LJParameter.GEOMETRIC)
        _pvt_test_equal(lj1, LJParameter.ARITHMETIC)
        _pvt_test_equal(lj1, LJParameter.GEOMETRIC)

        ari0 = lj0.combine(lj1, LJParameter.ARITHMETIC)
        geo0 = lj0.combine(lj1, LJParameter.GEOMETRIC)

        ari1 = lj1.combine(lj0, LJParameter.ARITHMETIC)
        geo1 = lj1.combine(lj0, LJParameter.GEOMETRIC)

        assert_almost_equal( ari0.sigma().value(), ari1.sigma().value(), 3 )
        assert_almost_equal( geo0.sigma().value(), geo1.sigma().value(), 3 )
        assert_almost_equal( ari0.epsilon().value(), ari1.epsilon().value(), 3 )
        assert_almost_equal( geo0.epsilon().value(), geo1.epsilon().value(), 3 )

        # compare against the manual calculation
        sig_ari = 0.5 * (lj0.sigma().value() + lj1.sigma().value())
        sig_geo = math.sqrt( lj0.sigma().value() * lj1.sigma().value() )

        eps = math.sqrt( lj0.epsilon().value() * lj1.epsilon().value() )

        assert_almost_equal( ari0.sigma().value(), sig_ari, 3 )
        assert_almost_equal( ari0.epsilon().value(), eps, 3 )

        assert_almost_equal( geo0.sigma().value(), sig_geo, 3 )
        assert_almost_equal( geo0.epsilon().value(), eps, 3 )

if __name__ == "__main__":
    test_params(True)
    test_combining(True)

