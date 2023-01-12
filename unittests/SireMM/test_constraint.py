
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.MM import *
from Sire.System import *
from Sire.Base import *
from Sire.CAS import *

from nose.tools import assert_almost_equal

def test_constraint(verbose = False):

    s = System()

    lam = Symbol("lambda")

    s.setConstant(lam, 0.0)

    ff1 = InterFF("shift")
    ff1.setCLJFunction(CLJSoftShiftFunction())
    ff1.setCLJFunction("f", CLJSoftShiftFunction())
    s.add(ff1)

    ff2 = InterFF("rf")
    ff2.setCLJFunction(CLJSoftRFFunction())
    ff2.setCLJFunction("b", CLJSoftRFFunction())
    s.add(ff2)

    s.add(PropertyConstraint("alpha", ff1.name(), 1 - lam))
    s.add(PropertyConstraint("alpha[f]", ff1.name(), 0.5*lam))

    s.add(PropertyConstraint("alpha", ff2.name(), lam))
    s.add(PropertyConstraint("alpha[b]", ff2.name(), 0.25*lam))

    for i in range(0,11):
        l = 0.1 * i
        if verbose:
            print("Setting lambda equal to %s" % l)

        s.setConstant(lam, l)

        f1 = s.forceField(ff1.name())
        f2 = s.forceField(ff2.name())

        assert_almost_equal( f1.property("alpha").value(), 1 - l )
        assert_almost_equal( f1.property("alpha[f]").value(), 0.5*l )
        assert_almost_equal( f2.property("alpha").value(), l )
        assert_almost_equal( f2.property("alpha[b]").value(), 0.25*l )
    


if __name__ == "__main__":
    test_constraint(True)

