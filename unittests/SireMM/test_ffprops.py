
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.MM import *
from Sire.Base import *

def _pvt_test_props(ff, verbose = False):

    try:
        ff.setCLJFunction( CLJShiftFunction() )
        ff.setCLJFunction("f", CLJRFFunction())
        ff.setCLJFunction("next", CLJSoftShiftFunction())
    except:
        ff.setCLJFunction( CLJIntraShiftFunction() )
        ff.setCLJFunction("f", CLJIntraRFFunction() )
        ff.setCLJFunction("next", CLJSoftIntraShiftFunction())

    if verbose:
        print("Testing %s" % ff)

    ff.setProperty("dielectric[f]", wrap(2))

    d = ff.property("dielectric[f]")

    assert( d.value() == 2 )

    ff.setProperty("dielectric[all]", wrap(3))

    d = ff.property("dielectric[f]")

    assert( d.value() == 3 )

    ff.setProperty("combiningRules[all]", wrap("geometric"))
    assert( ff.property("combiningRules").value() == "geometric" )
    assert( ff.property("combiningRules[default]").value() == "geometric" )
    assert( ff.property("combiningRules[f]").value() == "geometric" )

    ff.setProperty("combiningRules", wrap("arithmetic"))
    assert( ff.property("combiningRules").value() == "arithmetic" )
    assert( ff.property("combiningRules[default]").value() == "arithmetic" )
    assert( ff.property("combiningRules[f]").value() == "geometric" )

    ff.setProperty("combiningRules[default]", wrap("geometric"))
    assert( ff.property("combiningRules").value() == "geometric" )
    assert( ff.property("combiningRules[default]").value() == "geometric" )
    assert( ff.property("combiningRules[f]").value() == "geometric" )

    ff.setProperty("combiningRules[all]", wrap("arithmetic"))
    assert( ff.property("combiningRules").value() == "arithmetic" )
    assert( ff.property("combiningRules[default]").value() == "arithmetic" )
    assert( ff.property("combiningRules[f]").value() == "arithmetic" )

    ff.setProperty("alpha[all]", wrap(5.0))
    assert( ff.property("alpha[next]").value() == 5.0 )

    ff.setProperty("alpha[next]", wrap(2.0))
    assert( ff.property("alpha[next]").value() == 2.0 )

def test_props(verbose=False):
    ff = InterFF()

    _pvt_test_props(ff, verbose)

    ff = InterGroupFF()

    _pvt_test_props(ff, verbose)

    ff = IntraFF()

    _pvt_test_props(ff, verbose)

    ff = IntraGroupFF()

    _pvt_test_props(ff, verbose)

if __name__ == "__main__":
    test_props(True)

