try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Maths import *

from nose.tools import assert_equal


def _test_box(rstfile, boxdims, verbose):
    a = AmberRst7(rstfile)

    if verbose:
        print("%s vs %s" % (a.boxDimensions(), boxdims))

    assert_equal(a.boxDimensions(), boxdims)


def test_box(verbose=False):
    _test_box(
        "../io/boxtest.rst7",
        Vector(84.1587143, 85.2377853, 88.7002182),
        verbose,
    )
    _test_box(
        "../io/NA16.rst7", Vector(89.1457140, 90.4230920, 94.0961620), verbose
    )
    _test_box(
        "../io/SYSTEM.crd", Vector(65.6399495, 61.9639270, 61.3930179), verbose
    )
    _test_box(
        "../io/ala.crd", Vector(26.6842650, 28.9807890, 24.8783733), verbose
    )


if __name__ == "__main__":
    test_box(True)
