try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.Vol import *
from Sire.Maths import *
from Sire.Units import *

from nose.tools import assert_equal


def test_grid(verbose=True):

    r = RegularGrid(Vector(1, 2, 3), 5, 2 * angstrom)

    if verbose:
        print(r)
        print(r.center())
        print(r.gridSpacing())
        print(r.points())

    r2 = r.rotate(Quaternion(32 * degrees, Vector(1, 0, 0)), r.center())

    if verbose:
        print(r2)
        print(r2.center())
        print(r2.gridSpacing())
        print(r2.points())

    assert_equal(r.center(), r2.center())
    assert_equal(r.gridSpacing(), r2.gridSpacing())
    assert_equal(len(r.points()), len(r2.points()))


if __name__ == "__main__":
    test_grid(True)
