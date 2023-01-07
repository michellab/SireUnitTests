try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.System import *
from Sire.Vol import *
from Sire.Base import *
from Sire.FF import *
from Sire.Maths import *
from Sire.MM import *

from nose.tools import assert_equal, assert_almost_equal


def test_props(verbose=False):
    sys = System()

    box0 = PeriodicBox(Vector(10.0, 10.0, 10.0))
    box1 = PeriodicBox(Vector(20.0, 20.0, 20.0))

    if verbose:
        print(box0)
        print(box0.volume())
        print(box1.volume())

    assert not sys.containsProperty("space")

    sys.add(InterCLJFF("cljff"))

    if verbose:
        print(sys)
        print(sys.property("space"))
        print(sys.userProperties().propertyKeys())
        print(sys.builtinProperties().propertyKeys())

    assert sys.containsProperty("space")
    assert_equal(sys.property("space"), Cartesian())

    sys.setProperty("space0", LinkToProperty("space", FFIdx(0)))

    if verbose:
        print(sys.property("space0"))

    assert sys.containsProperty("space0")

    sys.setProperty("space0", box0)

    if verbose:
        print(sys.property("space"))

    assert_equal(sys.property("space0"), box0)

    sys.setProperty("space1", box1)

    sys.setProperty("combined_space", CombineSpaces("space0", "space1"))

    assert_equal(sys.property("space1"), box1)

    if verbose:
        print(sys.properties().propertyKeys())

        print(sys.property("combined_space"))
        print(sys.property("combined_space").volume())

    assert_almost_equal(
        sys.property("combined_space").volume().value(),
        sys.property("space0").volume().value()
        + sys.property("space1").volume().value(),
        5,
    )

    space3 = PeriodicBox(Vector(5, 5, 5))
    sys.setProperty("space0", space3)

    assert_equal(sys.property("space0"), space3)

    if verbose:
        print(sys.property("combined_space"))
        print(sys.property("combined_space").volume())

    assert_almost_equal(
        sys.property("combined_space").volume().value(),
        sys.property("space0").volume().value()
        + sys.property("space1").volume().value(),
        5,
    )

    sys.removeProperty("space0")

    if verbose:
        print(sys.properties().propertyKeys())

    assert not sys.containsProperty("space0")


if __name__ == "__main__":
    test_props(True)
