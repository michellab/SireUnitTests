
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.Base import *

from nose.tools import assert_true,assert_false,assert_equal

def test_hasvalue(verbose=False):
    p = PropertyName("hello")
    assert_false( p.hasValue() )
    assert_true( p.hasSource() )

    m = PropertyMap()
    m.set("welcome", "hello")
    assert_false( m["welcome"].hasValue() )
    assert_true( m["welcome"].hasSource() )

    s = wrap("hello")
    assert_equal( s.__class__, StringProperty )

    m = PropertyMap()
    m.set("welcome", wrap("hello"))
    assert_true( m["welcome"].hasValue() )
    assert_false( m["welcome"].hasSource() )

    s = StringProperty("hello")

    p = PropertyName(s)
    assert_false( p.hasSource() )
    assert_true( p.hasValue() )

    m = PropertyMap()
    m.set("welcome", s)

    assert_false( m["welcome"].hasSource() )
    assert_true( m["welcome"].hasValue() )


if __name__ == "__main__":
    test_hasvalue(True)

