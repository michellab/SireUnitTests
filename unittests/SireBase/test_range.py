
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.Base import *

from nose.tools import assert_equal,assert_true

def _compare_equal( r1, r2, verbose ):

    if verbose:
        print("\nCompare %s versus %s" % (r1,r2))

    for val2 in r2:
        val1 = r1.next()

        if verbose:
            print("%s vs %s" % (val1,val2))
 
        assert_equal( val1, val2 )
 
    assert_true( r1.atEnd() )

def test_simple_range(verbose=False):

    r = SimpleRange(0,10)

    _compare_equal( r.populate(100), range(0,10), verbose )
    _compare_equal( r.populate(50), range(0,10), verbose )

    r = SimpleRange(10,0,-1)

    _compare_equal( r.populate(100), range(10,0,-1), verbose)
    _compare_equal( r.populate(11), range(10,0,-1), verbose)

    r = SimpleRange(0,-1)

    _compare_equal( r.populate(100), range(0,100), verbose )
    _compare_equal( r.populate(20), range(0,20), verbose )

    r = SimpleRange(-1, -5, -1)

    _compare_equal( r.populate(11), range(10,6,-1), verbose ) 

if __name__ == "__main__":
    test_simple_range(True)

    
