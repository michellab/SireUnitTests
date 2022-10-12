from sire.legacy.Base import *
from sire.legacy.Stream import *
from sire.legacy.Qt import *

from nose.tools import assert_equal

def test_sharedstream(verbose=False):

    props = Properties()

    p = StringProperty("mieow")

    props.setProperty( "cat", p )
    props.setMetadata( "cat", p )
    props.setProperty( "tiger", p )
    props.setMetadata( "tiger", p )

    data = save(props)

    if verbose:
        print(props)
        print(("Properties takes up %d bytes" % data.size()))
        print("Loading the properties...")

    p2 = load(data)

    if verbose:
        print("...complete!")
        print(props)
        print(p2)

        print((props.property("cat")))
        print((p2.property("cat")))

    assert_equal( props.property("cat"), p2.property("cat") )
    assert_equal( props.metadata("cat"), p2.metadata("cat") )

if __name__ == "__main__":
    test_sharedstream(True)

