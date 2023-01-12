
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.MM import *
from Sire.Base import *
from Sire.Mol import *
from Sire.Units import *

from nose.tools import assert_equal

def test_cast(verbose=False):
    c = StringProperty( VariantProperty("hello") )
    assert_equal( c, StringProperty("hello") )

    c = NumberProperty( VariantProperty(1.0) )
    assert_equal( c, NumberProperty(1.0) )

    c = NumberProperty( VariantProperty(5) )
    assert_equal( c, NumberProperty(5) )

    c = BooleanProperty( VariantProperty(False) )
    assert_equal( c, BooleanProperty(False) )

    c = BooleanProperty( VariantProperty(True) )
    assert_equal( c, BooleanProperty(True) )

def test_ff_cast(verbose=False):
    ff = InternalFF()
    ff.setProperty("combiningRules", wrap("geometric"))
    
    if verbose:
        print(ff.property("combiningRules"))

    assert_equal( ff.property("combiningRules").asAString(), "geometric" )

def test_mol_property(verbose=False):

    m = Molecule()

    m = m.edit().setProperty("test1", wrap("true")) \
                .setProperty("test2", wrap("off")) \
                .setProperty("test3", wrap(5.4)) \
                .setProperty("test4", wrap(42)) \
                .setProperty("test5", wrap(42.0)) \
                .setProperty("test6", wrap("Hello World")) \
                .setProperty("test7", wrap( 5 * angstrom )) \
         .commit()

    assert_equal( m.property("test1").asABoolean(), True )
    assert_equal( m.property("test2").asABoolean(), False )
    assert_equal( m.property("test3").asADouble(), 5.4 )
    assert_equal( m.property("test4").asAnInteger(), 42 )
    assert_equal( m.property("test5").asAnInteger(), 42 )
    assert_equal( m.property("test6").asAString(), "Hello World" )

    try:
        assert_equal( m.property("test7").value(), 5 * angstrom )
    except Exception:
        assert_equal( m.property("test7"), 5 * angstrom )

if __name__ == "__main__":
    test_cast(True)
    test_ff_cast(True)
    test_mol_property(True)

