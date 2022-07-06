
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Vol import *
from Sire.Base import *
from Sire.CAS import *

try:
    from sr.legacy.Base import PropertyList
except Exception:
    pass

from nose.tools import assert_equal

def assert_array_equal( array0, array1 ):
    assert( len(array0) == len(array1) )

    for i in range(0, len(array0)):
        assert( array0[i] == array1[i] )

def test_wrap(verbose=False):
    water = PDB().readMolecule("../io/water.pdb")

    center = water.evaluate().center()

    dblarray = [ 1.0,2,3,4,5 ]
    intarray = [ 1,2,3,4,5 ]
    vecarray = [ Vector(1), Vector(2), Vector(3) ]
    strarray = [ "cat", "dog", "fish" ]
    x = Symbol("x")
    f = (x+5)**2
    mixarray = [ x, f, 5.3, "hello", Molecule(), [f, "cat", PeriodicBox()] ]

    water = water.edit().setProperty("center", wrap(center)) \
                        .setProperty("dblarray", wrap(dblarray)) \
                        .setProperty("intarray", wrap(intarray)) \
                        .setProperty("vecarray", wrap(vecarray)) \
                        .setProperty("strarray", wrap(strarray)) \
                        .setProperty("type", wrap("ligand")) \
                        .setProperty("alpha", wrap(0.5)) \
                        .setProperty("copies", wrap(1)) \
                        .setProperty("mix", wrap(mixarray)).commit()

    assert_equal( water.property("center").value(), center )
    assert_array_equal( water.property("dblarray").value(), dblarray )
    assert_array_equal( water.property("intarray").value(), intarray )
    assert_array_equal( water.property("vecarray").value(), vecarray )
    assert_array_equal( water.property("strarray").value(), strarray )
    assert_equal( water.property("type").value(), "ligand" )
    assert_equal( water.property("alpha").value(), 0.5 )
    assert_equal( water.property("copies").value(), 1 )

    p = water.property("mix")

    assert_equal( Expression(x), p[0].value() )
    assert_equal( f, p[1].value() )
    assert_equal( p[2].value(), 5.3 )
    assert_equal( p[3].value(), "hello" )
    assert_equal( p[4], Molecule() )
    assert_equal( p[5][0].value(), f )
    assert_equal( p[5][1].value(), "cat" )
    assert_equal( p[5][2], PeriodicBox() )

if __name__ == "__main__":
    test_wrap(True)
