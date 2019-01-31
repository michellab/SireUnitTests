from Sire.Maths import *
from Sire.Mol import *

from nose.tools import assert_equal

def test_unnamed(verbose=False):
    na = Molecule()
    na = na.edit().add( CGName("1") ).add( AtomName("Na") ).molecule().commit()

    na = na.edit().atom(AtomIdx(0)) \
                  .setProperty("coordinates", Vector(0,0,0)) \
                  .molecule() \
                  .commit()

    assert_equal(na.nAtoms(), 1)

if __name__ == "__main__":
    test_unnamed(True)


