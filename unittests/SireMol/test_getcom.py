try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.ID import *

from nose.tools import assert_almost_equal


def test_com(verbose=False):

    protein = PDB().readMolecule("../io/cox2.pdb")

    mass = 0
    com = Vector(0)

    for atom in protein.atoms(AtomName("CA", CaseInsensitive)):
        com += (
            atom.property("coordinates")
            * atom.property("element").mass().value()
        )
        mass += atom.property("element").mass().value()

    com /= mass

    if verbose:
        print("Center of mass = %s" % (com))

    t = (
        protein.selectAll(AtomName("CA", CaseInsensitive))
        .evaluate()
        .centerOfMass()
    )

    if verbose:
        print("Center of mass = %s" % (t))

    assert_almost_equal(com.x(), t.x(), 5)
    assert_almost_equal(com.y(), t.y(), 5)
    assert_almost_equal(com.z(), t.z(), 5)


if __name__ == "__main__":
    test_com(True)
