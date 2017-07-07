
from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.Maths import *

from nose.tools import assert_almost_equal

rangen = RanGenerator()

mol = MoleculeParser.read("../io/ose.top","../io/ose.crd") \
        [MolWithResID(ResName("OSE"))].molecule()

def _assert_expressions_equal(ex0, ex1, symbol):
    try:
        for i in range(0,100):
            values = {symbol:rangen.rand(-5.0,5.0)}

            assert_almost_equal( ex0.evaluate(values), ex1.evaluate(values) )
    except Exception as e:
        print("FAILED: %s != %s" % (ex0,ex1))
        raise e

def test_bonds(verbose=False):
    bonds = mol.property("bond").potentials()

    for bond in bonds:
        pot = bond.function()
        amberbond = AmberBond(pot, Symbol("r"))

        if verbose:
             print("%s == %s ??" % (amberbond,pot))

        _assert_expressions_equal( pot, amberbond.toExpression(Symbol("r")),
                                   Symbol("r") )

def test_angles(verbose=False):
    angles = mol.property("angle").potentials()

    for angle in angles:
        pot = angle.function()
        amberangle = AmberAngle(pot, Symbol("theta"))

        if verbose:
             print("%s == %s ??" % (amberangle,pot))

        _assert_expressions_equal( pot, amberangle.toExpression(Symbol("theta")),
                                   Symbol("theta") )

def test_dihedrals(verbose=False):
    dihedrals = mol.property("dihedral").potentials()

    for dihedral in dihedrals:
        pot = dihedral.function()
        amberdihedral = AmberDihedral(pot, Symbol("phi"))

        if verbose:
             print("%s == %s ??" % (amberdihedral,pot))
        
        _assert_expressions_equal( pot, amberdihedral.toExpression(Symbol("phi")),
                                   Symbol("phi") )

def test_impropers(verbose=False):
    impropers = mol.property("improper").potentials()

    for improper in impropers:
        pot = improper.function()
        amberdihedral = AmberDihedral(pot, Symbol("phi"))

        if verbose:
             print("%s == %s ??" % (amberdihedral,pot))

        _assert_expressions_equal( pot, amberdihedral.toExpression(Symbol("phi")),
                                   Symbol("phi") )

if __name__ == "__main__":
    test_bonds(True)
    test_angles(True)
    test_dihedrals(True)
    test_impropers(True)

