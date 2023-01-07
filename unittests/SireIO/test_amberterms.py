try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.Maths import *

from nose.tools import assert_almost_equal

rangen = RanGenerator()

mol = MoleculeParser.read("../io/ose.top", "../io/ose.crd")[
    MolWithResID(ResName("OSE"))
].molecule()


def _assert_expressions_equal(ex0, ex1, symbol):
    try:
        for i in range(0, 100):
            values = {symbol: rangen.rand(-5.0, 5.0)}

            assert_almost_equal(ex0.evaluate(values), ex1.evaluate(values))
    except Exception as e:
        print("FAILED: %s != %s" % (ex0, ex1))
        raise e


def test_bonds(verbose=False):
    bonds = mol.property("bond").potentials()

    for bond in bonds:
        pot = bond.function()
        amberbond = AmberBond(pot, Symbol("r"))

        if verbose:
            print("%s == %s ??" % (amberbond, pot))

        _assert_expressions_equal(
            pot, amberbond.toExpression(Symbol("r")), Symbol("r")
        )


def test_angles(verbose=False):
    angles = mol.property("angle").potentials()

    for angle in angles:
        pot = angle.function()
        amberangle = AmberAngle(pot, Symbol("theta"))

        if verbose:
            print("%s == %s ??" % (amberangle, pot))

        _assert_expressions_equal(
            pot, amberangle.toExpression(Symbol("theta")), Symbol("theta")
        )


def test_dihedrals(verbose=False):
    dihedrals = mol.property("dihedral").potentials()

    for dihedral in dihedrals:
        pot = dihedral.function()
        amberdihedral = AmberDihedral(pot, Symbol("phi"))

        if verbose:
            print("%s == %s ??" % (amberdihedral, pot))

        _assert_expressions_equal(
            pot, amberdihedral.toExpression(Symbol("phi")), Symbol("phi")
        )


def test_impropers(verbose=False):
    impropers = mol.property("improper").potentials()

    for improper in impropers:
        pot = improper.function()
        amberdihedral = AmberDihedral(pot, Symbol("phi"))

        if verbose:
            print("%s == %s ??" % (amberdihedral, pot))

        _assert_expressions_equal(
            pot, amberdihedral.toExpression(Symbol("phi")), Symbol("phi")
        )


def test_dihedral_forms(verbose=False):

    # The symbol for the expression.
    Phi = Symbol("Phi")

    # First test single terms expressions.

    # 1. Standard form.
    f = Expression(0.3 * (1 + Cos(4 * Phi - 2.0)))  # Functional form.
    d = AmberDihedral(f, Phi)  # The Amber dihedral object.
    t = d[0]  # The dihedral part (term).

    # Make sure k, periodicity, and phase are correct.
    assert_almost_equal(t.k(), 0.3)
    assert t.periodicity() == 4
    assert_almost_equal(t.phase(), 2.0)

    # 2. Negative k.
    f = Expression(-0.3 * (1 + Cos(4 * Phi - 2.0)))  # Functional form.
    d = AmberDihedral(f, Phi)  # The Amber dihedral object.
    t = d[0]  # The dihedral part (term).

    # Make sure k, periodicity, and phase are correct.
    assert_almost_equal(t.k(), -0.3)
    assert t.periodicity() == 4
    assert_almost_equal(t.phase(), 2.0)

    # 3. k [ 1 - cos(Phi - phase) ]
    f = Expression(0.3 * (1 - Cos(4 * Phi - 2.0)))  # Functional form.
    d = AmberDihedral(f, Phi)  # The Amber dihedral object.
    t = d[0]  # The dihedral part (term).

    # Make sure k, periodicity, and phase are correct.
    assert_almost_equal(t.k(), 0.3)
    assert t.periodicity() == 4
    assert_almost_equal(t.phase(), 2.0 - pi)

    # Now check series expressions with different types of term.

    # A regular AMBER style dihedral series. All terms have positive cosine factors.
    f = Expression(0.3 * (1 + Cos(Phi)) + 0.8 * (1 + Cos(4 * Phi)))
    d = AmberDihedral(f, Phi)

    # Assert that the expression is the same.
    assert d.toExpression(Phi) == f

    # Try to create an AmberDihedral from a using terms with positive and
    # negative cosine factors. These appear in the CHARMM force field.
    f = Expression(0.3 * (1 + Cos(Phi)) - 0.8 * (1 + Cos(4 * Phi)))
    d = AmberDihedral(f, Phi)

    # Assert that the expression is the same.
    assert d.toExpression(Phi) == f

    # Try to create an AmberDihedral from a using terms with positive and
    # negative cosine factors, with the negative arising from a representation
    # of the form k [ 1 - Cos(Phi) ], rather than -k [ 1 + Cos(Phi) ]. These
    # can occur in Gromacs.
    f = Expression(0.3 * (1 + Cos(Phi)) + 0.8 * (1 - Cos(4 * Phi)))
    d = AmberDihedral(f, Phi)

    # Create a value object: Phi = 2.0.
    val = Values(SymbolValue(Phi.ID(), 2.0))

    # Make sure both expressions evaluate to the same result.
    assert_almost_equal(f.evaluate(val), d.toExpression(Phi).evaluate(val))

    # Try a three-term expression that mixes all formats.
    f = Expression(
        0.3 * (1 + Cos(Phi))
        - 1.2 * (1 + Cos(3 * Phi))
        + 0.8 * (1 - Cos(4 * Phi))
    )
    d = AmberDihedral(f, Phi)
    assert_almost_equal(f.evaluate(val), d.toExpression(Phi).evaluate(val))

    # Try a three-term expression with equal factors.
    f = Expression(
        0.3 * (1 + Cos(Phi))
        - 0.3 * (1 + Cos(3 * Phi))
        + 0.3 * (1 - Cos(4 * Phi))
    )
    d = AmberDihedral(f, Phi)
    assert_almost_equal(f.evaluate(val), d.toExpression(Phi).evaluate(val))


if __name__ == "__main__":
    test_bonds(True)
    test_angles(True)
    test_dihedrals(True)
    test_dihedral_forms(True)
    test_impropers(True)
