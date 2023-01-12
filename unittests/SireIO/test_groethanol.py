try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.FF import *
from Sire.CAS import *
from Sire.Units import *

from nose.tools import assert_almost_equal


def _assert_almost_equal(a, b, tol):
    if abs(a - b) > 0.1:
        assert_almost_equal(a, b, tol)


def _getEnergies(s):
    intraclj = IntraFF("intraclj")
    intraclj.add(s.molecules())

    intraff = InternalFF("intraff")
    intraff.setUse14Calculation(True)
    # intraff.setCombiningRules( CLJFunction.GEOMETRIC )
    intraff.add(s.molecules())

    interff = InterFF("interff")
    interff.add(s.molecules())

    ffs = ForceFields()
    ffs.add(intraclj)
    ffs.add(intraff)
    ffs.add(interff)

    return ffs.energies()


def test_ethanol(verbose=False):

    s = MoleculeParser.read("../io/ethanol.grotop", "../io/ethanol.gro")

    nrgs = _getEnergies(s)

    # Expected energies from Gromacs 5.1.4
    # Energy                      Average   Err.Est.       RMSD  Tot-Drift
    # -------------------------------------------------------------------------------
    # Bond                       0.193339         --          0          0  (kJ/mol)
    # Angle                       23.5227         --          0          0  (kJ/mol)
    # Ryckaert-Bell.              5.73895         --          0          0  (kJ/mol)
    # LJ-14                       2.66035         --          0          0  (kJ/mol)
    # Coulomb-14                 -30.5237         --          0          0  (kJ/mol)
    # LJ (SR)                     5745.08         --          0          0  (kJ/mol)
    # Coulomb (SR)               -41509.9         --          0          0  (kJ/mol)

    gbond = 0.193339 * kJ_per_mol
    gangle = 23.5227 * kJ_per_mol
    gdihedral = 5.73895 * kJ_per_mol
    glj14 = 2.66035 * kJ_per_mol
    gcoul14 = -30.5237 * kJ_per_mol
    glj = 5745.08 * kJ_per_mol
    gcoul = -41509.9 * kJ_per_mol

    bond = nrgs(Symbol("E_{intraff}^{bond}")) * kcal_per_mol
    angle = nrgs(Symbol("E_{intraff}^{angle}")) * kcal_per_mol
    dihedral = nrgs(Symbol("E_{intraff}^{dihedral}")) * kcal_per_mol
    lj14 = nrgs(Symbol("E_{intraff}^{1-4[LJ]}")) * kcal_per_mol
    coul14 = nrgs(Symbol("E_{intraff}^{1-4[coulomb]}")) * kcal_per_mol

    if verbose:
        print("Bond: %s versus %s" % (bond, gbond))
        print("Angle: %s versus %s" % (angle, gangle))
        print("Dihedral: %s versus %s" % (dihedral, gdihedral))
        print("LJ14: %s versus %s" % (lj14, glj14))
        print("COUL14: %s versus %s" % (coul14, gcoul14))

    _assert_almost_equal(bond.value(), gbond.value(), 2)
    _assert_almost_equal(angle.value(), gangle.value(), 2)
    _assert_almost_equal(dihedral.value(), gdihedral.value(), 2)
    _assert_almost_equal(lj14.value(), glj14.value(), 2)
    _assert_almost_equal(coul14.value(), gcoul14.value(), 2)


if __name__ == "__main__":
    test_ethanol(True)
