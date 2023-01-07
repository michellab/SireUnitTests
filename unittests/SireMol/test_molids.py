try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.ID import *

from nose.tools import assert_equal, assert_almost_equal

s = MoleculeParser.read("../io/NA16.top", "../io/NA16.rst")

protein = s[MolWithResID("ALA")]


def _test_atomid(mol, id, verbose):
    natoms = mol.nAtoms()
    mass = mol.evaluate().mass().value()
    charge = mol.evaluate().charge().value()

    try:
        nats_id = mol.atoms(id).count()
        mass_id = mol.atoms(id).evaluate().mass().value()
        chg_id = mol.atoms(id).evaluate().charge().value()
    except:
        nats_id = 0
        mass_id = 0
        chg_id = 0

    try:
        nats_neg_id = mol.atoms(id.inverse()).count()
        mass_neg_id = mol.atoms(id.inverse()).evaluate().mass().value()
        chg_neg_id = mol.atoms(id.inverse()).evaluate().charge().value()
    except:
        nats_neg_id = 0
        mass_neg_id = 0
        chg_neg_id = 0

    if verbose:
        print(
            "\n%s : natoms = %s, mass = %s, charge = %s"
            % (id, nats_id, mass_id, chg_id)
        )
        print(
            "%s : natoms = %s, mass = %s, charge = %s"
            % (id.inverse(), nats_neg_id, mass_neg_id, chg_neg_id)
        )
        print(
            "Total : natoms = %s, mass = %s, charge = %s"
            % (natoms, mass, charge)
        )

    assert_equal(natoms, nats_id + nats_neg_id)
    assert_almost_equal(mass, mass_id + mass_neg_id, 3)
    assert_almost_equal(charge, chg_id + chg_neg_id, 3)


def test_atomids(verbose=False):

    ids = [
        AtomName("CA"),
        AtomName("CA") * AtomName("HA"),
        ResName("GLY") + AtomID.any(),
        (ResName("ASP") * ResName("ALA")) + AtomID.any(),
        ResName("HIS") + AtomName("CA"),
        ResName("ALA") + (AtomName("CA") * AtomName("CB")),
        AtomName("ca", CaseInsensitive),
        (ResName("ALA") * ResName("GLY")) + (AtomName("CA") * AtomName("HA")),
    ]

    for id in ids:
        _test_atomid(protein, id, verbose)


if __name__ == "__main__":
    test_atomids(True)
