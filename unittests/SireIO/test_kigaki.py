try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.CAS import *

from nose.tools import assert_equal, assert_almost_equal


def _same_function(f1, f2, x, start, end):
    delta = (end - start) / 25

    for i in range(0, 26):
        val = {x: (start + i * delta)}
        v1 = f1.evaluate(val)
        v2 = f2.evaluate(val)

        if abs(v2 - v1) > 0.001:
            return False

    return True


def _same_atoms(a1, a2):
    if a1.atom0() != a2.atom0() or a1.atom1() != a1.atom1():
        return False

    try:
        if a1.atom2() != a2.atom2():
            return False

        if a1.atom3() != a2.atom3():
            return False
    except:
        pass

    return True


def _assert_almost_equal(a, b, tol):
    if abs(a - b) > 0.1:
        assert_almost_equal(a, b, tol)


def _getEnergies(s):
    intraclj = IntraFF("intraclj")
    intraclj.add(s.molecules())

    intraff = InternalFF("intraff")
    intraff.setUse14Calculation(True)
    intraff.add(s.molecules())

    interff = InterFF("interff")
    interff.add(s.molecules())

    ffs = ForceFields()
    ffs.add(intraclj)
    ffs.add(intraff)
    ffs.add(interff)

    return ffs.energies()


def _assert_all_almost_equal(oldnrgs, newnrgs):
    oldkeys = list(oldnrgs.keys())
    newkeys = list(newnrgs.keys())

    oldkeys.sort()
    newkeys.sort()

    assert_equal(oldkeys, newkeys)

    for key in oldkeys:
        _assert_almost_equal(oldnrgs[key], newnrgs[key], 3)


def _printCompareEnergies(oldnrgs, newnrgs):
    keys = list(oldnrgs.keys())
    keys.sort()

    for key in keys:
        print("%s: %s  %s" % (key, oldnrgs[key], newnrgs[key]))


def _compare_dihedrals(d1, d2, verbose):
    ps1 = []
    ps2 = []

    # only compare non-zero dihedrals
    for p in d1.potentials():
        if not p.function().isZero():
            ps1.append(p)

    for p in d2.potentials():
        if not p.function().isZero():
            ps2.append(p)

    assert_equal(len(ps1), len(ps2))

    for i in range(0, len(ps1)):
        p1 = ps1[i]
        p2 = ps2[i]

        # find the matching atomids...
        if not _same_atoms(p1, p2):
            p2 = None

            for i in range(0, len(ps2)):
                if _same_atoms(p1, ps2[i]):
                    p2 = ps2[i]
                    break

        if p2 is None:
            if verbose:
                print("Cannot find the dihedral matching '%s'" % p1)
            assert_equal(p1, p2)

        if not _same_function(
            p1.function(), p2.function(), Symbol("phi"), 0, 3.141
        ):
            if verbose:
                print("FAILED: %s != %s" % (p1, p2))

            assert_equal(p1.function(), p2.function())


def _atom(m, i):
    a = m.atom(AtomIdx(i))
    r = a.residue()
    return "%s:%s:%s" % (
        a.name().value(),
        r.name().value(),
        r.number().value(),
    )


def _print_different_scls(scl1, scl2, mol1):
    for i in range(0, scl1.nAtoms()):
        for j in range(0, scl2.nAtoms()):
            s1 = scl1.get(AtomIdx(i), AtomIdx(j))
            s2 = scl2.get(AtomIdx(i), AtomIdx(j))
            if s1 != s2:
                print(
                    "%s %s : %s %s" % (_atom(mol1, i), _atom(mol1, j), s1, s2)
                )


def _compare_nbscl(scl1, scl2, m, verbose):
    assert_equal(scl1.nAtoms(), scl2.nAtoms())

    for i in range(0, scl1.nAtoms()):
        for j in range(0, scl2.nAtoms()):
            s1 = scl1.get(AtomIdx(i), AtomIdx(j))
            s2 = scl2.get(AtomIdx(i), AtomIdx(j))

            if s1 != s2:
                if verbose:
                    print("FAIL: %s %s" % (i, j))
                    _print_different_scls(scl1, scl2, m)

                assert_equal(s1, s2)


def test_kigaki(verbose=False):

    if verbose:
        print("Reading files...")

    s = MoleculeParser.read(
        "../io/kigaki.top",
        "../io/kigaki.gro",
        {"GROMACS_PATH": "../io/gromacs"},
    )

    if verbose:
        print("...done! Calculating the energies...")

    # Gromacs 5.1.4 gives the energies as
    # Energy                      Average   Err.Est.       RMSD  Tot-Drift
    # -------------------------------------------------------------------------------
    # Bond                        1537.58         --          0          0  (kJ/mol)
    # Angle                       304.341         --          0          0  (kJ/mol)
    # Proper Dih.                 598.914         --          0          0  (kJ/mol)
    # Improper Dih.               14.9692         --          0          0  (kJ/mol)
    # LJ-14                       312.687         --          0          0  (kJ/mol)
    # Coulomb-14                  5687.13         --          0          0  (kJ/mol)
    # LJ (SR)                     87642.7         --          0          0  (kJ/mol)
    # Coulomb (SR)                -183099         --          0          0  (kJ/mol)

    nrgs = _getEnergies(s)
    bond = nrgs(Symbol("E_{intraff}^{bond}")) * kcal_per_mol
    angle = nrgs(Symbol("E_{intraff}^{angle}")) * kcal_per_mol
    dihedral = nrgs(Symbol("E_{intraff}^{dihedral}")) * kcal_per_mol
    lj14 = nrgs(Symbol("E_{intraff}^{1-4[LJ]}")) * kcal_per_mol
    coul14 = nrgs(Symbol("E_{intraff}^{1-4[coulomb]}")) * kcal_per_mol

    gbond = 1537.58 * kJ_per_mol
    gangle = 304.341 * kJ_per_mol
    gdihedral = (598.914 + 14.9692) * kJ_per_mol
    glj14 = 312.687 * kJ_per_mol
    gcoul14 = 5687.13 * kJ_per_mol

    if verbose:
        print("Bond: %s versus %s" % (bond, gbond))
        print("Angle: %s versus %s" % (angle, gangle))
        print("Dihedral: %s versus %s" % (dihedral, gdihedral))
        print(
            "Dihedral: %s kJ mol-1 versus %s kJ mol-1"
            % (dihedral.to(kJ_per_mol), gdihedral.to(kJ_per_mol))
        )
        print("LJ14: %s versus %s" % (lj14, glj14))
        print("COUL14: %s versus %s" % (coul14, gcoul14))

    _assert_almost_equal(bond.value(), gbond.value(), 2)
    _assert_almost_equal(angle.value(), gangle.value(), 2)
    _assert_almost_equal(dihedral.value(), gdihedral.value(), 2)
    _assert_almost_equal(lj14.value(), glj14.value(), 2)
    _assert_almost_equal(coul14.value(), gcoul14.value(), 2)

    files = MoleculeParser.write(s, ["test.prm7", "test.rst", "test.pdb"])

    if verbose:
        print("Written files %s" % files)
        print("Loading back up from amber...")

    s2 = MoleculeParser.read("test.prm7", "test.rst")

    if verbose:
        print("Comparing intramolecular terms...")

    for i in range(0, s.nMolecules()):
        mol1 = s[MolIdx(i)]
        mol2 = s2[MolIdx(i)]
        _compare_dihedrals(
            mol1.property("dihedral"), mol2.property("dihedral"), verbose
        )
        _compare_nbscl(
            mol1.property("intrascale"),
            mol2.property("intrascale"),
            mol1,
            verbose,
        )

    if verbose:
        print("Calculating energies...")

    nrgs2 = _getEnergies(s2)

    if verbose:
        _printCompareEnergies(nrgs, nrgs2)

    _assert_all_almost_equal(nrgs, nrgs2)

    # write back to gromacs...
    if verbose:
        print("Writing back to gromacs format...")

    g = GroTop(s2)
    g.writeToFile("test.grotop")

    # read back from gromacs...
    if verbose:
        print("Reading back from gromacs format...")

    s3 = MoleculeParser.read("test.grotop", "test.rst")

    if verbose:
        print("Recalculating energies...")

    nrgs3 = _getEnergies(s3)

    if verbose:
        _printCompareEnergies(nrgs, nrgs3)

    _assert_all_almost_equal(nrgs, nrgs3)


if __name__ == "__main__":
    test_kigaki(True)
