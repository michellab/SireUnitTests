
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.Mol import *
from Sire.MM import *
from Sire.IO import *

from nose.tools import assert_equal

mol = MoleculeParser.read("../io/kigaki.top",
              {"GROMACS_PATH":"../io/gromacs"})[MolIdx(0)]

scl14 = CLJScaleFactor( 1.0/1.2, 0.5 )

def test_cljpairs(verbose=False):
    if verbose:
        print("Creating a default set...")

    scl = CLJNBPairs( mol.info(), CLJScaleFactor(1,1) )

    assert_equal( scl.nAtoms(), mol.info().nAtoms() )
    assert_equal( scl.nGroups(), mol.info().nCutGroups() )

    # everything should be 1,1
    for i in range(0,scl.nAtoms()):
        for j in range(0,scl.nAtoms()):
            s1 = scl.get( AtomIdx(i), AtomIdx(j) )
            s2 = scl.get( AtomIdx(j), AtomIdx(i) )

            assert_equal( s1, s2 )
            assert_equal( s1, CLJScaleFactor(1) )

    # now set a value...
    scl.set( AtomIdx(0), AtomIdx(24), scl14 )

    for i in range(0,scl.nAtoms()):
        for j in range(0,scl.nAtoms()):
            s1 = scl.get( AtomIdx(i), AtomIdx(j) )
            s2 = scl.get( AtomIdx(j), AtomIdx(i) )

            assert_equal( s1, s2 )

            if (i == 0 and j == 24) or (i == 24 and j == 0):
                assert_equal( s1, scl14 )
            else:
                assert_equal( s1, CLJScaleFactor(1) )

def test_autopairs(verbose=False):
    if verbose:
        print("Checking the auto-generated pairs...")

    conn = mol.property("connectivity")

    scl = CLJNBPairs(conn, scl14)

    assert_equal( scl.nAtoms(), mol.info().nAtoms() )
    assert_equal( scl.nGroups(), mol.info().nCutGroups() )

    # everything should be 1,1
    for i in range(0,scl.nAtoms()):
        for j in range(0,scl.nAtoms()):
            atm0 = AtomIdx(i)
            atm1 = AtomIdx(j)

            s1 = scl.get( atm0, atm1 )
            s2 = scl.get( atm1, atm0 )

            assert_equal( s1, s2 )

            bonded = conn.areBonded(atm0,atm1)
            assert_equal( bonded, conn.areBonded(atm1,atm0) )

            angled = conn.areAngled(atm0,atm1)
            assert_equal( angled, conn.areAngled(atm1,atm0) )

            dihedraled = conn.areDihedraled(atm0,atm1)
            assert_equal( dihedraled, conn.areDihedraled(atm1,atm0) )

            typ = conn.connectionType(atm0,atm1)

            if atm0 == atm1:
                assert_equal( typ, 1 )
            elif bonded:
                assert_equal( typ, 2 )
            elif angled:
                assert_equal( typ, 3 )
            elif dihedraled:
                assert_equal( typ, 4 )
            else:
                assert_equal( typ, 0 )

            if atm0 == atm1 or bonded or angled:
                if verbose and s1 != CLJScaleFactor(0):
                    print("%s %s : %s %s %s : %s" % \
                     (atm0,atm1,bonded,angled,dihedraled,s1))

                assert_equal( s1, CLJScaleFactor(0) )

            elif dihedraled:
                assert_equal( s1, scl14 )

            else:
                assert_equal( s1, CLJScaleFactor(1) )


def test_molpairs(verbose=False):
    if verbose:
        print("Checking the molecule scale factors...")

    scl = mol.property("intrascale")

    assert_equal( scl.nAtoms(), mol.info().nAtoms() )
    assert_equal( scl.nGroups(), mol.info().nCutGroups() )

    # everything should be 1,1
    for i in range(0,scl.nAtoms()):
        for j in range(0,scl.nAtoms()):
            s1 = scl.get( AtomIdx(i), AtomIdx(j) )
            s2 = scl.get( AtomIdx(j), AtomIdx(i) )

            assert_equal( s1, s2 )

    assert_equal( scl, CLJNBPairs(mol.property("connectivity"), scl14) )

if __name__ == "__main__":
    test_autopairs(True)
    test_molpairs(True)
    test_cljpairs(True)

