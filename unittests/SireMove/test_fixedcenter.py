
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

import Sire.Stream

from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Move import *
from Sire.IO import *
from Sire.Base import *
from Sire.Maths import *
from Sire.Units import *

from nose.tools import assert_almost_equal

def _pvt_test_fixed_center(center_method, verbose = False):
    ligand = Sire.Stream.load("../io/osel.s3")

    if center_method:
        old_center = center_method.getPoint(ligand)
    else:
        ligand = ligand.edit().setProperty("center", wrap(Vector(1.0,2.0,3.0))).commit()
        old_center = ligand.property("center")

    intraff = InternalFF("intraff")
    intraff.add(ligand)

    intraclj = IntraCLJFF("intraclj")
    intraclj.add(ligand)

    system = System()
    system.add(intraff)
    system.add(intraclj)

    mols = MoleculeGroup("mols")
    mols.add(ligand)
    system.add(mols)

    intramove = InternalMove(mols)

    rbmove = RigidBodyMC(mols)
    rbmove.setMaximumTranslation(0*angstrom)

    if center_method:
        rbmove.setCenterOfRotation(center_method)
        intramove.setCenterOfMolecule(center_method)

    moves = WeightedMoves()
    moves.add(intramove, 1)
    moves.add(rbmove, 1)

    for i in range(0,10):
        system = moves.move(system, 25, False)

        if verbose:
            print("Completed 25 moves...")

        ligand = system[ligand.number()][0].molecule()

        if center_method:
            new_center = center_method.getPoint(ligand)
        else:
            new_center = ligand.property("center")

        if verbose:
            print("Old center = %s" % old_center)
            print("New center = %s" % new_center)

        assert_almost_equal( old_center.x(), new_center.x(), 1 )
        assert_almost_equal( old_center.y(), new_center.y(), 1 )
        assert_almost_equal( old_center.z(), new_center.z(), 1 )

def test_fixed_center(verbose=False):
    if verbose:
        print("\nTesting with center of mass...")

    _pvt_test_fixed_center( GetCOMPoint(), verbose )

    if verbose:
        print("\nTesting with centroid...")

    _pvt_test_fixed_center( GetCentroidPoint(), verbose )

    if verbose:
        print("\nTesting with fixed point...")

    _pvt_test_fixed_center( None, verbose )

if __name__ == "__main__":
    test_fixed_center(True)
