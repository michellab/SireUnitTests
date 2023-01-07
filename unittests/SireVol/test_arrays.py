try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.Vol import *
from Sire.IO import *

import copy


def test_arrays(verbose=False):

    coords = (
        PDB()
        .readMolecule("../io/sb2.pdb")
        .molecule()
        .property("coordinates")
        .array()
    )

    coords2 = copy.copy(coords)

    if verbose:
        print((coords == coords2))

        print(coords.nCoordGroups(), coords.nCoords())

    coords.remove(1, 2)

    assert coords.nCoordGroups() == coords2.nCoordGroups() - 2
    assert (
        coords.nCoords()
        == coords2.nCoords() - coords2[1].count() - coords2[2].count()
    )

    if verbose:
        print(coords.nCoordGroups(), coords.nCoords())
        print(coords2.nCoordGroups(), coords2.nCoords())

        print((coords == coords2))

    ngroups = coords.nCoordGroups()
    ncoords = coords.nCoords()

    coords.append(coords2)

    assert coords.nCoordGroups() == ngroups + coords2.nCoordGroups()
    assert coords.nCoords() == ncoords + coords2.nCoords()

    if verbose:
        print(coords.nCoordGroups(), coords.nCoords())

    while coords.nCoordGroups() != coords2.nCoordGroups():
        coords.remove(0)

    assert coords == coords2

    if verbose:
        print(coords.nCoordGroups(), coords.nCoords())

    coords.append(coords)
    coords.append(coords)

    assert coords.nCoordGroups() == coords2.nCoordGroups() * 4
    assert coords.nCoords() == coords2.nCoords() * 4

    for i in range(0, coords2.nCoordGroups()):
        assert coords[i] == coords[i + coords2.nCoordGroups()]
        assert coords[i] == coords[i + 2 * coords2.nCoordGroups()]
        assert coords[i] == coords[i + 3 * coords2.nCoordGroups()]

    if verbose:
        print(coords.nCoordGroups(), coords.nCoords())

    coords.remove(0, coords.nCoordGroups())

    if verbose:
        print(coords.nCoordGroups(), coords.nCoords())

    assert coords.isEmpty()


if __name__ == "__main__":
    test_arrays(True)
