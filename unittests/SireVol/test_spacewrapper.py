try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.Move import *
from Sire.System import *
from Sire.MM import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.CAS import *
from Sire.Units import *


def test_wrapper(verbose=False):
    s = MoleculeParser.read("../io/kcl.crd", "../io/kcl.top")

    space = s.property("space")

    all = s[MGName("all")]

    ff = InterFF("cljff")
    ff.setCLJFunction(CLJShiftFunction(7.5 * angstrom))
    ff.add(all)

    ff2 = InterFF("cljff2")
    ff2.setCLJFunction(CLJShiftFunction(7 * angstrom))
    ff2.add(all)

    s.add(ff)
    s.add(ff2)

    s.setProperty("space", space)

    slow = ff.components().total()
    fast = ff2.components().total()

    slow_s = Symbol("E_slow")
    fast_s = Symbol("E_fast")

    s.setComponent(slow_s, slow)
    s.setComponent(fast_s, fast)
    s.setComponent(s.totalComponent(), slow_s)

    s.add(SpaceWrapper(Vector(0), all))

    move = RigidBodyMC(all)
    vmove = VolumeMove(all)
    moves = WeightedMoves()
    moves.add(move, 10)
    moves.add(vmove, 2)

    mtsmc = MTSMC(moves, 10)
    mtsmc.setFastEnergyComponent(fast_s)
    mtsmc.setSlowEnergyComponent(slow_s)

    moves = SameMoves(mtsmc)

    if verbose:
        print("Moving the molecules...")

    for i in range(0, 10):
        if verbose:
            print("Move %d" % (i + 1))

        s = moves.move(s, 1)

        if verbose:
            print(moves)

        p = PDB2(s)
        p.writeToFile("test%03d.pdb" % (i + 1))

    if verbose:
        print("Moves complete")


if __name__ == "__main__":
    test_wrapper(True)
