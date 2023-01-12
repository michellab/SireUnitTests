try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

import Sire.Stream

from Sire.Move import *
from Sire.System import *

(system, moves) = Sire.Stream.load("../io/lsrc_system.s3")


def _pvt_test(sys, move, verbose):

    s = System(sys)
    g = move.generator()
    g.seed(42)
    move.setGenerator(g)

    s.mustNowRecalculateFromScratch()

    oldnrgs = s.energies()

    if verbose:
        print("\nTESTING %s" % move)

    for i in range(0, 10):
        if verbose:
            print("MOVE %s" % i)

        move.move(s, 10, False)

    if verbose:
        print("Recalculating energies...")

    newnrgs = s.energies()

    s.mustNowRecalculateFromScratch()

    r_newnrgs = s.energies()

    keys = list(r_newnrgs.keys())
    keys.sort()

    all_agree = True

    if verbose:
        print("MOVE COMPLETE: %s" % move)

    for key in keys:
        if abs(newnrgs[key] - r_newnrgs[key]) > 0.0001:
            if verbose:
                print(
                    "BROKEN %s : %s versus %s"
                    % (key, oldnrgs[key], newnrgs[key])
                )

            all_agree = False

    assert all_agree


def test_moves(verbose=False):

    for move in moves:
        _pvt_test(system, move, verbose)


if __name__ == "__main__":
    test_moves(True)
