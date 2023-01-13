####################################################################################################
#                                                                                                  #
#          Sire script to test that LJ correction works correctly                                  #
#                                                                                                  #
####################################################################################################


try:
    import sire
    sire.use_old_api()
except ImportError:
    pass

from Sire.Tools import LJcutoff
from Sire.Tools import readParams
from Sire.Units import *

from nose.tools import assert_almost_equal

params = {}
params = readParams("../io/lj_cor/sim.cfg")

def test_lj_cor(verbose=False):
    """Check that LJ correction gives expected dG and SD"""

    dg, _ = LJcutoff.runLambda(params) # Ignore SD as this will be 0 for single frame

    if verbose:
        print("#######################################")
        print("Testing LJ Correction...")
        print(f"LJ correction deltaG = {dg}")
        print("#######################################")

    assert_almost_equal(dg.value(), -428.6141143872228)

if __name__ == '__main__':
    test_lj_cor(True)