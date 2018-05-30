
from Sire.IO import *
from Sire.Base import *

from nose.tools import assert_equal

def test_broken_prm(verbose=False):
    if verbose:
        print("Reading the file...")

    a = AmberPrm("../io/hsp90.top")

    if verbose:
        print("Creating the system in serial...")

    s1 = a.toSystem( {"parallel":wrap(False)} )

    if verbose:
        print("Creating the system in parallel...")

    s2 = a.toSystem()

    if verbose:
        print("Comparing the systems...")

    
    assert_equal( s1, s2 )

if __name__ == "__main__":
    test_broken_prm(True)
