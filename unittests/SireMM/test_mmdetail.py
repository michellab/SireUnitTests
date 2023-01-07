
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.MM import *
from Sire.FF import *

from nose.tools import assert_equal

amber14coul = 1.0 / 1.2
amber14lj = 0.5

try:
    amber = MMDetail()
    supported = True
except:
    supported = False

if supported:
    amber = MMDetail( "amber::ff", "arithmetic", amber14coul, amber14lj,
                      "coulomb", "lj", "harmonic", "harmonic", "cosine" )

def test_amber(verbose=False):
    if not supported:
        return

    if verbose:
        print("Testing 'amber' variable\n%s" % amber)

    assert_equal( amber.name(), "amber::ff" )
    assert_equal( amber.combiningRules(), "arithmetic" )
    assert_equal( amber.usesArithmeticCombiningRules(), True )
    assert_equal( amber.usesGeometricCombiningRules(), False )
    assert_equal( amber.electrostatic14ScaleFactor(), amber14coul )
    assert_equal( amber.vdw14ScaleFactor(), amber14lj )
    assert_equal( amber.electrostaticStyle(), "coulomb" )
    assert_equal( amber.usesCoulombCharges(), True )
    assert_equal( amber.vdwStyle(), "lj" )
    assert_equal( amber.usesLJTerm(), True )
    assert_equal( amber.bondStyle(), "harmonic" )
    assert_equal( amber.usesHarmonicBonds(), True )
    assert_equal( amber.angleStyle(), "harmonic" )
    assert_equal( amber.usesHarmonicAngles(), True )
    assert_equal( amber.dihedralStyle(), "cosine" )
    assert_equal( amber.usesCosineDihedrals(), True )

    assert_equal( "amber::ff" in FFDetail.forcefields(), True )
    assert_equal( FFDetail.get("amber::ff"), amber )


if __name__ == "__main__":
    test_amber(True)

