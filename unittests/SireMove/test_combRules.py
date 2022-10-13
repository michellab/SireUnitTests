#
# Sire script to compare Sire and SOMD single point energies
#

from nose.tools import assert_almost_equal

import sire as sr
sr.use_old_api()

from sire.legacy.Tools import OpenMMMD

# Overwrite default parameters.

OpenMMMD.topfile = OpenMMMD.Parameter("topfile", "../io/combRules/toluene_methane.top",
                    """File name of the topology file containing the system to be simulated.""")

OpenMMMD.crdfile = OpenMMMD.Parameter("crdfile", "../io/combRules/toluene_methane.crd",
                    """File name of the coordinate file containing the coordinates of the
                       system to be simulated.""")

OpenMMMD.morphfile = OpenMMMD.Parameter("morphfile", "../io/combRules/toluene_methane.pert",
                      """Name of the morph file containing the perturbation to apply to the system.""")

OpenMMMD.andersen = OpenMMMD.Parameter("thermostat", False,
      """Whether or not to use the Andersen thermostat (needed for NVT or NPT simulation).""")

OpenMMMD.barostat = OpenMMMD.Parameter("barostat", False, """Whether or not to use a barostat (needed for NPT simulation).""")

OpenMMMD.cutoff_type = OpenMMMD.Parameter("cutoff type", "nocutoff", """The cutoff method to use during the simulation.""")

OpenMMMD.constraint = OpenMMMD.Parameter("constraint", "none", """The constraint model to use during dynamics.""")

OpenMMMD.platform = OpenMMMD.Parameter("platform", "Reference", """Which OpenMM platform should be used to perform the dynamics.""")

def test_combining_rules(rule, verbose=False):
    # Set the combining rule.
    OpenMMMD.combining_rules = OpenMMMD.Parameter("combining rules", rule,
                            """Combining rules to use for the non-bonded interactions.""")

    amber = OpenMMMD.Amber()

    (molecules, space) = amber.readCrdTop(OpenMMMD.crdfile.val, OpenMMMD.topfile.val)
    system = OpenMMMD.createSystemFreeEnergy(molecules)

    system = OpenMMMD.setupForceFieldsFreeEnergy(system, space)
    if OpenMMMD.debug_seed.val:
        ranseed = OpenMMMD.debug_seed.val
    else:
        ranseed = OpenMMMD.RanGenerator().randInt(100000, 1000000)

    moves = OpenMMMD.setupMovesFreeEnergy(system, ranseed, 
            OpenMMMD.gpu.val, OpenMMMD.lambda_val.val)

    # Get energy from Sire
    nrg_sire = system.energy()
    # Get energy from SOMD
    mdmoves = moves.moves()[0]
    integrator = mdmoves.integrator()
    nrg_somd = integrator.getPotentialEnergy(system)

    nrg_diff = (nrg_sire - nrg_somd).value()

    if verbose:
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        if rule == 'arithmetic':
            print ("TESTING ARITHMETIC COMBINING RULES")
        else:
            print ("TESTING GEOMETRIC COMBINING RULES")
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print ("* The energy from Sire is: %s" %nrg_sire)
        print ("* The energy from SOMD is: %s" % nrg_somd)

        if rule == 'arithmetic':
            print (" For the arithmetic combining rules the single point energy difference between sire and somd at lambda %s is %s " % (OpenMMMD.lambda_val.val, nrg_diff) )
        else:
            print (" For the geometric combining rules the single point energy difference between sire and somd at lambda %s is %s " % (OpenMMMD.lambda_val.val, nrg_diff) )

    assert_almost_equal(nrg_diff, 0.0, 4)

if __name__ == '__main__':
    test_combining_rules("arithmetic", True)
    test_combining_rules("geometric", True)
