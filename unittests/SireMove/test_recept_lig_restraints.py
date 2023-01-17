####################################################################################################
#                                                                                                  #
#          Sire script to test that receptor-ligand restraints are working correctly               #
#                                                                                                  #
####################################################################################################

try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

import os
from nose.tools import assert_almost_equal
from Sire.Tools import OpenMMMD

# Overwrite default parameters

OpenMMMD.topfile = OpenMMMD.Parameter("topfile", "../io/receptor_ligand_restraints/SYSTEM.top",
                    """File name of the topology file containing the system to be simulated.""")

OpenMMMD.crdfile = OpenMMMD.Parameter("crdfile", "../io/receptor_ligand_restraints/SYSTEM.crd",
                    """File name of the coordinate file containing the coordinates of the
                       system to be simulated.""")

OpenMMMD.distance_restraints_dict = OpenMMMD.Parameter("distance restraints dictionary", {(21, 4950): (2.72, 19.55, 0), (11, 4946): (5.92, 
                                    12.75, 0), (3, 4909): (3.63, 9.64, 0), (5, 4949): (6.87, 9.28, 0), (4, 971): (4.17,
                                     7.42, 0), (2, 1613): (8.67, 7.07, 0), (14, 963): (5.56, 6.81, 0), (20, 
                                     959): (4.73, 6.07, 0), (10, 950): (6.68, 5.77, 0), (12, 51): (9.15, 5.64, 0),(18, 
                                     4951): (9.52, 5.31, 0), (13, 47): (6.0, 5.02, 0), (19, 49): (8.75, 4.94, 0), (6,
                                      53): (9.8, 3.92, 0), (17, 34): (4.24, 3.64, 0), (16, 45): (7.8, 3.59, 0), 
                                      (9, 548): (8.14, 2.81, 0), (15, 48): (9.26, 2.74, 0), (7, 584): (8.25, 1.63, 
                                      0), (8, 1633): (10.93, 1.6, 0), (0, 1737): (6.64, 1.28, 0), (1, 4914): (10.4, 1.13, 0)},
                                     """Dictionary of pair of atoms whose distance is restrained, and restraint
                                     parameters. Syntax is {(atom0,atom1):(reql, kl, Dl)} where atom0, atom1 are atomic
                                     indices. reql the equilibrium distance. Kl the force constant of the restraint.
                                     D the flat bottom radius.""")

OpenMMMD.boresch_restraints_dict = OpenMMMD.Parameter("boresch restraints dictionary", {"anchor_points":{"r1":4946, "r2":4944, "r3":4949,
                                     "l1":11, "l2":2, "l3":3},"equilibrium_values":{"r0":5.92, "thetaA0":1.85, 
                                     "thetaB0":1.59,"phiA0":-0.30, "phiB0":-1.55, "phiC0":2.90},
                                     "force_constants":{"kr":12.75, "kthetaA":33.37, "kthetaB":19.20, "kphiA":107.68,
                                      "kphiB":24.62, "kphiC":24.90}},
                                    """Dictionary of four dictionaries: anchor points in ligand, anchor points in receptor,
                                    equilibrium values for 6 Boresch-style external degrees of freedom, and associated force
                                    constants. Syntax is:
                                    {
                                    "anchor_points":{"r1":r1, "r2":r2, "r3":r3, "l1":l1, "l2":l2, "l3":l3},
                                    "equilibrium_values":{"r0":r0, "thetaA0": thetaA0, "thetaB0": thetaB0,
                                                          "phiA0":phiA0, "phiB0": phiB0, "phiC0":phiC0},
                                    "force_constants":{"kr":kr, "kthetaA": kthetaA, "kthetaB": kthetaB,
                                                       "kphiA":kphiA, "kphiB": kphiB, "kphiC":kphiC}
                                    } 
                                    r1 - 3 and l1 - 3 are the anchor points in the receptor and ligand, respectively, 
                                    given by atomic indices. r is | l1 - r1 | (A). thetaA, and thetaB are the angles
                                    (r2, r1, l1) and (r1, l1, l2) (rad). phiA, phiB, and phiC are the dihedral angles
                                    (r3, r2, r1, l1), (r2, r1, l1, l2), and (r1, l1, l2, l3), respectively (rad). A first 
                                    character of k indicates a force constant (kcal mol^-1 A^-2 for the distance and 
                                    kcal mol^-1 rad^-2 for the angles) and a final character of 0 indicates an
                                    equilibrium value (A or rad). To use Boresch restraints, "use boresch restraints" 
                                    must be set equal to True in the config file. 
                                    """)

OpenMMMD.morphfile = OpenMMMD.Parameter("morphfile", "../io/receptor_ligand_restraints/MORPH.dummy.pert",
                      """Name of the morph file containing the perturbation to apply to the system.""")

OpenMMMD.platform = OpenMMMD.Parameter("platform", "CPU", "Override the GPU platform")

######## MAIN SCRIPTS  #############

def getPotEnergyRestr(boresch_on=False, mult_dist_on=False, verbose=False):
    """Function to set up system with or without Boresch/ multiple
    distance restraints.

    Args:
        boresch_on (bool, optional): Whether or not to set up Boresch restraints. Defaults to False.
        mult_dist_on (bool, optional): Whether or not to set up multiple distance restraints. Defaults to False.
        verbose (bool, optional): Whether to print output to screen. Defaults to False.

    Returns:
        float: The energy of the system (kcal / mol)
    """
    # Silence Sire output
    os.environ['SIRE_SILENT_PHONEHOME'] = "1"
    os.environ['OPENMM_DEFAULT_PLATFORM'] = "CPU"

    amber = OpenMMMD.Amber()
    (molecules, space) = amber.readCrdTop(OpenMMMD.crdfile.val, OpenMMMD.topfile.val)
    system = OpenMMMD.createSystemFreeEnergy(molecules)

    if mult_dist_on:
        system = OpenMMMD.setupDistanceRestraints(system, restraints=OpenMMMD.distance_restraints_dict.val)

    if boresch_on:
        system = OpenMMMD.setupBoreschRestraints(system)

    system = OpenMMMD.setupForceFieldsFreeEnergy(system, space)
    moves = OpenMMMD.setupMovesFreeEnergy(system, OpenMMMD.debug_seed.val, 
                                          OpenMMMD.gpu.val, OpenMMMD.lambda_val.val)
    mdmoves = moves.moves()[0]
    integrator = mdmoves.integrator()
    nrg = integrator.getPotentialEnergy(system)
    return nrg


def test_boresch_restraints(verbose=False):
    """Check that Boresch restraints give expected increase in energy"""

    boresch_on_nrg = getPotEnergyRestr(boresch_on=True, verbose=verbose)
    boresch_off_nrg = getPotEnergyRestr(boresch_on=False, verbose=verbose)
    nrg_diff = boresch_on_nrg - boresch_off_nrg

    if verbose:
        print("#######################################")
        print("Testing Boresch restraints...")
        print(f"Energy of system with Boresch restraints on = {boresch_on_nrg}")
        print(f"Energy of system with Boresch restraints off = {boresch_off_nrg}")
        print(f"Energy change when restraints turned on = {nrg_diff}")
        print("#######################################")

    assert_almost_equal(nrg_diff.value(), 1.2708966768404935, places=3)


def test_multiple_distance_restraints(verbose=False):
    """Check that multiple distance restraints give expected increase in energy"""

    mdr_on_nrg = getPotEnergyRestr(mult_dist_on=True, verbose=verbose)
    mdr_off_nrg = getPotEnergyRestr(mult_dist_on=False, verbose=verbose)
    nrg_diff = mdr_on_nrg - mdr_off_nrg

    if verbose:
        print("#######################################")
        print("Testing multiple distance restraints...")
        print(f"Energy of system with multiple distance restraints on = {mdr_on_nrg}")
        print(f"Energy of system with multiple distance restraints off = {mdr_off_nrg}")
        print(f"Energy difference when restraints turned on = {nrg_diff}")
        print("#######################################")

    assert_almost_equal(nrg_diff.value(), 5.6443886268, places=3)


if __name__ == '__main__':
    test_boresch_restraints(True)
    test_multiple_distance_restraints(True)
