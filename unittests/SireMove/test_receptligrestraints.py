####################################################################################################
#                                                                                                  #
#          Sire script to test that receptor-ligand restraints are working correctly               #
#                                                                                                  #
####################################################################################################

####################################################################################################
#
#   IMPORTS
#
####################################################################################################

import os
import re
import sys

from Sire.Base import *

# Make sure that the OPENMM_PLUGIN_DIR enviroment variable is set correctly if unset.
try:
    # The user has already set the plugin location.
    os.environ["OPENMM_PLUGIN_DIR"]
except KeyError:
    # Set to the default location of the bundled OpenMM package.
    os.environ["OPENMM_PLUGIN_DIR"] = getLibDir() + "/plugins"

from nose.tools import assert_almost_equal
from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *
from Sire.Tools.DCDFile import *
from Sire.Tools import Parameter, resolveParameters, readParams
import Sire.Stream
import time
import numpy as np


MIN_MASSES = {'C': 5.96, 'N': 7.96}
HMR_MIN = 1.0
HMR_MAX = 4.0


####################################################################################################
#
#   Config file parameters
#
####################################################################################################
gpu = Parameter("gpu", 0, """The device ID of the GPU on which to run the simulation.""")

rf_dielectric = Parameter("reaction field dielectric", 78.3,
                          """Dielectric constant to use if the reaction field cutoff method is used.""")

temperature = Parameter("temperature", 25 * celsius, """Simulation temperature""")

pressure = Parameter("pressure", 1 * atm, """Simulation pressure""")

topfile = Parameter("topfile", "../io/receptor_ligand_restraints/SYSTEM.top",
                    """File name of the topology file containing the system to be simulated.""")

crdfile = Parameter("crdfile", "../io/receptor_ligand_restraints/SYSTEM.crd",
                    """File name of the coordinate file containing the coordinates of the
                       system to be simulated.""")

s3file = Parameter("s3file", "SYSTEM.s3",
                   """Filename for the system state file. The system state after topology and and coordinates
                   were loaded are saved in this file.""")

restart_file = Parameter("restart file", "sim_restart.s3",
                         """Filename of the restart file to use to save progress during the simulation.""")

dcd_root = Parameter("dcd root", "traj", """Root of the filename of the output DCD trajectory files.""")

nmoves = Parameter("nmoves", 1000, """Number of Molecular Dynamics moves to be performed during the simulation.""")

debug_seed = Parameter("debug seed", 0, """Debugging seed number seed. Set this if you
                         want to reproduce a single cycle. Don't use this seed for production simulations
                         since the same seed will be used for all cycles! A value of zero means that a unique
                         seed will be generated for each cycle.""")

ncycles = Parameter("ncycles", 1,
                    """The number of MD cycles. The total elapsed time will be nmoves*ncycles*timestep""")

maxcycles = Parameter("maxcycles",99999,
                      """The maximum number of MD cycles to carry out. Useful to restart simulations from a checkpoint""")

ncycles_per_snap = Parameter("ncycles_per_snap", 1, """Number of cycles between saving snapshots""")

save_coords = Parameter("save coordinates", True, """Whether or not to save coordinates.""")

buffered_coords_freq = Parameter("buffered coordinates frequency", 1,
                                 """The number of time steps between saving of coordinates during
                                 a cycle of MD. 0 disables buffering.""")
minimal_coordinate_saving = Parameter("minimal coordinate saving", False, "Reduce the number of coordinates writing for states"
                                                                    "with lambda in ]0,1[")

time_to_skip = Parameter("time to skip", 0 * picosecond, """Time to skip in picoseconds""")

minimise = Parameter("minimise", False, """Whether or not to perform minimization before the simulation.""")

minimise_tol = Parameter("minimise tolerance", 1, """Tolerance used to know when minimization is complete.""")

minimise_max_iter = Parameter("minimise maximum iterations", 1000, """Maximum number of iterations for minimization.""")

equilibrate = Parameter("equilibrate", False , """Whether or not to perform equilibration before dynamics.""")

equil_iterations = Parameter("equilibration iterations", 2000, """Number of equilibration steps to perform.""")

equil_timestep = Parameter("equilibration timestep", 0.5 * femtosecond, """Timestep to use during equilibration.""")

combining_rules = Parameter("combining rules", "arithmetic",
                            """Combining rules to use for the non-bonded interactions.""")

timestep = Parameter("timestep", 2 * femtosecond, """Timestep for the dynamics simulation.""")

platform = Parameter("platform", "CUDA", """Which OpenMM platform should be used to perform the dynamics.""")

precision = Parameter("precision", "mixed", """The floating point precision model to use during dynamics.""")

constraint = Parameter("constraint", "hbonds", """The constraint model to use during dynamics.""")

cutoff_type = Parameter("cutoff type", "cutoffperiodic", """The cutoff method to use during the simulation.""")

cutoff_dist = Parameter("cutoff distance", 10 * angstrom,
                        """The cutoff distance to use for the non-bonded interactions.""")

integrator_type = Parameter("integrator", "leapfrogverlet", """The integrator to use for dynamics.""")

inverse_friction = Parameter("inverse friction", 0.1 * picosecond,
                             """Inverse friction time for the Langevin thermostat.""")

andersen = Parameter("thermostat", True,
                     """Whether or not to use the Andersen thermostat (needed for NVT or NPT simulation).""")

barostat = Parameter("barostat", True, """Whether or not to use a barostat (needed for NPT simulation).""")

andersen_frequency = Parameter("andersen frequency", 10.0, """Collision frequency in units of (1/ps)""")

barostat_frequency = Parameter("barostat frequency", 25,
                               """Number of steps before attempting box changes if using the barostat.""")

lj_dispersion = Parameter("lj dispersion", False, """Whether or not to calculate and include the LJ dispersion term.""")

cmm_removal = Parameter("center of mass frequency", 10,
                        """Frequency of which the system center of mass motion is removed.""")

center_solute = Parameter("center solute", False,
                          """Whether or not to centre the centre of geometry of the solute in the box.""")

use_restraints = Parameter("use restraints", False, """Whether or not to use harmonic restaints on the solute atoms.""")

k_restraint = Parameter("restraint force constant", 100.0, """Force constant to use for the harmonic restaints.""")

heavy_mass_restraint = Parameter("heavy mass restraint", 1.10,
                                 """Only restrain solute atoms whose mass is greater than this value.""")

unrestrained_residues = Parameter("unrestrained residues", ["WAT", "HOH"],
                                  """Names of residues that are never restrained.""")

freeze_residues = Parameter("freeze residues", False, """Whether or not to freeze certain residues.""")

frozen_residues = Parameter("frozen residues", ["LGR", "SIT", "NEG", "POS"],
                            """List of residues to freeze if 'freeze residues' is True.""")


use_distance_restraints = Parameter("use distance restraints", False,
                                    """Whether or not to use restraints distances between pairs of atoms.""")

distance_restraints_dict = Parameter("distance restraints dictionary", {(21, 4950): (2.72, 19.55, 0), (11, 4946): (5.92, 
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

turn_on_restraints_mode = Parameter("turn on receptor-ligand restraints mode", False,
                                  """If true, lambda will be used to scale the receptor-ligand restraint strength. A dummy
                                  pert file mapping all original ligand atom parameters to themselves must be supplied.""")

use_boresch_restraints = Parameter("use boresch restraints", False, 
                                    """Whether or not to use Boresch restraints between the ligand and receptor""")

boresch_restraints_dict = Parameter("boresch restraints dictionary", {"anchor_points":{"r1":4946, "r2":4944, "r3":4949,
                                     "l1":11, "l2":2, "l3":3},"equilibrium_values":{"r0":5.92, "thetaA0":1.85, 
                                     "thetaB0":1.59,"phiA0":-0.30, "phiB0":-1.55, "phiC0":2.90},
                                     "force_constants":{"kr":25.49, "kthetaA":66.74, "kthetaB":38.39, "kphiA":215.36,
                                      "kphiB":49.23, "kphiC":49.79}},
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

## Free energy specific keywords
morphfile = Parameter("morphfile", "../io/receptor_ligand_restraints/MORPH.dummy.pert",
                      """Name of the morph file containing the perturbation to apply to the system.""")

lambda_val = Parameter("lambda_val", 0.0,
                       """Value of the lambda parameter at which to evaluate free energy gradients.""")

delta_lambda = Parameter("delta_lambda", 0.001,
                         """Value of the lambda interval used to evaluate free energy gradients by finite difference.""")

lambda_array = Parameter("lambda array",[] ,
                        """Array with all lambda values lambda_val needs to be part of the array. """)

shift_delta = Parameter("shift delta", 2.0,
                        """Value of the Lennard-Jones soft-core parameter.""")

coulomb_power = Parameter("coulomb power", 0,
                          """Value of the Coulombic soft-core parameter.""")

energy_frequency = Parameter("energy frequency", 1,
                             """The number of time steps between evaluation of free energy gradients.""")

simfile = Parameter("outdata_file", "simfile.dat", """Filename that records all output needed for the free energy analysis""")

perturbed_resnum = Parameter("perturbed residue number",1,"""The residue number of the molecule to morph.""")


####################################################################################################
#
#   Helper functions
#
####################################################################################################

def getSolute(system):
    """Find the solute molecule based on the perturbed residue number.

    Args:
        system (system): The Sire system

    Returns:
        molecule: Molecule matching perturbed residue number assumed to be solvent
    """

    # Search the system for a single molcule containing a residue
    # matching the perturbed_resnum.val.

    # Create the query string.
    query = f"mol with resnum {perturbed_resnum.val}"

    # Perform the search.
    search = system.search(query)

    # Make sure there is only one result.
    if len(search) != 1:
        msg = ("FATAL! Could not find a solute to perturb with residue "
              f"number {perturbed_resnum.val} in the input! Check the value of "
               "your config keyword 'perturbed residue number' The system should "
               "contain a single molecule with this residue number.")
        raise Exception(msg)

    # Return the matching molecule, i.e. the solute.
    return search[0]

def atomNumListToProperty(list):

    prop = Properties()
    i = 0
    for value in list:
        prop.setProperty(str(i), VariantProperty(value.value()))
        i += 1
    return prop


def atomNumVectorListToProperty(list):
    prop = Properties()

    i = 0

    for value in list:
        prop.setProperty("AtomNum(%d)" % i, VariantProperty(value[0].value()))
        prop.setProperty("x(%d)" % i, VariantProperty(value[1].x()))
        prop.setProperty("y(%d)" % i, VariantProperty(value[1].y()))
        prop.setProperty("z(%d)" % i, VariantProperty(value[1].z()))
        prop.setProperty("k(%d)" % i, VariantProperty(value[2].val ) )
        i += 1

    prop.setProperty("nrestrainedatoms", VariantProperty(i));

    return prop


def linkbondVectorListToProperty(list):

    prop = Properties()

    i = 0

    for value in list:
        prop.setProperty("AtomNum0(%d)" % i, VariantProperty(value[0]))
        prop.setProperty("AtomNum1(%d)" % i, VariantProperty(value[1]))
        prop.setProperty("reql(%d)" % i, VariantProperty(value[2]))
        prop.setProperty("kl(%d)" % i, VariantProperty(value[3]))
        prop.setProperty("dl(%d)" % i, VariantProperty(value[4]))
        i += 1

    prop.setProperty("nbondlinks", VariantProperty(i));

    return prop


def boreschDistRestraintToProperty(boresch_dict):
    """Generates properties to store information needed to set up the single
    Boresch distance restraint.

    Args:
        boresch_dict (dict): Containts the information required to set up all
        Boresch restraints

    Returns:
        class 'Sire.Base._Base.Properties': The properties required to
        set up the Boresch distance restraint
    """

    prop = Properties()

    prop.setProperty("AtomNum0", VariantProperty(boresch_dict['anchor_points']['l1']))
    prop.setProperty("AtomNum1", VariantProperty(boresch_dict['anchor_points']['r1']))
    prop.setProperty("equil_val", VariantProperty(boresch_dict['equilibrium_values']['r0']))
    prop.setProperty("force_const", VariantProperty(boresch_dict['force_constants']['kr']))

    return prop


def boreschAngleRestraintsToProperty(boresch_dict): 
    """Generates properties to store information needed to set up the two
    Boresch angle restraints.

    Args:
        boresch_dict (dict): Containts the information required to set up all
        Boresch restraints

    Returns:
        class 'Sire.Base._Base.Properties': The properties required to
        set up the Boresch angle restraints 
    """

    prop = Properties()

    angle_anchor_dict = {"thetaA":["r2", "r1", "l1"], "thetaB":["r1", "l1", "l2"]}

    i = 0
    for angle in ["thetaA", "thetaB"]:
        if boresch_dict["force_constants"][f"k{angle}"] != 0:
            for j in range(3): 
                prop.setProperty(f"AtomNum{j}-{i}", 
                    VariantProperty(boresch_dict['anchor_points'][angle_anchor_dict[angle][j]]))
            prop.setProperty(f"equil_val-{i}", 
                VariantProperty(boresch_dict["equilibrium_values"][f"{angle}0"]))
            prop.setProperty(f"force_const-{i}", 
                VariantProperty(boresch_dict["force_constants"][f"k{angle}"]))
            
            i += 1
            
    prop.setProperty("n_boresch_angle_restraints", VariantProperty(i));

    return prop


def boreschDihedralRestraintsToProperty(boresch_dict): 
    """Generates properties to store information needed to set up the three
    Boresch dihedral restraints.

    Args:
        boresch_dict (dict): Containts the information required to set up all
        Boresch restraints

    Returns:
        class 'Sire.Base._Base.Properties': The properties required to
        set up the Boresch dihedral restraints 
    """

    prop = Properties()

    dihedral_anchor_dict = {"phiA":["r3", "r2", "r1", "l1"], "phiB":["r2", "r1", "l1", "l2"],
                            "phiC":["r1", "l1", "l2", "l3"]}

    i = 0
    for dihedral in ["phiA", "phiB", "phiC"]:
        if boresch_dict["force_constants"][f"k{dihedral}"] != 0:
            for j in range(4): 
                prop.setProperty(f"AtomNum{j}-{i}", 
                    VariantProperty(boresch_dict['anchor_points'][dihedral_anchor_dict[dihedral][j]]))
            prop.setProperty(f"equil_val-{i}", 
                VariantProperty(boresch_dict["equilibrium_values"][f"{dihedral}0"]))
            prop.setProperty(f"force_const-{i}", 
                VariantProperty(boresch_dict["force_constants"][f"k{dihedral}"]))
            
            i += 1
            
    prop.setProperty("n_boresch_dihedral_restraints", VariantProperty(i));

    return prop


def propertyToAtomNumList(prop):
    list = []
    i = 0
    try:
        while True:
            list.append(AtomNum(prop[str(i)].toInt()))
            i += 1
    except:
        pass
    return list

def propertyToAtomNumVectorList(prop):
    list = []
    i = 0
    try:
        while True:
            num = AtomNum(prop["AtomNum(%d)" % i].toInt())
            x = prop["x(%d)" % i].toDouble()
            y = prop["y(%d)" % i].toDouble()
            z = prop["z(%d)" % i].toDouble()
            k = prop["k(%d)" % i].toDouble()

            list.append((num, Vector(x, y, z), k ))

            i += 1
    except:
        pass

    return list


def setupDistanceRestraints(system, restraints=None):
    prop_list = []

    molecules = system[MGName("all")].molecules()

    if restraints is None:
        #dic_items = list(distance_restraints_dict.val.items())
        dic_items = list(dict(distance_restraints_dict.val).items())
    else:
        dic_items = list(restraints.items())

    molecules = system[MGName("all")].molecules()
    moleculeNumbers = molecules.molNums()

    for moleculeNumber in moleculeNumbers:
        mol = molecules.molecule(moleculeNumber)[0].molecule()
        atoms_mol = mol.atoms()
        natoms_mol = mol.nAtoms()
        for j in range(0, natoms_mol):
            at = atoms_mol[j]
            atnumber = at.number()
            for k in range(len(dic_items)):
                if dic_items[k][0][0] == dic_items[k][0][1]:
                    print ("Error! It is not possible to place a distance restraint on the same atom")
                    sys.exit(-1)
                if atnumber.value() - 1 in dic_items[k][0]:
                    # atom0index atom1index, reql, kl, dl
                    prop_list.append((
                        dic_items[k][0][0] + 1, dic_items[k][0][1] + 1, dic_items[k][1][0], dic_items[k][1][1],
                        dic_items[k][1][2]))

    unique_prop_list = []

    [unique_prop_list.append(item) for item in prop_list if item not in unique_prop_list]
    # The solute will store all the information related to the receptor-ligand restraints
    solute = getSolute(system)
    solute = solute.edit().setProperty("linkbonds", linkbondVectorListToProperty(unique_prop_list)).commit()
    system.update(solute)

    return system


def setupBoreschRestraints(system):
    """Takes initial system and adds information specifying the Boresch
    restraints. The distance, angle, and torsional restraints are stored as
    properties in solute molecule.

    Args:
        system (System): The initial system

    Returns:
        System: The updated system with
        Boresch restraint properties stored in the solute.
    """
    # Get Boresch restraint dict in dict form
    boresch_dict = dict(boresch_restraints_dict.val)

    # Correct atom numbers by + 1
    for key in boresch_dict["anchor_points"].keys():
        boresch_dict["anchor_points"][key] += 1

    # Get anchor points dicts
    anchors_dict = boresch_dict["anchor_points"]
    
    molecules = system[MGName("all")].molecules()
    moleculeNumbers = molecules.molNums()
    
    # The solute will store all the information related to the Boresch restraints in the system
    solute = getSolute(system)
    solute = solute.edit().setProperty("boresch_dist_restraint", boreschDistRestraintToProperty(boresch_dict)).commit()
    solute = solute.edit().setProperty("boresch_angle_restraints", boreschAngleRestraintsToProperty(boresch_dict)).commit()
    solute = solute.edit().setProperty("boresch_dihedral_restraints", boreschDihedralRestraintsToProperty(boresch_dict)).commit()
    system.update(solute)

    return system

def getDummies(molecule):
    natoms = molecule.nAtoms()
    atoms = molecule.atoms()

    from_dummies = None
    to_dummies = None

    for x in range(0, natoms):
        atom = atoms[x]
        if atom.property("initial_ambertype") == "du":
            if from_dummies is None:
                from_dummies = molecule.selectAll(atom.index())
            else:
                from_dummies += molecule.selectAll(atom.index())
        elif atom.property("final_ambertype") == "du":
            if to_dummies is None:
                to_dummies = molecule.selectAll(atom.index())
            else:
                to_dummies += molecule.selectAll(atom.index())

    return to_dummies, from_dummies


def createSystemFreeEnergy(molecules):
    r"""creates the system for free energy calculation
    Parameters
    ----------
    molecules : Sire.molecules
        Sire object that contains a lot of information about molecules
    Returns
    -------
    system : Sire.system

    """
    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber)[0].molecule()
        moleculeList.append(molecule)

    # Scan input to find a molecule with passed residue number
    # The residue name of the first residue in this molecule is
    # used to name the solute. This is used later to match
    # templates in the flex/pert files.

    solute = None
    for molecule in moleculeList:
        if ( molecule.residue(ResIdx(0)).number() == ResNum(perturbed_resnum.val) ):
            solute = molecule
            moleculeList.remove(molecule)
            break

    if solute is None:
        msg = ("FATAL! Could not find a solute to perturb with residue "
              f"number {perturbed_resnum.val} in the input! Check the value of "
               "your config keyword 'perturbed residue number' The system should "
               "contain a single molecule with this residue number.")
        raise Exception(msg)

    #solute = moleculeList[0]

    lig_name = solute.residue(ResIdx(0)).name().value()

    solute = solute.edit().rename(lig_name).commit()

    perturbations_lib = PerturbationsLibrary(morphfile.val)
    solute = perturbations_lib.applyTemplate(solute)

    perturbations = solute.property("perturbations")

    lam = Symbol("lambda")

    initial = Perturbation.symbols().initial()
    final = Perturbation.symbols().final()

    solute = solute.edit().setProperty("perturbations",
                                       perturbations.recreate((1 - lam) * initial + lam * final)).commit()

    # We put atoms in three groups depending on what happens in the perturbation
    # non dummy to non dummy --> the hard group, uses a normal intermolecular FF
    # non dummy to dummy --> the todummy group, uses SoftFF with alpha = Lambda
    # dummy to non dummy --> the fromdummy group, uses SoftFF with alpha = 1 - Lambda
    # We start assuming all atoms are hard atoms. Then we call getDummies to find which atoms
    # start/end as dummies and update the hard, todummy and fromdummy groups accordingly

    solute_grp_ref = MoleculeGroup("solute_ref", solute)
    solute_grp_ref_hard = MoleculeGroup("solute_ref_hard")
    solute_grp_ref_todummy = MoleculeGroup("solute_ref_todummy")
    solute_grp_ref_fromdummy = MoleculeGroup("solute_ref_fromdummy")

    solute_ref_hard = solute.selectAllAtoms()
    solute_ref_todummy = solute_ref_hard.invert()
    solute_ref_fromdummy = solute_ref_hard.invert()

    to_dummies, from_dummies = getDummies(solute)

    if to_dummies is not None:
        ndummies = to_dummies.count()
        dummies = to_dummies.atoms()

        for x in range(0, ndummies):
            dummy_index = dummies[x].index()
            solute_ref_hard = solute_ref_hard.subtract(solute.select(dummy_index))
            solute_ref_todummy = solute_ref_todummy.add(solute.select(dummy_index))

    if from_dummies is not None:
        ndummies = from_dummies.count()
        dummies = from_dummies.atoms()

        for x in range(0, ndummies):
            dummy_index = dummies[x].index()
            solute_ref_hard = solute_ref_hard.subtract(solute.select(dummy_index))
            solute_ref_fromdummy = solute_ref_fromdummy.add(solute.select(dummy_index))

    solute_grp_ref_hard.add(solute_ref_hard)
    solute_grp_ref_todummy.add(solute_ref_todummy)
    solute_grp_ref_fromdummy.add(solute_ref_fromdummy)

    solutes = MoleculeGroup("solutes")
    solutes.add(solute)

    molecules = MoleculeGroup("molecules")
    molecules.add(solute)

    solvent = MoleculeGroup("solvent")

    #for molecule in moleculeList[1:]:
    for molecule in moleculeList:
        molecules.add(molecule)
        solvent.add(molecule)

    all = MoleculeGroup("all")

    all.add(molecules)
    all.add(solvent)

    all.add(solutes)
    all.add(solute_grp_ref)
    all.add(solute_grp_ref_hard)
    all.add(solute_grp_ref_todummy)
    all.add(solute_grp_ref_fromdummy)

    # Add these groups to the System
    system = System()

    system.add(solutes)
    system.add(solute_grp_ref)
    system.add(solute_grp_ref_hard)
    system.add(solute_grp_ref_todummy)
    system.add(solute_grp_ref_fromdummy)

    system.add(molecules)

    system.add(solvent)

    system.add(all)

    return system


def setupForceFieldsFreeEnergy(system, space):
    r"""sets up the force field for the free energy calculation
    Parameters
    ----------
    system : Sire.system
    space : Sire.space
    Returns
    -------
    system : Sire.system
    """
    solutes = system[MGName("solutes")]

    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]

    solvent = system[MGName("solvent")]

    all = system[MGName("all")]

    # ''solvent'' is actually every molecule that isn't perturbed !
    solvent_intraff = InternalFF("solvent_intraff")
    solvent_intraff.add(solvent)

    # Solute bond, angle, dihedral energy
    solute_intraff = InternalFF("solute_intraff")
    solute_intraff.add(solute)

    # Solvent-solvent coulomb/LJ (CLJ) energy
    solventff = InterCLJFF("solvent:solvent")
    if (cutoff_type.val != "nocutoff"):
        solventff.setUseReactionField(True)
        solventff.setReactionFieldDielectric(rf_dielectric.val)
    solventff.add(solvent)

    #Solvent intramolecular CLJ energy
    solvent_intraclj = IntraCLJFF("solvent_intraclj")
    if (cutoff_type.val != "nocutoff"):
        solvent_intraclj.setUseReactionField(True)
        solvent_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solvent_intraclj.add(solvent)

    # Solute intramolecular CLJ energy
    solute_hard_intraclj = IntraCLJFF("solute_hard_intraclj")
    if (cutoff_type.val != "nocutoff"):
        solute_hard_intraclj.setUseReactionField(True)
        solute_hard_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_intraclj.add(solute_hard)

    solute_todummy_intraclj = IntraSoftCLJFF("solute_todummy_intraclj")
    solute_todummy_intraclj.setShiftDelta(shift_delta.val)
    solute_todummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_intraclj.setUseReactionField(True)
        solute_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_intraclj.add(solute_todummy)

    solute_fromdummy_intraclj = IntraSoftCLJFF("solute_fromdummy_intraclj")
    solute_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_intraclj.setUseReactionField(True)
        solute_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_intraclj.add(solute_fromdummy)

    solute_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_hard:todummy_intraclj")
    solute_hard_todummy_intraclj.setShiftDelta(shift_delta.val)
    solute_hard_todummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_hard_todummy_intraclj.setUseReactionField(True)
        solute_hard_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_todummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intraclj.add(solute_todummy, MGIdx(1))

    solute_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_hard:fromdummy_intraclj")
    solute_hard_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_hard_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_hard_fromdummy_intraclj.setUseReactionField(True)
        solute_hard_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_fromdummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))

    solute_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intraclj")
    solute_todummy_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_todummy_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_fromdummy_intraclj.setUseReactionField(True)
        solute_todummy_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_fromdummy_intraclj.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))

    #Solute-solvent CLJ energy
    solute_hard_solventff = InterGroupCLJFF("solute_hard:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_hard_solventff.setUseReactionField(True)
        solute_hard_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_solventff.add(solute_hard, MGIdx(0))
    solute_hard_solventff.add(solvent, MGIdx(1))

    solute_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_solventff.setUseReactionField(True)
        solute_todummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_todummy_solventff.add(solvent, MGIdx(1))

    solute_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_solventff.setUseReactionField(True)
        solute_fromdummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_solventff.add(solvent, MGIdx(1))


    # TOTAL
    forcefields = [solute_intraff,
                   solute_hard_intraclj, solute_todummy_intraclj, solute_fromdummy_intraclj,
                   solute_hard_todummy_intraclj, solute_hard_fromdummy_intraclj,
                   solute_todummy_fromdummy_intraclj,
                   solvent_intraff,
                   solventff, solvent_intraclj,
                   solute_hard_solventff, solute_todummy_solventff, solute_fromdummy_solventff]


    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)

    if (cutoff_type.val != "nocutoff"):
        system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff_dist.val))
    else:
        system.setProperty("switchingFunction", NoCutoff())

    system.setProperty("combiningRules", VariantProperty(combining_rules.val))
    system.setProperty("coulombPower", VariantProperty(coulomb_power.val))
    system.setProperty("shiftDelta", VariantProperty(shift_delta.val))

    # TOTAL
    total_nrg = solute_intraff.components().total() + solute_hard_intraclj.components().total() + \
                solute_todummy_intraclj.components().total(0) + solute_fromdummy_intraclj.components().total(0) + \
                solute_hard_todummy_intraclj.components().total(
                    0) + solute_hard_fromdummy_intraclj.components().total(0) + \
                solute_todummy_fromdummy_intraclj.components().total(0) + \
                solvent_intraff.components().total() + solventff.components().total() + \
                solvent_intraclj.components().total() + \
                solute_hard_solventff.components().total() + \
                solute_todummy_solventff.components().total(0) + \
                solute_fromdummy_solventff.components().total(0)

    e_total = system.totalComponent()

    lam = Symbol("lambda")

    system.setComponent(e_total, total_nrg)

    system.setConstant(lam, 0.0)

    system.add(PerturbationConstraint(solutes))

    # NON BONDED Alpha constraints for the soft force fields

    system.add(PropertyConstraint("alpha0", FFName("solute_todummy_intraclj"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy_intraclj"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:todummy_intraclj"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:fromdummy_intraclj"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:fromdummy_intraclj"), Max(lam, 1 - lam)))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:solvent"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy:solvent"), 1 - lam))

    system.setComponent(lam, lambda_val.val)

    # printEnergies( system.componentValues() )

    return system


def setupMovesFreeEnergy(system, debug_seed, GPUS, lam_val):

    molecules = system[MGName("molecules")]
    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]
    #import pdb ; pdb.set_trace()
    Integrator_OpenMM = OpenMMFrEnergyST(molecules, solute, solute_hard, solute_todummy, solute_fromdummy)
    Integrator_OpenMM.setRandomSeed(debug_seed)
    Integrator_OpenMM.setIntegrator(integrator_type.val)
    Integrator_OpenMM.setFriction(inverse_friction.val)  # Only meaningful for Langevin/Brownian integrators
    Integrator_OpenMM.setPlatform(platform.val)
    Integrator_OpenMM.setConstraintType(constraint.val)
    Integrator_OpenMM.setCutoffType(cutoff_type.val)
    Integrator_OpenMM.setFieldDielectric(rf_dielectric.val)
    Integrator_OpenMM.setAlchemicalValue(lambda_val.val)
    Integrator_OpenMM.setAlchemicalArray(lambda_array.val)
    Integrator_OpenMM.setDeviceIndex(str(GPUS))
    Integrator_OpenMM.setCoulombPower(coulomb_power.val)
    Integrator_OpenMM.setShiftDelta(shift_delta.val)
    Integrator_OpenMM.setDeltatAlchemical(delta_lambda.val)
    Integrator_OpenMM.setPrecision(precision.val)
    Integrator_OpenMM.setTimetoSkip(time_to_skip.val)
    Integrator_OpenMM.setBufferFrequency(buffered_coords_freq.val)

    if cutoff_type != "nocutoff":
        Integrator_OpenMM.setCutoffDistance(cutoff_dist.val)

    Integrator_OpenMM.setCMMremovalFrequency(cmm_removal.val)

    Integrator_OpenMM.setEnergyFrequency(energy_frequency.val)

    if use_restraints.val:
        Integrator_OpenMM.setRestraint(True)

    if andersen.val:
        Integrator_OpenMM.setTemperature(temperature.val)
        Integrator_OpenMM.setAndersen(andersen.val)
        Integrator_OpenMM.setAndersenFrequency(andersen_frequency.val)

    if barostat.val:
        Integrator_OpenMM.setPressure(pressure.val)
        Integrator_OpenMM.setMCBarostat(barostat.val)
        Integrator_OpenMM.setMCBarostatFrequency(barostat_frequency.val)

    seed = RanGenerator().randInt(100000, 1000000)

    #This calls the OpenMMFrEnergyST initialise function
    Integrator_OpenMM.initialise()
    velocity_generator = MaxwellBoltzmann(temperature.val)
    velocity_generator.setGenerator(RanGenerator(seed))

    mdmove = MolecularDynamics(molecules, Integrator_OpenMM, timestep.val,
                              {"velocity generator":velocity_generator})

    moves = WeightedMoves()
    moves.add(mdmove, 1)

    moves.setGenerator(RanGenerator(seed))

    return moves

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

    amber = Amber()
    (molecules, space) = amber.readCrdTop(crdfile.val, topfile.val)
    system = createSystemFreeEnergy(molecules)

    if mult_dist_on:
        system = setupDistanceRestraints(system, restraints=distance_restraints_dict.val)

    if boresch_on:
        system = setupBoreschRestraints(system)

    system = setupForceFieldsFreeEnergy(system, space)
    moves = setupMovesFreeEnergy(system, debug_seed.val, gpu.val, lambda_val.val)
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
        print("Testing Boresch restraints...")
        print(f"Energy of system with Boresch restraints on = {boresch_on_nrg}")
        print(f"Energy of system with Boresch restraints off = {boresch_off_nrg}")
        print(f"Energy difference = {nrg_diff}")

    assert_almost_equal(nrg_diff.value(), 1.27076849)


def test_multiple_distance_restraints(verbose=False):
    """Check that multiple distance restraints give expected increase in energy"""

    mdr_on_nrg = getPotEnergyRestr(mult_dist_on=True, verbose=verbose)
    mdr_off_nrg = getPotEnergyRestr(mult_dist_on=False, verbose=verbose)
    nrg_diff = mdr_on_nrg - mdr_off_nrg

    if verbose:
        print("Testing multiple distance restraints...")
        print(f"Energy of system with multiple distance restraints on = {mdr_on_nrg}")
        print(f"Energy of system with multiple distance restraints off = {mdr_off_nrg}")
        print(f"Energy difference = {nrg_diff}")

    assert_almost_equal(nrg_diff.value(), 5.6443886268)


if __name__ == '__main__':
    test_boresch_restraints(True)
    test_multiple_distance_restraints(True)
