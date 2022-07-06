
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

import Sire

from Sire.IO import *
from Sire.Move import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Units import *
from Sire.Maths import *
import os.path
import numpy as np

from nose.tools import assert_almost_equal
plugins = os.path.join('lib','plugins')
try:
    os.environ["OPENMM_PLUGIN_DIR"]
except KeyError:
    os.environ["OPENMM_PLUGIN_DIR"] = os.path.join(Sire.Base.getInstallDir(),plugins)


(mols, space) = Amber().readCrdTop("../io/ala.crd", "../io/ala.top")

coul_cutoff = 10 * angstrom
lj_cutoff = 10 * angstrom

rf_diel = 78.3

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff)

interff = InterCLJFF("interclj")
interff.setProperty("space", space)
interff.setProperty("switchingFunction", switchfunc)
interff.setUseReactionField(True)
interff.setReactionFieldDielectric(rf_diel)
interff.add(mols)

intraclj = IntraCLJFF("intraclj")
intraclj.setProperty("space", space)
intraclj.setProperty("switchingFunction", switchfunc )
intraclj.setUseReactionField(True)
intraclj.setReactionFieldDielectric(rf_diel)
intraclj.add(mols)

intraff = InternalFF("intraff")
intraff.add(mols)

system = System()
system.add(interff)
system.add(intraff)
system.add(intraclj)
system.add(mols)

system.setProperty("space", space)


def _pvt_create_npt_OpenMM(mols, temperature, pressure):
    openmm = OpenMMMDIntegrator(mols)
    openmm.setPlatform("CPU")
    openmm.setConstraintType("none")
    openmm.setCutoffType("cutoffperiodic")
    openmm.setIntegrator("leapfrogverlet")
    openmm.setFriction(0.1 * picosecond)
    openmm.setPrecision("double")
    openmm.setTimetoSkip(0*picosecond)
    openmm.setDeviceIndex("0")
    openmm.setLJDispersion(False)
    openmm.setFieldDielectric(rf_diel)
    openmm.setCMMremovalFrequency(0)
    openmm.setBufferFrequency(0)
    openmm.setRestraint(False)
    openmm.setTemperature(temperature)
    openmm.setAndersen(True)
    openmm.setAndersenFrequency(10)
    openmm.setPressure(pressure)
    openmm.setMCBarostat(True)
    openmm.setMCBarostatFrequency(25)
    openmm.initialise()
    return openmm

def _pvt_create_nve_OpenMM(mols):
    openmm = OpenMMMDIntegrator(mols)
    openmm.setPlatform("CPU")
    openmm.setConstraintType("none")
    openmm.setCutoffType("cutoffperiodic")
    openmm.setIntegrator("leapfrogverlet")
    openmm.setFriction(0.1 * picosecond)
    openmm.setPrecision("double")
    openmm.setTimetoSkip(0*picosecond)
    openmm.setDeviceIndex("0")
    openmm.setLJDispersion(False)
    openmm.setFieldDielectric(rf_diel)
    openmm.setCMMremovalFrequency(0)
    openmm.setBufferFrequency(0)
    openmm.setRestraint(False)
    openmm.setAndersen(False)
    openmm.setAndersenFrequency(10)
    openmm.setMCBarostat(False)
    openmm.setMCBarostatFrequency(25)
    openmm.initialise()
    return openmm

def _pvt_create_nvt_OpenMM(mols, temperature):
    openmm = OpenMMMDIntegrator(mols)
    openmm.setPlatform("CPU")
    openmm.setConstraintType("none")
    openmm.setCutoffType("cutoffperiodic")
    openmm.setIntegrator("leapfrogverlet")
    openmm.setFriction(0.1 * picosecond)
    openmm.setPrecision("double")
    openmm.setTimetoSkip(0*picosecond)
    openmm.setDeviceIndex("0")
    openmm.setLJDispersion(False)
    openmm.setFieldDielectric(rf_diel)
    openmm.setCMMremovalFrequency(0)
    openmm.setBufferFrequency(0)
    openmm.setRestraint(False)
    openmm.setTemperature(temperature)
    openmm.setAndersen(True)
    openmm.setAndersenFrequency(10)
    openmm.setMCBarostat(False)
    openmm.setMCBarostatFrequency(25)
    print("Initialising NVT integrator")
    openmm.initialise()
    return openmm

def test_getters_setters(verbose=False):
    openmm = _pvt_create_nve_OpenMM(system[mols.number()])
    assert (openmm.getAndersen()==False)
    assert (openmm.getPlatform()=="CPU")
    assert (openmm.getConstraintType()=="none")
    assert (openmm.getCutoffType()=="cutoffperiodic")
    assert (openmm.getIntegrator() == "leapfrogverlet")
    assert (openmm.getPrecision() == "double")
    assert (openmm.getLJDispersion() == False)


def test_npt_setup(verbose=False):
    if verbose:
        print ("=========NPT energy test============")

    sire_nrg = system.energy().value()

    if verbose:
        print("\nInitial Sire energy = %s kcal mol-1" % sire_nrg)

    openmm = _pvt_create_npt_OpenMM(system[mols.number()], 25*celsius, 1*atm)
    openmm_nrg = openmm.getPotentialEnergy(system).value()

    if verbose:
        print("\nInitial OpenMM energy = %s kcal mol-1" % openmm_nrg)

    assert_almost_equal(sire_nrg, openmm_nrg, 1)

    if verbose:
        print ("========NPT energy test done========")

def test_nvt_setup(verbose=False):
    if verbose:
        print ("=========NVT energy test============")

    sire_nrg = system.energy().value()

    if verbose:
        print("\nInitial Sire energy = %s kcal mol-1" % sire_nrg)

    openmm = _pvt_create_nvt_OpenMM(system[mols.number()], 25*celsius)

    openmm_nrg = openmm.getPotentialEnergy(system).value()

    if verbose:
        print("\nInitial OpenMM energy = %s kcal mol-1" % openmm_nrg)

    assert_almost_equal(sire_nrg, openmm_nrg, 1)

    if verbose:
        print ("========NVT energy test done========")

def test_nve_setup(verbose=False):
    if verbose:
        print ("=========NVE energy test============")

    sire_nrg = system.energy().value()

    if verbose:
        print("\nInitial Sire energy = %s kcal mol-1" % sire_nrg)

    openmm = _pvt_create_nve_OpenMM(system[mols.number()])

    openmm_nrg = openmm.getPotentialEnergy(system).value()

    if verbose:
        print("\nInitial OpenMM energy = %s kcal mol-1" % openmm_nrg)

    assert_almost_equal(sire_nrg, openmm_nrg, 1)

    if verbose:
        print ("========NVE energy test done========")

def _pvt_test_nve(verbose=False):

    openmm = _pvt_create_nve_OpenMM(system[mols.number()])
    openmm_nrg = openmm.getPotentialEnergy(system).value()

    sire_nrg = system.energy().value()
    if verbose:
        print("\nInitial Sire energy = %s kcal mol-1" % sire_nrg)
    if verbose:
        print("\nInitial OpenMM energy = %s kcal mol-1" % openmm_nrg)
    # now let's integrate a few steps
    #system_int = openmm.integrate(system, system.totalComponent(), 2*femtosecond, 10, True)
    #if verbose:
        #print ("Sire Energy after 10 steps of dynamics = %s kcal mol-1" %system_int.energ().value())

    #MD move integration
    velocity_generator = MaxwellBoltzmann(25*celsius)
    velocity_generator.setGenerator(RanGenerator(1234))

    mdmove = MolecularDynamics(mols, openmm, 1*femtosecond,
                              {"velocity generator":velocity_generator})
    tot_e = []
    for i in range(10):
        print("++++++++++++++++++++++++%d++++++++++++++++++++"%i)
        mdmove.move(system, 1, True)
        total_energy = system.energy()+mdmove.kineticEnergy()
        print ("Potential energy is %s and kinetic energy is %s" %(system.energy(),mdmove.kineticEnergy()))
        print ("Total energy from sire is: %s" %total_energy)
        tot_e.append(total_energy.value())
    percent_fluc = np.std(tot_e)/np.abs(np.mean(tot_e))*100
    print ("Percent fluctuations are: %s "%(percent_fluc))
    assert(percent_fluc<0.1)
        #print("Total energy is %s " %total_energy)




if __name__ == "__main__":
    test_npt_setup(True)
    test_nvt_setup(True)
    test_nve_setup(True)
    test_getters_setters(True)
    #test_nve(True)
