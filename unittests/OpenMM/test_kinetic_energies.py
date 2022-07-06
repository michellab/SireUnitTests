
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

np = Sire.try_import("numpy")

from nose.tools import assert_almost_equal
plugins = os.path.join('lib','plugins')
try:
    os.environ["OPENMM_PLUGIN_DIR"]
except KeyError:
    os.environ["OPENMM_PLUGIN_DIR"] = os.path.join(Sire.Base.getInstallDir(),plugins)


(mols, space) = Amber().readCrdTop("../io/ethane.crd", "../io/ethane.top")

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

def _pvt_create_nve_OpenMM(mols):
    openmm = OpenMMMDIntegrator(mols)
    openmm.setPlatform("CPU")
    openmm.setConstraintType("none")
    openmm.setCutoffType("cutoffnonperiodic")
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

def test_kinetic_energy():
    openmm = _pvt_create_nve_OpenMM(system[mols.number()])
    #MD move integration
    velocity_generator = MaxwellBoltzmann(25*celsius)
    velocity_generator.setGenerator(RanGenerator(1234))

    mdmove = MolecularDynamics(mols, openmm, 0.001*picosecond,
                              {"velocity generator":velocity_generator})
    tot_e = []
    for i in range(10):
        print("++++++++++++++++++++++++%d++++++++++++++++++++"%i)
        mdmove.move(system, 1, True)
    #    total_energy = system.energy()+openmm.getKineticEnergy(system[mols.number()])
    #    print ("Potential energy is %s and kinetic energy is %s" %(system.energy(),mdmove.kineticEnergy()))
    #    print ("Total energy from sire is: %s" %total_energy)
    #    tot_e.append(total_energy.value())
    #percent_fluc = np.std(tot_e)/np.abs(np.mean(tot_e))*100
    #print ("Percent fluctuations are: %s "%(percent_fluc))
    #assert(percent_fluc<0.1)

if __name__ == "__main__":
    test_kinetic_energy()
