
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.MM import *
from Sire.FF import *
from Sire.System import *
from Sire.Move import *
from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Vol import *
from Sire.Base import *
from Sire.Units import *
from Sire.Qt import *

import Sire.Config
import Sire.Stream

from nose.tools import assert_almost_equal

import os

protein = Sire.Stream.load("../io/protein.s3")

zmat_maker = ZmatrixMaker()
zmat_maker.loadTemplates( os.path.join(Sire.Config.parameter_directory, "amber.zmatrices") )
protein = zmat_maker.applyTemplates(protein)

nmoves = 100
ranseed = 3989
temperature = 25 * celsius

space = Cartesian()

coul_cutoff = 15 * angstrom
lj_cutoff = 10 * angstrom

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff)

mols = MoleculeGroup("protein")

for i in range(0,protein.nResidues()):
    mols.add( protein.residue(ResIdx(i)) )

oldff = IntraCLJFF("oldff")
oldff.setSwitchingFunction(switchfunc)
oldff.setSpace(space)
oldff.setShiftElectrostatics(True)

oldff.add(protein)

newff = IntraFF("newff")
newff.setProperty("switchingFunction",switchfunc)
newff.setProperty("space",space)
newff.enableParallelCalculation()

newff.add(protein)

oldff14 = InternalFF("oldff14")
oldff14.disable14Calculation()
oldff14.add(protein)

newff14 = InternalFF("newff14")
newff14.enable14Calculation()
newff14.add(protein)

def test_energy(verbose = False):
    t = QElapsedTimer()

    newff.mustNowRecalculateFromScratch()
    oldff.mustNowRecalculateFromScratch()

    t.start()
    oldnrgs = oldff.energies()
    oldns = t.nsecsElapsed()
    oldcnrg = oldff.energy(oldff.components().coulomb()).value()
    oldljnrg = oldff.energy(oldff.components().lj()).value()

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()
    newcnrg = newff.energy(newff.components().coulomb()).value()
    newljnrg = newff.energy(newff.components().lj()).value()

    t.start()
    oldnrgs = oldff14.energies()
    old14ns = t.nsecsElapsed()

    oldcnrg += oldff14.energy( oldff14.components().intra14Coulomb() ).value()
    oldljnrg += oldff14.energy( oldff14.components().intra14LJ() ).value()

    t.start()
    newnrgs = newff14.energies()
    new14ns = t.nsecsElapsed()

    newcnrg += newff14.energy( newff14.components().intra14Coulomb() ).value()
    newljnrg += newff14.energy( newff14.components().intra14LJ() ).value()

    if verbose:
        print("\nTotal energy")
        print("OLD FF :  %s  %s  %s  : %s | %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns, 0.000001*old14ns))
        print("NEW FF :  %s  %s  %s  : %s | %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns, 0.000001*new14ns))

    assert_almost_equal( oldcnrg, newcnrg, 1 )
    assert_almost_equal( oldljnrg, newljnrg, 1 )

def test_sim(verbose = False):

    oldsys = System()
    newsys = System()

    oldsys.add(mols)
    newsys.add(mols)

    oldsys.add(oldff)
    oldsys.add(oldff14)
    newsys.add(newff)
    newsys.add(newff14)

    t = QElapsedTimer()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()

    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    oldcnrg = oldsys.energy( oldff.components().coulomb() ).value() + \
              oldsys.energy( oldff14.components().intra14Coulomb() ).value()
    oldljnrg = oldsys.energy( oldff.components().lj() ).value() + \
               oldsys.energy( oldff14.components().intra14LJ() ).value()

    newcnrg = newsys.energy( newff.components().coulomb() ).value() + \
              newsys.energy( newff14.components().intra14Coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value() + \
               newsys.energy( newff14.components().intra14LJ() ).value()

    if verbose:
        print("\nStarting energy")
        print("OLD SYS:  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))

    assert_almost_equal( oldcnrg, newcnrg, 1 )
    assert_almost_equal( oldljnrg, newljnrg, 1 )

    moves = ZMatMove(mols)
    moves.setGenerator( RanGenerator(ranseed) )
    moves.setTemperature(temperature)

    t.start()
    moves.move(oldsys, nmoves, False)
    move_oldns = t.nsecsElapsed()

    old_naccepted = moves.nAccepted()
    old_nrejected = moves.nRejected()

    moves.setGenerator( RanGenerator(ranseed) )
    moves.clearStatistics()

    t.start()
    moves.move(newsys, nmoves, False)
    move_newns = t.nsecsElapsed()

    new_naccepted = moves.nAccepted()
    new_nrejected = moves.nRejected()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()

    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    oldcnrg = oldsys.energy( oldff.components().coulomb() ).value() + \
              oldsys.energy( oldff14.components().intra14Coulomb() ).value()
    oldljnrg = oldsys.energy( oldff.components().lj() ).value() + \
               oldsys.energy( oldff14.components().intra14LJ() ).value()

    newcnrg = newsys.energy( newff.components().coulomb() ).value() + \
              newsys.energy( newff14.components().intra14Coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value() + \
               newsys.energy( newff14.components().intra14LJ() ).value()

    if verbose:
        print("\nMoves: %s ms vs. %s ms" % (0.000001*move_oldns, 0.000001*move_newns))
        print("OLD SYS:  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("nAccepted() = %s, nRejected() = %s" % (old_naccepted, old_nrejected))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("nAccepted() = %s, nRejected() = %s" % (new_naccepted, new_nrejected))

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()

    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    r_oldcnrg = oldsys.energy( oldff.components().coulomb() ).value() + \
                oldsys.energy( oldff14.components().intra14Coulomb() ).value()
    r_oldljnrg = oldsys.energy( oldff.components().lj() ).value() + \
                 oldsys.energy( oldff14.components().intra14LJ() ).value()

    r_newcnrg = newsys.energy( newff.components().coulomb() ).value() + \
                newsys.energy( newff14.components().intra14Coulomb() ).value()
    r_newljnrg = newsys.energy( newff.components().lj() ).value() + \
                 newsys.energy( newff14.components().intra14LJ() ).value()

    if verbose:
        print("\nRecalculated energy")
        print("OLD SYS:  %s  %s  %s  : %s ms" % (r_oldcnrg+r_oldljnrg,r_oldcnrg,r_oldljnrg,
                                                 0.000001*oldns))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (r_newcnrg+r_newljnrg,r_newcnrg,r_newljnrg,
                                                 0.000001*newns))

    assert_almost_equal( oldcnrg, r_oldcnrg, 3 )
    assert_almost_equal( oldljnrg, r_oldljnrg, 3 )

    assert_almost_equal( newcnrg, r_newcnrg, 3 )
    assert_almost_equal( newljnrg, r_newljnrg, 3 )

def test_property_map(verbose = False):

    # Load the molecular system.
    s = MoleculeParser.read(["../io/ala.top", "../io/ala.crd"])

    # Get the first molecule.
    mol = s.molecule(MolIdx(0))

    # Make the molecule editable.
    edit_mol = mol.edit()

    # Copy the coordinates property to something else.
    edit_mol.setProperty("backup", mol.property("coordinates"))

    # Delete the coordinates property.
    edit_mol.removeProperty("coordinates")

    # Commit the changes.
    mol = edit_mol.commit()

    # Map the coordinates property.
    prop_map = { "coordinates" : "backup" }

    # Create an FF object.
    intraclj = IntraFF("intraclj")

    if verbose:
        print("\nTrying to add molecule with mapped coordinates: %s" % prop_map)

    # Try to add the molecule using the property map.
    try:
        intraclj.add(mol, prop_map)

        if verbose:
            print("[SUCCESS] Added molecule.")
    except:
        raise

if __name__ == "__main__":
    test_energy(True)
    test_sim(True)
    test_property_map(True)
