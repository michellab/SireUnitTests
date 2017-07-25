
from Sire.IO import *

from nose.tools import assert_equal,assert_almost_equal

def _assert_vec_equal(vec0, vec1):
    assert_almost_equal( vec0.x(), vec1.x(), 5 )
    assert_almost_equal( vec0.y(), vec1.y(), 5 )
    assert_almost_equal( vec0.z(), vec1.z(), 5 )

def _assert_vecs_equal(vals0, vals1):
    assert_equal( len(vals0), len(vals1) )

    for i in range(0,len(vals0)):
        _assert_vec_equal( vals0[i], vals1[i] )

def test_amberrst(verbose=False):

    if verbose:
        print("Reading...")

    s = MoleculeParser.read("../io/NA16.top", "../io/NA16.rst")

    if verbose:
        print("Converting...")

    newfile = AmberRst(s)

    if verbose:
        print(newfile)
        print(newfile.creatorApplication())

    if verbose:
        print("Loading original...")

    oldfile = AmberRst("../io/NA16.rst")

    if verbose:
        print("Comparing...")

    assert_equal( newfile.nFrames(), oldfile.nFrames() )
    _assert_vecs_equal( newfile.coordinates(), oldfile.coordinates() )
    _assert_vecs_equal( newfile.velocities(), oldfile.velocities() )
    _assert_vecs_equal( newfile.forces(), oldfile.forces() )
    _assert_vec_equal( newfile.boxDimensions(), oldfile.boxDimensions() )
    _assert_vec_equal( newfile.boxAngles(), oldfile.boxAngles() )
    assert_equal( newfile.createdFromRestart(), oldfile.createdFromRestart() )
    assert_almost_equal( newfile.time(), oldfile.time() )

    # save the file and reload
    if verbose:
        print("Writing to a temporary file...")

    newfile.writeToFile("test.rst")

    if verbose:
        print("Reloading from the temporary file...")

    new2file = AmberRst("test.rst")

    if verbose:
        print("Comparing the data...")
    
    assert_equal( new2file.nFrames(), oldfile.nFrames() )
    _assert_vecs_equal( new2file.coordinates(), oldfile.coordinates() )
    _assert_vecs_equal( new2file.velocities(), oldfile.velocities() )
    _assert_vecs_equal( new2file.forces(), oldfile.forces() )
    _assert_vec_equal( new2file.boxDimensions(), oldfile.boxDimensions() )
    _assert_vec_equal( new2file.boxAngles(), oldfile.boxAngles() )
    assert_equal( new2file.createdFromRestart(), oldfile.createdFromRestart() )
    assert_almost_equal( new2file.time(), oldfile.time() )

if __name__ == "__main__":
    test_amberrst(True)

