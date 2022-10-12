from sire.legacy.IO import *

def test_netcdf(verbose=False):

    r = AmberRst("../io/NA16.rst")

    print(r)

    print("Writing back to netcdf...")
    r.writeToFile("test.rst")

    print("...done")

    print(r.creatorApplication())
    print(r.time())

    s = MoleculeParser.read("test.rst", "../io/NA16.top")

    r7 = AmberRst7(s)
    r7.writeToFile("test.rst7")

    print(r7)

    PDB().write(s.molecules(), "test.pdb")

if __name__ == "__main__":
    test_netcdf(True)

