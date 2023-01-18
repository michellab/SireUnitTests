try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *


def test_netcdf(verbose=False):

    r = AmberRst("../io/NA16.rst")

    print(r)

    print("Writing back to netcdf...")
    r.writeToFile("file_netcdf_test.rst")

    print("...done")

    print(r.creatorApplication())
    print(r.time())

    s = MoleculeParser.read("file_netcdf_test.rst", "../io/NA16.top")

    r7 = AmberRst7(s)
    r7.writeToFile("file_netcdf_test2.rst7")

    print(r7)

    PDB().write(s.molecules(), "file_netcdf_test.pdb")


if __name__ == "__main__":
    test_netcdf(True)
