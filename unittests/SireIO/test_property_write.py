
from Sire.IO import MoleculeParser
from Sire.Base import PropertyMap, StringProperty


def test_property_write():
    mols = MoleculeParser.load("../io/cholesterol.sdf")

    p = PropertyMap()
    p.set("fileformat", StringProperty("pd"))

    print(p)

    f = MoleculeParser.save(mols, "test", map=p)
    print(f)


if __name__ == "__main__":
    test_property_write()

