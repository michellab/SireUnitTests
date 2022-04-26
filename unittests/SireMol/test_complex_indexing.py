
import pytest

def test_complex_indexing():
    try:
        import Sire as sr
        mols = sr.load("../io/p38.pdb")
        mol = mols[0]
    except Exception:
        print("Skipping as older version of Sire")
        return

    assert mol.atom(0).index() == sr.Mol.AtomIdx(0)
    assert mol[0] == mol.atom(0)
    
    with pytest.raises(KeyError):
        mol.atom("CA")

    s = mol.atoms("CA")

    assert len(s) == 351
    assert type(s) is not list
    assert s.what() == "SireMol::Selector<SireMol::Atom>"

    assert s[0].name().value() == "CA"
    assert mol[s[0].index()].index() == s[0].index()

    assert s == mol["atomname CA"]
    assert s == mol["CA"]
    assert s == mol[sr.atomid("CA")]

    s = mol.residues("resnum <= 10")
    
    assert len(s) == 10
    assert type(s) is not list
    assert s.what() == "SireMol::Selector<SireMol::Residue>"

    assert s[0].name().value() == "GLU"

    for i in range(0, 10):
        assert s[i].number().value() == i+1

    assert s == mol["resnum <= 10"]
    assert s == mol["resnum 1, 2, 3, 4, 5, 6, 7, 8, 9, 10"]

    assert mol["atoms in resnum 1"] == mol.residue(0).atoms()

    s = mol["residues with atomname HZ"]

    assert len(s) == 13

    for res in s:
        assert res.name().value() == "PHE"
        assert res["HZ"].name().value() == "HZ"
        assert res.atom("HZ").residue() == res

    from Sire.Mol import AtomName, AtomNum

    names = [AtomName('N'), AtomName('CA'), AtomName('C'), AtomName('O'), AtomName('CB')]

    assert mol[0:5].names() == names
    assert mol[4::-1].names() == names[4::-1]

    numbers = [AtomNum(7), AtomNum(9), AtomNum(11)]

    assert mol[6:12:2].numbers() == numbers
    assert mol[10:4:-2].numbers() == numbers[2::-1]

    with pytest.raises(KeyError):
        mol["X"]

    with pytest.raises(IndexError):
        mol[10000]

    s = mol["resnum >= 5 and resnum <= 10"]

    assert len(s) == 6

    for i in range(0, 6):
        assert s[i].number().value() == i+5

    s = mol[sr.resid(idx=0)]["atomname C, CA"]

    assert len(s) == 2

    assert s["CA"].name().value() == "CA"
    assert s["C"].name().value() == "C"

    count = 0

    for res in mol["residx < 3"]:
        for atom in res["atomname C, N"]:
            count += 1

    assert count == 6


if __name__ == "__main__":
    test_complex_indexing()

