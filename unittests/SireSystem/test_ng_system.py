
import pytest

def test_ng_system():
    try:
        import Sire as sr
        mols = sr.load("../io/ala.*")
    except Exception:
        print("Skipping as using an older version of Sire")
        return

    assert mols["residues with element O"] == mols["element O"].residues()
    assert mols["element O"] == mols.atoms("element O")
    assert mols.num_atoms() == len(mols.atoms())
    assert mols.num_residues() == len(mols.residues())

    with pytest.raises(KeyError):
        mols.chains()

    assert mols.num_chains() == 0

    with pytest.raises(KeyError):
        mols.segments()

    assert mols.num_segments() == 0


if __name__ == "__main__":
    test_ng_system()

