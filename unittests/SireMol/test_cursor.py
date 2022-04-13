
import Sire as sr

try:
    c = sr.Mol.Cursor()
    has_cursor = True
    mols = sr.load("../io/3NSS.pdb")
except Exception:
    print("Skipping as no Cursor support in this version of Sire")
    has_cursor = False

def test_cursor():
    if not has_cursor:
        return

    mol = mols[0]

    assert mol.nChains() == 3

    c = mol.cursor()

    assert c.type() == "molecule"

    assert "cat" not in c

    c["cat"] = "meow"

    assert "cat" in c
    assert c["cat"].asAString() == "meow"

    c = c.atom(0)

    assert c.type() == "atom"

    assert "dog" not in c

    c["dog"] = "woof"

    assert "dog" in c
    assert c["dog"] == "woof"

    c = c.parent().atom(1)

    assert c.parent().type() == "residue"
    assert c.parent().parent().type() == "chain"
    assert c.type() == "atom"

    assert "dog" in c
    assert c["dog"] == ""

    c = c.parent().atom(0)

    assert "dog" in c
    assert c["dog"] == "woof"

    c = c.parent().atom(1)
    
    try:
        c["dog"] = ["bark", "growl"]
    except TypeError:
        pass

    del c["dog"]

    assert "dog" not in c
    
    c["dog"] = ["bark", "growl"]

    assert "dog" in c
    assert c["dog"][0] == "bark"
    assert c["dog"][1] == "growl"

    c = c.parent().atom(0)

    assert len(c["dog"]) == 0


def test_nav():
    if not has_cursor:
        return

    c = mols[0].cursor()

    mol = mols[0]

    assert c.type() == "molecule"

    assert len(c.chains()) == mol.nChains()

    ids = []

    for chain in c.chains():
        assert chain.type() == "chain"
        assert chain.name() == mol.chain(chain.ID()).name()
        
        assert len(chain.residues()) == mol.nResidues(chain.ID())

        chain["fish"] = f"bubble-{chain.ID()}"
        ids.append((chain.ID(), "fish", "bubble"))

        # too slow (and not needed) to do all of the residues
        for residue in chain.residues()[0::200]:
            assert residue.type() == "residue"
            assert residue.name() == mol.residue(residue.ID()).name()
            assert residue.number() == mol.residue(residue.ID()).number()
            assert residue.parent().type() == "chain"
            assert residue.parent().ID() == chain.ID()

            assert len(residue.atoms()) == mol.nAtoms(residue.ID())

            residue["dog"] = f"woof-{residue.ID()}"
            ids.append((residue.ID(), "dog", "woof"))

            for i, atom in enumerate(residue.atoms()):
                assert atom.type() == "atom"
                assert atom.name() == mol.atom(atom.ID()).name()
                assert atom.number() == mol.atom(atom.ID()).number()
                assert atom.parent().type() == "residue"
                assert atom.parent().ID() == residue.ID()
                assert atom.parent().atom(i).ID() == atom.ID()

                atom["cat"] = f"meow-{atom.ID()}"
                ids.append((atom.ID(), "cat", "meow"))

    mol = c.commit()

    assert mol.hasProperty("cat")
    assert mol.hasProperty("dog")
    assert mol.hasProperty("fish")

    for id, animal, sound in ids:
        c = mol[id].cursor()
        assert c[animal] == f"{sound}-{id}"

if __name__ == "__main__":
    test_cursor()
    test_nav()

