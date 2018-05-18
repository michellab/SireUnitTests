
from Sire.Mol import *

try:
    s = Select()
    have_select = True
except:
    have_select = False

def test_selections(verbose=False):
    if not have_select:
        return

    if verbose:
        print("\nBase name selections\n")

    s = Select("atomnam hello")
    s = Select("resname goodbye")
    s = Select("chainname world")
    s = Select("(atomnam bracket)")
    s = Select("((atomnam doublebracket))")
    s = Select("(atomnam one and resnam two) and (resnam three)")
    s = Select("(atomnam one and resnam two) and (resnam three and atomname four)")
    s = Select("atomnam hello; resname goodbye")
    s = Select("atomnam hello and resname goodbye")
    s = Select("atomnam hello and resname goodbye and chainname world")
    s = Select("atomnam hello and resname goodbye and (chainname world and resname one)")
    #s = Select("atomnam hello and (resname goodbye and (chainname world and (resname one and atomname two)))")

    if verbose:
        print("\nMulti name\n")

    s = Select("atomname C,CA , N,O")
    s = Select("atomname 'C 1',O")
    s = Select("atomname /hello/")
    s = Select("atomname /goodbye/i")
    s = Select("atomname C,/hello/,'C 1',/goodbye/i")
    s = Select("(atomname /CA/i or atomname /C/i) and resname /ALA/i")

    if verbose:
        print("\nNumbers...\n")

    s = Select("atomnum 5")
    s = Select("resnum 1:10")
    s = Select("cgidx 1:10:2, 5, 2:10")
    s = Select("resnum 100 and chainname A")
    s = Select("resnum > 5 and resnum <= 100")

    if verbose:
        print("\nWith/In/not...\n")

    s = Select("molecules with (resname /ala/i,/asp/i,/gly/i and atomname /ca/i)")
    s = Select("atoms in residx 1:10")
    s = Select("not resname ALA")
    s = Select("chainname B and not atomname CA,C,N,O")

    if verbose:
        print("\nSubscripting\n")

    s = Select("{atoms in resname ALA}[0:10]")
    s = Select("atoms in {resname ALA}[0:10]")
    s = Select("{resname ALA and atomname /ca/i}[5]")
    s = Select("{molecules with resname /ala/i}[-1]")

    if verbose:
        print("\nComments\n")

    s = Select("resname /* This is a comment */ ALA")
    s = Select("resname /ala/i;\n//another comment\n resname /gly/i")
    s = Select("{molecules /*comment*/ with resname /ala/i /*comment */}/*comment*/[-1]")

    if verbose:
        print("\nWithin distance\n")

    s = Select("molecules within 5.0 of resname ALA")
    s = Select("{atoms within 10 nm of molecules with resname /ala/i}[0:-1:100]")

    if verbose:
        print("\nUser-supplied tokens\n")

    Select.setToken("protein", "molecules with resname /ala/i,/asp/i,/arg/i,/leu/i")

    s = Select("protein")

    Select.setToken("ligand", "resname /lig/i")

    s = Select("protein or ligand")
    s = Select("{protein}[0]")

    if verbose:
        print("\nwhere parsing\n")

    s = Select("molecules where coords.min > 5,3,2")
    s = Select("atoms where center within 3 of resname /ala/i")
    s = Select("atoms where center <= (1,2,3)")

if __name__ == "__main__":
    test_selections(True)

