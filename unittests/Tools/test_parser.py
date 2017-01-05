
from Sire.Tools import resolveParameters,readParams
from Sire.Units import *

params = readParams("configfile")

params["test1"] = "5 angstrom"
params["test2"] = "0.5 * gram/(centimeter*centimeter*centimeter)"
params["test3"] = 500

@resolveParameters
def test_resolveparams():
    print("Hello")

if __name__ == "__main__":
    test_resolveparams(params)

