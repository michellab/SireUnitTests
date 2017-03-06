
import os

import Sire.Config

def test_analyse_freenrg(verbose=False):

    if verbose:
        print("Analysing free energy...")

    output = os.popen("%s/analyse_freenrg -i ../io/freenrgs.s3" % Sire.Config.binary_directory).readlines()

    if verbose:
        print("Complete!")

    for line in output:
        if line.find("Bennetts = 0.7754") != -1:
            if verbose:
                print("Energy evaluated ok!")
                return

    if verbose:
        print("Incorrect output?")

        for line in output:
            print(line.lstrip().rstrip())

    assert(False)

if __name__ == "__main__":
    test_analyse_freenrg(True)


