try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

import os
import subprocess

import Sire.Config


def test_analyse_freenrg(verbose=False):

    if verbose:
        print("Analysing free energy...")

    if os.path.exists("%s/analyse_freenrg" % Sire.Config.binary_directory):
        p = subprocess.Popen(
            (
                "%s/analyse_freenrg -i ../io/freenrgs.s3"
                % Sire.Config.binary_directory
            ).split(),
            stdout=subprocess.PIPE,
        )
    else:
        p = subprocess.Popen(
            (
                "%s/sire_python %s/scripts/analyse_freenrg.py -i ../io/freenrgs.s3"
                % (Sire.Config.binary_directory, Sire.Config.share_directory)
            ).split(),
            stdout=subprocess.PIPE,
        )
    output, _ = p.communicate()
    output = output.decode("UTF-8").split("\n")

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

    assert False


if __name__ == "__main__":
    test_analyse_freenrg(True)
