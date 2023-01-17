
import os
import sys
import tarfile

try:
    import sire
    sire.use_old_api()
except ImportError:
    pass

from Sire import try_import

import Sire.Config

nose = try_import("nose")

import Sire.Base

#print("\nRunning C++ unit tests...\n")
#try:
#    Sire.Base.UnitTest.runAll()
#except:
#    pass
#
#print("\nNow running Python-based unit tests...\n")

old_cwd = os.path.abspath(".")
testdir = "."

if len(sys.argv) > 1:
    testdirs = sys.argv[1:]
    sys.argv = [sys.argv[0]]

    drop = []

    for dir in testdirs:
        if dir.startswith("-"):
            drop.append(dir)

    if len(drop) > 0:
        testdirs = []

        dirs = os.listdir(testdir)

        for dir in dirs:
            for d in drop:
                if d[1:] != dir:
                    testdirs.append(dir)

else:
    testdirs = os.listdir(testdir)

failures = []

for file in testdirs:
    subdir = os.path.join(testdir, file)

    if not file.startswith(".") and os.path.isdir(subdir):
        print("\n\nRunning tests in directory %s..." % subdir)
        print("############################################")
        os.chdir(subdir)
        success = nose.run()

        if not success:
            failures.append(file)

        os.chdir(old_cwd)

if len(failures) > 0:
    print("\n\nWARNING: SOME OF THE TEST JOBS FAILED!!!")
    print("#########################################")

    s = []

    for failure in failures:
        print("One of more jobs in %s failed!" % failure)
        s.append("One of more jobs in %s failed!" % failure)

    print("\nEXITING AS FAIL")
    error_message = "SOME OF THE UNIT TESTS FAILED!\n%s" % \
                         "\n".join(s)
    raise AssertionError(error_message)
else:
    print("\n\nHOORAY - ALL OF THE UNIT TESTS PASSED!!!")
    print("\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/")
