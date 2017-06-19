#!/usr/bin/python
#
from __future__ import print_function
import sys
import subprocess

ret_val = subprocess.call([sys.executable, "UnitTesting/unitTest_singleVessel.py"])
if ret_val != 0:
    msg = "ERROR: 1 or more tests failed."
    msg += "\n\nUse the '--no-verify' option to by-pass this hook"
    sys.exit(msg)
sys.exit(ret_val)


