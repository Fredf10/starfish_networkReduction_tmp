#!/usr/bin/python
#
import sys
import subprocess

ret_val = subprocess.call([sys.executable, "UnitTesting/unitTest_singleVessel.py"])

if ret_val != 0:

    yes = ('y', 'Y', 'yes', 'Yes')
    no = ('n', 'N', 'no', 'No', 'NO')
    value = 'K'
    stdin = sys.stdin
    sys.stdin = open('/dev/tty')
    while value not in yes+no:
        value = str(raw_input("\n Tests Failed, do you still want to commit the changes? y/[N]: "))
        if str(value) in yes:
            ret_val = 1

    sys.stdin = stdin

exit(ret_val)


