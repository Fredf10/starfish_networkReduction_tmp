#
# Copyright 2014 Vinzenz Eck (vinzenz.g.eck@ntnu.no)
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

STARFiSh: 
is a shell-based scientific simulation program for blood flow in mamals!

vascular network creator - vnc - is a shell-based tool to create 
arterial networks for the simuation tool STARFiSh.

For programmers going to work on this code:

The standard convention for writing Python code:
https://www.python.org/dev/peps/pep-0008/

Using Sphinx with some extensions, docstrings in the code 
will be used to autogenerate documentation. To generate it, 
enter the AutoDocumentation folder and run "make html". If 
you've altered modules, run "make clean" first. If you've 
altered module names, folders, or top level modules, edit 
AutoDocumentation/sources/index.rst accordingly, and then 
run "sh fullclean.sh".

We will use a slightly modified Google standard 
for writing docstrings. (return is different)

def foo(input1, input2):
	"""
	describe function here
	
	Args:
		input1 (type): description of input1
		input2 (type): description of input2
	
	Returns:
		type of output1
			description of output1
	
		type of output2
			description of output2
	
	Raises:
		IOError: An error occured loading myClass.myStuff
	"""

Example can be found in UtilityLib/moduleCSV.py
in the function readBCFromCSV



