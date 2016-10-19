# List of confirmed tasks
* update Main.py to be a universal interface for simulation, network creation and visualisation (Vinz)
* Add confirmation of overwriting any networks in the working directory (Even simulation cases)  (Vinz)
* Add confirmation for rerunning samples/simulations of a UQSA case if they are already existing (Vinz)
* Update and verify the list of dependencies on Fedora (Jacob)
* Verify that template networks are working (Jacob)

# Unassigned tasks
* Check for (einar) TODOs
* Note that the WD must be given as an absolute path (or enable relative to the user's home folder (Could provide a sensible default folder in Documents a la MATLAB).
* Check TODOs
* Delete central/minimum VenousPressure tags from xmlTemplates
* Test venous system
* Test BRX
* Test gravity
* Test heart models
* Test UQSA
* Test Anastomosis

# List of components to fully separate
* Baroreflex
* Venous system
* Heart model data class
* Gravity

# List of ideas
* update path manipulations to use `os.path` 
* Make sure logging can be "quiet"
* Replace printouts with logging 
* (VNC) always able to "cancel" operation in VNC
* (VNC) Inform about geometry of vessel/how to change this
* (VNC) Missing As in Hayashi default


# List of done tasks
* Venous system
* Update and verify the list of dependencies on Ubuntu (Jacob)
* Fix optional initiation of venous system (Jacob)
* Get website privileges for Fredrik, Jacob and Leif (Assigned to:Vinz)
* Fix WD manager so that WD's are successfully created from command line (Vinz)
* Ensure that the config file is generated when first running the code (Vinz)
* Fix vnc to not printout options for CSV BCs and RVs (Assigned to:Vinz)
* VNC never should crash while loading networks and vessels from xml and csv (delimiter still problematic). (Vinz)
* VNC shall detect delimiter and either warn/autofix rather than crash. (Vinz)
* flow from file inflowfile is needed only in the "networkXMLFile" directory (Fredrik)
