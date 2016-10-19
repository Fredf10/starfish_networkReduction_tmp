# List of confirmed tasks
* update Main.py to be a universal interface for simulation, network creation and visualisation (Vinz)
* Add confirmation of overwriting any networks in the working directory (Even simulation cases)  (Vinz)
* Add confirmation for rerunning samples/simulations of a UQSA case if they are already existing (Vinz)
* Update and verify the list of dependencies on Fedora (Jacob)
* Verify that template networks are working (Jacob)
* Fix issue with deleting WD
* rename vascularPolynomialChaos to vascularUQSA ?;w


# Unassigned tasks
* Check for (einar) TODOs
* Note that the WD must be given as an absolute path (or enable relative to the user's home folder (Could provide a sensible default folder in Documents a la MATLAB). But need system variables to that as OS specific
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
* flow from file inflowfile is set to absolute path. If template files are copied, a absoulte path is set (Fredrik)



# Wild TODOs
(Jacob Sturdy   2016-10-19 11:36:08 +0200   7) # created by Vinzenz Eck TODO:vinzenz.eck@mytum.de
(Jacob Sturdy   2016-10-19 11:36:08 +0200   9) # TODO: Hallvard M Nydal, Yappi, ... Fredrik Fossan, Yvan Gugler, Jacob Sturdy, Einar ..., 
(Jacob Sturdy   2016-10-19 11:36:08 +0200  10) # TODO: ADD LICENSE (MIT) and COPYRIGHT
(Jacob Sturdy   2015-02-17 15:22:54 +0000   8) ## TODO: needsto be imported as modules
(Jacob Sturdy   2016-02-09 13:31:12 +0100  59)         self.geometryType       = 'uniform'         # #TODO: fix this geometry type: 'uniform', 'cone', 'Cons'
(Jacob Sturdy   2015-02-17 15:22:54 +0000  84)         # FLUID properties TODO: annotate units
(Jacob Sturdy   2016-02-09 13:31:12 +0100 226)                 #TODO: Not sure if this is important for anything.
(Jacob Sturdy   2015-11-06 19:44:33 +0100 244)         # TODO: THESE MUST BE CORRECTED TO THE RIGHT SHAPE AS THE ADAPTIVE GRID DOESN'T RESIZE WHEN CHANGING N
(Jacob Sturdy   2015-11-06 19:44:33 +0100 245)         # TODO: The memory estimate is therefore off as well
(Jacob Sturdy   2015-11-16 10:52:58 +0100 322)         # TODO Implement h5py direct_read method to improve speed
NetworkLib/classBaroreceptor.py (Jacob Sturdy 2015-11-06 19:44:33 +0100   70)         # TODO: is this safe with inheritance?
NetworkLib/classBaroreceptor.py (Jacob Sturdy 2015-11-06 11:44:39 +0100   73)         # TODO: Extract to link Boundary condtions
SolverLib/classBaroreceptor.py  (Jacob Sturdy 2015-09-10 20:48:39 +0200  163)             self.voi, self.states, self.algebraic = self.solveCellML(strain)  # TODO Results in two solves
SolverLib/classBaroreceptor.py  (Jacob Sturdy 2015-09-10 20:48:39 +0200  186)             self.voi, self.states, self.algebraic = self.solveCellML(strain)  # TODO Results in two solves
SolverLib/classBaroreceptor.py  (Jacob Sturdy 2015-07-07 18:22:56 +0200  497)         # TODO xml parser converts these to sec-1 if I give units min-1
NetworkLib/classBaroreceptor.py (Jacob Sturdy 2015-09-21 16:59:47 +0200  619)         # TODO Automated initialization of initial conditions
NetworkLib/classBaroreceptor.py (Jacob Sturdy 2015-11-06 11:44:39 +0100  644)         # TODO: How to handle this efferent delay in the chunking method?
SolverLib/classBaroreceptor.py  (Jacob Sturdy 2015-09-10 20:48:39 +0200 1058)         # TODO: HARD CODED ONLY TO USE THORACIC AORTA!!!
(Jacob Sturdy   2015-10-02 13:58:54 +0200   94)         # TODO: Generalize residual wrappers?
(Vinzenz G. Eck 2015-10-01 15:00:00 +0200  393)         except Exception:          #TODO This should get some comment to explain what is happening.
(einarnyb       2015-07-07 13:29:17 +0200  447) # TODO wrote this before convert to spaces
(einarnyb       2015-07-07 13:29:17 +0200 1158)         # TODO: this should refer to BoundaryCondition's update function instead of doing the exact same thing.
(einarnyb       2015-07-07 13:29:17 +0200 1213)         except Exception: #TODO: should be explained.
(Jacob Sturdy   2015-11-11 13:44:34 +0100 1338)         # TODO: Polymorphic scalar/array variables?
(Jacob Sturdy   2015-11-11 14:46:18 +0100 1425)         # TODO: Polymorphic scalar/array variables?
(Jacob Sturdy   2015-12-03 14:45:54 +0100 2496)         self.atriumPressure0 = 7.5 * 133.32  # TODO: Fix this: Pressure in the atrium ## venouse pressure?!
(Jacob Sturdy   2015-10-14 11:11:22 +0200 2720)     TODO: Convert this to Napoleon compliant docstring
(Jacob Sturdy   2015-11-06 11:44:39 +0100 2787)         self.atriumPressure = np.array([7.5 * 133.32])  # TODO: Fix this: Pressure in the atrium ## venouse pressure?!
(Jacob Sturdy   2015-10-14 11:11:22 +0200 2938)         TODO: What is the physical interpretation of this? Is it nonconservative in volume? momentum?
(Jacob Sturdy   2015-10-14 11:11:22 +0200 2957)         # TODO: fix how this inherits from generalizedPQ_BC
(Jacob Sturdy   2015-10-14 11:11:22 +0200 2959)             # TODO: Ignore return values as they are stored as class members
(Jacob Sturdy   2015-10-14 11:11:22 +0200 2968)             # TODO: Anyway to factor these functions so they are only active when the heart is "coupled" to the arteries?
(Jacob Sturdy   2015-10-30 16:00:38 +0100   65)         # TODO: Remove when refactored
(Jacob Sturdy   2015-11-06 12:22:50 +0100  136)         # random variables TODO: move out of here to global class
(Jacob Sturdy   2016-10-19 11:37:17 +0200  319)                     except Exception: #TODO: Should htis exception be propagated?
(Jacob Sturdy   2016-10-19 11:37:17 +0200  347)                 if vesselId != self.root: logger.error("Error Wrong Root found") #TODO: should this stop something?
(Jacob Sturdy   2015-11-06 13:54:25 +0100  415)         # TODO: Can this be moved?
(Jacob Sturdy   2015-11-06 18:57:26 +0100  520)         # TODO: Integrate precalculated data into data saving framework
(Fredrik        2016-09-11 16:24:07 +0200  590)                 # TODO: Better way to return this to normal, while clearing the data?
(Jacob Sturdy   2015-09-02 17:48:00 +0200  625)         # TODO: Is vessels[1] the root?
(Jacob Sturdy   2015-08-28 18:12:08 +0200  630)         # TODO determine appropriate behaviour if simulation time is shorter that head up tilt time
(Jacob Sturdy   2015-09-03 10:51:57 +0200  639)         # TODO: Why is the key "1" here?
(Jacob Sturdy   2015-08-28 18:12:08 +0200  644)         # TODO: Do these belong here? and do they need to happen every simulation?
(Jacob Sturdy   2015-11-19 15:05:56 +0100  654)         # TODO: Pressure update assumes happening last
(Jacob Sturdy   2015-11-06 18:57:26 +0100  664)         # TODO: Volume calculation assumes all other objects have been updated for the current time step!!!!
(Jacob Sturdy   2015-02-17 15:22:54 +0000  674)             # TODO: Is there a better way to define these in the vessel class
(Jacob Sturdy   2015-11-06 18:57:26 +0100  687)         # TODO: add heart handling
(Jacob Sturdy   2015-10-29 16:39:42 +0100  705)         # TODO: Integrate this better with the chunking mechanism
(Jacob Sturdy   2015-01-23 15:45:05 +0000  724)         # TODO, what if this fails? do we know?
(Jacob Sturdy   2015-11-07 23:22:20 +0100  784)         # TODO: should these be errors?
(Jacob Sturdy   2016-10-19 11:37:17 +0200  832)         #TODO: return full non interpolated solution
(Fredrik        2016-09-11 16:24:07 +0200 1429)             # TODO: 1) both flow and pressure BC at inlet
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 1604)                 ## TODO: uncomment again after MC simulations are done
(Jacob Sturdy   2015-08-28 18:12:08 +0200 2143)         # TODO: what is this?
NetworkLib/classVenousPool.py (Jacob Sturdy 2015-11-19 15:05:56 +0100  30)         # TODO: Make this more efficient
SolverLib/classVenousPool.py  (Jacob Sturdy 2015-09-14 11:07:43 +0200 122)         # TODO figure out how to have BRX not need these as maximal values
SolverLib/classVenousPool.py  (Jacob Sturdy 2015-09-10 15:51:30 +0200 172)                     # TODO: Add other types
NetworkLib/classVenousPool.py (Jacob Sturdy 2015-11-06 11:44:39 +0100 207)             # TODO: Concurrency issue as Vusv may be modified by BRX on same time step?
NetworkLib/classVenousPool.py (Jacob Sturdy 2015-11-06 13:54:25 +0100 277)     # TODO: This is a hack until the top level structure is resolved fully
SolverLib/class1DflowSolver.py (einarnyb       2015-07-07 13:29:17 +0200  45)         assert isinstance(vascularNetwork, VascularNetwork) #TODO asserts get removed automatically
SolverLib/class1DflowSolver.py (einarnyb       2015-07-07 13:29:17 +0200  46)         # the vascular network to solve                     #TODO when compiled for release
SolverLib/class1DflowSolver.py (Vinzenz G. Eck 2015-11-19 16:34:59 +0100 269)                 # TODO: comment grid correction in again !!
SolverLib/class1DflowSolver.py (Jacob Sturdy   2015-11-06 18:57:26 +0100 555)         # TODO: This depends on sequential iteration through the numerical objects list!!!!
SolverLib/class1DflowSolver.py (Jacob Sturdy   2015-09-10 15:51:30 +0200 639)                         raise # TODO: why does self.exception() not force the program to quit?
SolverLib/classFields.py (Jacob Sturdy 2015-09-02 17:48:00 +0200 109)         # TODO: make sure the areas used are correct
SolverLib/classFields.py (einarnyb     2015-07-08 09:18:54 +0200 151)         #TODO: Please explain this if statement in a comment.
SolverLib/classConnections.py (einarnyb     2015-07-08 09:18:54 +0200  215)         # TODO: This was set to test sumPError, not sumPErrorNonLin
SolverLib/classConnections.py (einarnyb     2015-07-08 09:18:54 +0200  216)         # TODO: When set correctly, program wouldn't run. Commented out.
(Jacob Sturdy   2016-10-19 16:01:03 +0200  9) * Check for (einar) TODOs
(Jacob Sturdy   2016-10-19 16:01:03 +0200 11) * Check TODOs
(Jacob Sturdy   2015-11-07 23:02:27 +0100   5)     # TODO: This produces an error!! with the pretty_print argument
(Jacob Sturdy   2015-08-31 11:25:18 +0200  40)         # TODO:Should we add a field to write out non SI to xml?
(einarnyb       2015-06-23 14:39:26 +0200  82)     except: #TODO: find out what errors etree.Element raises and fix this except statement.
(Jacob Sturdy   2015-05-06 08:12:24 +0000 220)         # TODO: Check with Vinz if there's a reason to use ' ', as sep=None seems better
(einarnyb       2015-06-23 14:39:26 +0200 276)                 #TODO: fix exception handling here
(einarnyb       2015-06-23 14:39:26 +0200 282)                 #TODO: fix exception handling here
(Vinzenz G. Eck 2015-07-07 09:23:33 +0200 399)                                     ## TODO: fix problem when loading an absolute path
(Jacob Sturdy 2015-11-06 18:58:31 +0100  62)             ValueError: (TODO: NOT IMPLEMENTED) raise error if an unappropriate dictObjType is passed
(Jacob Sturdy 2015-11-06 18:58:31 +0100  72)             #TODO: approptiate testing for dictObjType  
(Jacob Sturdy 2015-11-06 18:58:31 +0100 274)         ## TODO: check if there is not-needed information in the xml file
(Jacob Sturdy 2015-11-06 18:58:31 +0100 503)             # TODO: exchange with appropriate warning exception
(einarnyb       2015-06-23 16:03:03 +0200 230) # TODO: (einar) fix exception variable
(Vinzenz G. Eck 2015-08-24 15:14:36 +0200 130)     # TODO: add units
(Vinzenz G. Eck 2015-08-24 15:14:36 +0200 174)                 # TODO: add units
(Vinzenz Eck    2015-04-17 08:24:43 +0000 193)     # TODO: read write polynomial chaos variables
(Vinzenz Eck    2015-04-17 08:24:43 +0000 298)         # TODO: read write polynomial chaos variables
(einarnyb       2015-06-23 16:03:03 +0200  28) # TODO: (einar) Rename input variable exception and corresponding strings
(einarnyb       2015-06-23 16:03:03 +0200 137) # TODO: (einar) fix exception variable in function
(Jacob Sturdy 2015-11-06 18:58:31 +0100 158)         # TODO: This depends on being the last object called for the time step!!!!
(Jacob Sturdy   2015-11-11 14:46:18 +0100  59) # TODO what's the deal with this?
(Jacob Sturdy   2015-08-28 15:48:36 +0200 282)                         'communicators'         : communicatorReference, #TODO why is this Reference?
(Vinzenz G. Eck 2015-11-24 16:08:32 +0100  306)         # TODO!!
(Vinzenz G. Eck 2015-12-15 16:00:49 +0100  439)         # TODO!!
(Vinzenz G. Eck 2015-12-16 14:24:07 +0100  578)         # TODO!!
(Vinzenz        2015-12-17 16:37:59 +0100  720)         # TODO!!
(Vinzenz        2015-11-27 15:04:33 +0100  841)         # TODO!!
(Vinzenz        2015-12-17 16:37:59 +0100  965)         # TODO!!
(Vinzenz        2015-12-17 16:37:59 +0100 1082)         # TODO!!
(Vinzenz G. Eck 2015-10-09 14:36:31 +0200 1272)             # TODO: implement dependent dist methods for MC
(Vinzenz        2015-12-17 16:37:59 +0100 1368)                 # TODO: implement dependent dist methods for MC
(Vinzenz G. Eck 2015-11-18 10:36:15 +0100 144)         ## TODO: fix data saving methods as it will not work now
(Vinzenz G. Eck 2015-10-05 16:20:54 +0200  73)             ValueError: (TODO: NOT IMPLEMENTED) raise error if an unappropriate dictObjType is passed
(Jacob Sturdy   2015-11-06 13:54:25 +0100  83)             #TODO: appropriate testing for dictObjType  
(Vinzenz G. Eck 2015-09-03 11:52:21 +0200 284)         ## TODO: check if there is not-needed information in the xml file
(Vinzenz G. Eck 2015-11-18 10:36:15 +0100 553)                 # TODO: use externalVariable definitions dictionary ... #
(Vinzenz G. Eck 2015-11-18 10:36:15 +0100 36)             # check for trajectory stuff TODO: do it nicer
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 44)             #TODO: implement abc sampling hash if monte carlo is used
(Vinzenz G. Eck 2015-10-21 17:32:52 +0200 101)                     #TODO: create own class tha xVal is not missued as number of space points per vessel
(Vinzenz        2015-11-27 15:04:33 +0100 278)                     #TODO: create own class tha xVal is not missued as number of space points per vessel
(Vinzenz G. Eck 2015-06-15 15:58:31 +0200 356)             ## TODO: Jacob add here your stuffls
(Vinzenz Eck    2015-04-14 14:20:45 +0000 359)             ## TODO: add more locations as necessary baroreceptor etc.
(Vinzenz G. Eck 2015-10-21 17:32:52 +0200 599)             #TODO: save the figures and close ..
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100  48)         # TODO: if random sample, ensure that sub set samples of sample space are randomly choosen
(Vinzenz        2015-12-17 16:37:59 +0100  63)         #TODO: check if samplesSize and offsett+sampleSize is within range of samples
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100  70)             #TODO: implement abc sampling hash if monte carlo is used
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 114)             #TODO: implement abc sampling hash if monte carlo is used
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 193)         # TODO: type conversion should be moved to load function
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 211)         # TODO: if random sample, ensure that sub set samples of sample space are randomly choosen          
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 215)             # TODO: dependent Samples
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 231)         ## TODO: rewrite in more base-class manner
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 232)         # TODO: datatypeconversion!!
(Vinzenz G. Eck 2015-11-05 11:25:00 +0100 233)         # TODO: optional bool for object types
(Vinzenz G. Eck 2015-10-08 15:07:52 +0200  76)         self.localEvaluation        = True #TODO: add functions for server
(Vinzenz        2015-11-27 15:04:33 +0100 153)         #TODO: replace create evaluation files or save them to disc!!!
(Vinzenz G. Eck 2015-11-18 10:36:15 +0100 210)             #TODO rename caseName!!
(Vinzenz G. Eck 2015-10-09 14:36:31 +0200 125)             # TODO: server simulations not implemented yet
VisualisationLib/class3dVisualisation.py (einarnyb     2015-07-08 09:33:38 +0200  698)                     #TODO: Try Except Pass should be fixed
VisualisationLib/class3dVisualisation.py (einarnyb     2015-07-08 09:33:38 +0200  707)                     #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200  143)                     #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200  152)                         #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200  169)                 #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200  596)             #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200  617)         #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200  635)         #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200  690)                     #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200 1275)             #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200 1343)                 #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200 1425)             #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200 1438)             #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200 1591)         #TODO: Try Except Pass should be fixed
VisualisationLib/class2dVisualisation.py (einarnyb       2015-07-08 09:33:38 +0200 1887)         #TODO: Try Except Pass should be fixed
VisualisationLib/class3dControlGUI.py (einarnyb    2015-07-08 09:33:38 +0200 596)             #TODO: Try Except Pass should be fixed
VisualisationLib/classRealTimeVisualisation.py (einarnyb     2015-07-08 09:33:38 +0200  55)         #TODO: Try Except Pass should be fixed
VisualisationLib/classRealTimeVisualisation.py (einarnyb     2015-07-08 09:33:38 +0200  68)         #TODO: Try Except Pass should be fixed
VisualisationLib/classRealTimeVisualisation.py (einarnyb     2015-07-08 09:33:38 +0200  79)         #TODO: Try Except Pass should be fixed
VisualisationLib/classRealTimeVisualisation.py (einarnyb     2015-07-08 09:33:38 +0200  90)         #TODO: Try Except Pass should be fixed
(Vinzenz Eck    2014-08-04 10:28:49 +0000  964)         # TODO: handle charset
(Vinzenz Eck    2014-08-04 10:28:49 +0000 1203)                 # TODO: handle start/end points
