def systemCheck():
    '''
    Check for necessary moduls and 3thrd party libarys needed for STARFiSh and vnc:
        matplotlib
        numpy
        scipy
        #mayavi2
        #tvtk
        lxml2
        pydot
        gtk
        pyglet
        mencoder
        OpenGL
    '''
    matplotlibInstalled = False
    installed = " ... installed!"
    print ""
    print "Check for necessary moduls and 3thrd party libarys needed for vascular1DFlow and vnc"
    print ""
    
    print "check for matplotlib ",
    try:
        import matplotlib.pyplot 
        print installed
        matplotlibInstalled = True
    except:
        print " IMPORT ERROR; no version of matplotlib.pyplot found"
    
    print "check for mencoder   ",
    if matplotlibInstalled:
        import matplotlib.animation as animation
        
        if animation.MencoderWriter.isAvailable():
            print installed
        else:
            print " IMPORT ERROR; no version of mencoder found"     
    else:
        print " IMPORT ERROR; no version of matplotlib.pyplot found thus cannot check for mencoder"     
        
    print "check for numpy      ",
    try:
        import numpy
        print installed
        
    except:
        print " IMPORT ERROR; no version of numpy found"
    
    print "check for scipy      ",
    try:
        import scipy
        print installed
    except:
        print " IMPORT ERROR; no version of scipy found"
    
#     print "check for  mayavi2 ..."
#     try:
#         import mayavi
#         print installed
#     except:
#         try: 
#             import enthought.mayavi
#             print installed," [Old Version !!]"
#         except:
#             print "IMPORT ERROR; no version of mayavi2 found"
#     
#     print "check for  tvtk ..."
#     try:
#         import tvtk
#         print installed
#     except:
#         try:
#             import enthought.tvtk.api
#             print installed," [Old Version !!]"
#         except:
#             print "IMPORT ERROR; no version of tvtk found"
    
    print "check for lxml2      ",
    try:
        import lxml
        print installed
    except:
        print " IMPORT ERROR; no version of lxml2 found"
    
    print "check for pydot      ",
    try:
        import pydot
        print installed
    except:
        print " IMPORT ERROR; no version of pydot found"
    
    print "check for gtk        ",
    try:
        import gtk
        print installed
    except:
        print " IMPORT ERROR; no version of gtk found" 
    
    print "check for pyglet     ",
    try:
        import pyglet
        print installed
    except:
        print " IMPORT ERROR; no version of pyglet found" 
        
    print "check for pyOpenGL   ",
    try:
        import OpenGL
        print installed
    except:
        print " IMPORT ERROR; no version of pyOpenGL found" 
# 
#     print "check for xdot ...",
#     try:
#         import xdot
#         print installed
#     except:
#         print "IMPORT ERROR; no version of xdot found (http://pypi.python.org/pypi/xdot )" 
        
    print "check for pydot      ",
    try:
        import pydot
        print installed
    except:
        print "IMPORT ERROR; no version of found (http://pypi.python.org/pypi/pydot/1.0.2)" 


if __name__ == '__main__':
    systemCheck()