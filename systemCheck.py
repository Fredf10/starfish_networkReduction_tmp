from __future__ import print_function
import os
def systemCheck():
    """
    Check for necessary modules and 3thrd party libarys needed for STARFiSh: 
        matplotlib
        numpy
        scipy
        lxml2
        pydot
        gtk
        pyglet
        mencoder
        OpenGL
    """
    installed = " ... installed!"
    print("")
    print("Check for necessary modules and 3thrd party libarys needed for STARFiSh")
    print("")
    
    print("check for matplotlib ", end=' ')
    try:
        import matplotlib.pyplot 
        print(installed)
    except:
        print(" IMPORT ERROR; no version of matplotlib.pyplot found")
    
    print("check for numpy      ", end=' ')
    try:
        import numpy
        print(installed)
    except:
        print(" IMPORT ERROR; no version of numpy found")
    
    print("check for scipy      ", end=' ')
    try:
        import scipy
        print(installed)
    except:
        print(" IMPORT ERROR; no version of scipy found")
    
    print("check for psutil")
    try:
        import psutil
        print(installed)
    except:
        print(" IMPORT ERROR; no version of psutil found")
    
    print("check for lxml2      ", end=' ')
    try:
        import lxml
        print(installed)
    except:
        print(" IMPORT ERROR; no version of lxml2 found")
    
    print("check for pydot      ", end=' ')
    try:
        import pydot
        print(installed)
    except:
        print(" IMPORT ERROR; no version of pydot found")
    
    print("check for gtk        ", end=' ')
    try:
        import gtk
        print(installed)
    except:
        print(" IMPORT ERROR; no version of gtk found") 
    
    print("check for pyglet     ", end=' ')
    try:
        import pyglet
        print(installed)
    except:
        print(" IMPORT ERROR; no version of pyglet found") 
        
    print("check for pyOpenGL   ", end=' ')
    try:
        import OpenGL
        print(installed)
    except:
        print(" IMPORT ERROR; no version of pyOpenGL found") 
        
    print("check for h5py       ", end=' ')
    try:
        import h5py
        print(installed)
    except:
        print(" IMPORT ERROR; no version of h5py found")  
    
    
if __name__ == '__main__':
    systemCheck()
