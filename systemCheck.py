import os
def systemCheck():
    """
    Check for necessary moduls and 3thrd party libarys needed for STARFiSh and vnc:
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
    print ""
    print "Check for necessary moduls and 3thrd party libarys needed for vascular1DFlow and vnc"
    print ""
    
    print "check for matplotlib ",
    try:
        import matplotlib.pyplot 
        print installed
    except:
        print " IMPORT ERROR; no version of matplotlib.pyplot found"
    
    print "check for ffmpeg     ",
    if os.path.isfile("/bin/ffmpeg"):
        print installed
    else:
        print " ERROR could not find ffmpeg installation; this needs to be manual"
        print " (not supported by the automatic installer)"
        print " for Fedora:"
        print "            1.) add RPMfusion free repositories to yum:"
        print "                 http://rpmfusion.org/Configuration "
        print "                 or"
        print "                 sudo yum install --nogpgcheck http://download1.rpmfusion.org/free/fedora/rpmfusion-free-release-$(rpm -E %fedora).noarch.rpm http://download1.rpmfusion.org/nonfree/fedora/rpmfusion-nonfree-release-$(rpm -E %fedora).noarch.rpm"
        print "            2.) sudo yum install ffmpeg"     
        
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
    
    print "check for psutil"
    try:
        import psutil
        print installed
    except:
        print " IMPORT ERROR; no version of psutil found"
    
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
        
    print "check for h5py       ",
    try:
        import h5py
        print installed
    except:
        print " IMPORT ERROR; no version of h5py found"  
    
    print "check for sphinx     ",
    try:
        import sphinx
        print installed
    except:
        print " IMPORT ERROR; no version of sphinx found"  
    
    print "check for reportlab  ",
    try:
        import reportlab
        print installed
    except:
        print " IMPORT ERROR; no version of reportlab found"         

if __name__ == '__main__':
    systemCheck()