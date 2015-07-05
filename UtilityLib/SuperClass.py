# SuperClass.py (name should maybe be changed)
import sys, os
import traceback
import time

cur = os.path.dirname( os.path.realpath( __file__ ) )
root = ''.join([cur,'/../'])
logFilePath = root + 'warninglog.txt'

class STARFiSh(object):
    """Super class every other class in STARFiSh 
    will inherit from. 
    
    Contains:
    Global Warning function (for handling exceptions)
    
    Planned features: 
    Global Update function for dictionaries.
    """

    def warn(self, infoString = None, exceptionHappened = True, 
             quiet = False, verbose = False, saveToFile = False):
        """ 
        Global Warning Function
        
        Args: 
            infoString (String): String describing why the warning was made.
            exceptionHappened (Optional Bool): Boolean determining if it 
                should use the last exception's info or not. 
                Defaults to True.
            quiet (Opt. Bool): Whether to print warning at all. 
                Defaults to False.
            verbose (Opt. Bool): Whether or not to print full traceback.
                Defaults to False.
            saveToFile (Opt. Bool): Whether or not to save info to file.
                Will save info in STARFiSh/warninglog.txt
                Defaults to False.
        """

        completeString = "Warning: "
        if infoString != None:
            completeString = completeString + infoString
        else: 
            completeString = completeString + "no string passed to warn()"

        if exceptionHappened:
            try: raise
            except:
                exceptionInfo = "\n" + traceback.format_exc()
        else: exceptionInfo = " "

        if verbose:
            completeString = completeString + exceptionInfo

        if not quiet:
            print completeString
        
        if saveToFile:
            if not verbose: completeString = completeString + exceptionInfo
            try: f = open(logFilePath, 'a')
            except: f = open(logFilePath, 'w')
            f.write("\n======================")
            f.write("\n{}".format(
                time.strftime("%Y.%b.%d, %H:%M:%S")))
            f.write("\n----------------------\n")
            f.write(completeString)
            f.write("======================\n\n")
            f.close()

    def errorAppendInfo(self, infoString = None):
        """
        Global Exception information appending. 
        Warning: May mess up debugger, should be checked if this is true
        before employing it throughout code.
        Warning: Using this function will make a traceback appear at the 
        line this function was called. This is unavoidable behaviour.

        Args:
            infoString (string): String containing additional information.
        """
        if infoString == None:
            raise
        else: 
            try: raise
            except Exception as e:
                raise type(e), type(e)(e.message + "\nAppended : " + infoString),sys.exc_info()[2]










