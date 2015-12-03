import sys
import StringIO

class ProgressBar(object):
    """
    
    
    Examples:
        create a progress bar with 20 elements for a loop with 200 iterations 
        and loop over the 200 iterations
    
        >>> progressBar = ProgressBar(20, 200)
        [                    ]\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b
        >>> for i in xrange(200): progressBar.progress(i)
        ####################
        <BLANKLINE>
    
        create a progress bar with 2 elements for a loop with 12 iterations 
        and loop over the 12 iterations
    
        >>> progressBar = ProgressBar(2, 12)
        [  ]\b\b\b
        >>> for i in xrange(12): progressBar.progress(i)
        ##
        <BLANKLINE>
        
        create a progress bar with 4 elements for a loop with 10 iterations 
        check the print out after iteration 5, 7 and 10 respectively
        
        >>> progressBar = ProgressBar(4, 10)
        [    ]\b\b\b\b\b
        >>> progressBar.progress(5)
        ##
        >>> progressBar.progress(7)
        #
        >>> progressBar.progress(10)
        #
    """
    def __init__(self, nElements, nIterations):
        """
        Initialization of the the progress bar, prints out the brackets of 
        the progress bar immediately
        
        All outputs to stdout will be captured and displayed after the 
        progress bar is fully loaded.

        Args:
            nElements (int): number of elements of the progress bar
            nIterations (int): number of iterations of the process which 
            progress should be displayed
        """
        # handle for stdout
        self.stdout = sys.stdout
        # buffer for stdout outputs
        self.outputBuffer = StringIO.StringIO()
        # redirect system stdout to buffer
        sys.stdout = self.outputBuffer
        self.finished = False

        self.nElements = nElements
        self.nIterations = float(nIterations)
        self.eCount = 1.
        
        # write empty brackets of the loading bar
        spaces = ''.join([' ' for i in xrange(self.nElements)])
        loadingBar = ''.join(['[',spaces,']'])
        self.stdout.write(loadingBar)
        self.stdout.flush()
        # write backspaces to the start of the loading bar
        backspacing = ''.join(['\b' for i in xrange(nElements+1)])
        self.stdout.write(backspacing)
        
    def progress(self, currentIteration):
        """
        checks the progress of the adjoin process, and advances the progress 
        bar # if necessary
        
        Args:
            currentIteration (int): number of the current iterations of the 
            adjoin process
        """
        if self.finished == False:

            # completed tasks of the process in precent
            loopPercentage = currentIteration/(self.nIterations-1.) 
            # completed shown progress in percent
            barPercentage = self.eCount/self.nElements
            
            while (barPercentage <= loopPercentage):
                self.stdout.write("#")
                self.stdout.flush()
                self.eCount += 1.
                barPercentage = self.eCount/self.nElements
            
            if loopPercentage == 1:
                self.stdout.write("\n")
                self.stdout.flush()
                sys.stdout = self.stdout
                sys.stdout.write(self.outputBuffer.getvalue())
                sys.stdout.write("\n")
                self.outputBuffer.close() 
                self.finished = True
if __name__ == '__main__':
    
    import doctest
    doctest.testmod()