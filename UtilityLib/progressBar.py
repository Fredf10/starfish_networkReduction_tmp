import sys

class ProgressBar(object):
    """
    
    
    Examples:
        create a progress bar with 20 elements for a loop with 200 iterations and loop over the 200 iterations
    
        >>> progressBar = ProgressBar(20, 200)
        [                    ]\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b
        >>> for i in xrange(200): progressBar.progress(i)
        ####################
    
        create a progress bar with 2 elements for a loop with 12 iterations and loop over the 12 iterations
    
        >>> progressBar = ProgressBar(2, 12)
        [  ]\b\b\b
        >>> for i in xrange(12): progressBar.progress(i)
        ##
    
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
        Initialization of the the progress bar, prints out the brackets of the progress bar immediately
        
        Args:
            nElements (int): number of elements of the progress bar
            nIterations (int): number of iterations of the process which progress should be displayed
        """
        self.nElements = nElements
        self.nIterations = float(nIterations)

        self.eCount = 1.
       
        self.write = sys.stdout.write
        spaces = ''.join([' ' for i in xrange(self.nElements)])
        loadingBar = ''.join(['[',spaces,']'])
        self.write(loadingBar)
        sys.stdout.flush()
    
        backspacing = ''.join(['\b' for i in xrange(nElements+1)])
        self.write(backspacing)
        
    def progress(self, currentIteration):
        """
        checks the progress of the adjoin process, and advances the progress bar # if necessary
        
        Args:
            currentIteration (int): number of the current iterations of the adjoin process
        """
        
        loopPercentage = currentIteration/(self.nIterations-1.) 
        
        barPercentage = self.eCount/self.nElements
        
        #print "percentage",loopPercentage, barPercentage
        while (barPercentage <= loopPercentage):
            self.write("#")
            sys.stdout.flush()
            self.eCount += 1.
            barPercentage = self.eCount/self.nElements
        
        if loopPercentage == 1:
            self.write("\n")
            sys.stdout.flush()
        

if __name__ == '__main__':
    
    import doctest
    doctest.testmod()
    
        