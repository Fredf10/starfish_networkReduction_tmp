import psutil, os
        
def memoryUsagePsutil():
    '''
    evaluate the current memory usage from psutils
    
    Returns: memory currently used in MB
    '''
    process = psutil.Process(os.getpid())
    return process.memory_info()[0] / float(2 ** 20)

def getGitHash():
    '''
    This function evaluates the actual branch and leatest git commit hash
    to create an unique version reference for xml files
    
    
    Returns: git hash (str): latest git commit hash
    '''
   
    #branchName   = subprocess.check_output(["git", "describe",'--all'])
    gitHash = str(os.system('git rev-parse --quiet HEAD'))
    print "TODO: DB mHF replace os.sytem with subprocess and reove \enter from hash str!"
    # if 0 no changes otherwise uncommitted changes exist
    dirtyRepos = os.system('git diff --quiet --exit-code')
    
    if dirtyRepos != 0:
        print """WARNING: moduleHelperFunctions.getGitHash(): uncommited changes detected,
         git hash does not correspond to actual code version!"""
    
    return gitHash