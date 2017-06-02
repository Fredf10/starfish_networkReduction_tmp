import os,sys,shutil

cur = os.path.dirname( os.path.realpath( __file__ ) )

  

for directory in os.listdir(cur):
    try:
        print(directory)
        for root, dirs, files in os.walk(''.join([directory,'/SolutionData'])):
            for file in files:
                if '.hdf5' in file:
                    filename = file.split('.')[0]
                    # get folder name and check if exists already otherwise create it
                    solDataFolderName = ''.join(['SolutionData_',filename.split('_')[-1]])
                    solDataFolderPath = '/'.join([cur,directory,solDataFolderName])
                    
                    if not os.path.isdir(solDataFolderPath):
                        print(" create folder for {} and moved data into it".format(solDataFolderName))
                        os.mkdir(solDataFolderPath)
                        # move data into it
                        shutil.move(''.join([cur,'/',directory,'/SolutionData/',file]),solDataFolderPath)
                        shutil.move(''.join([cur,'/',directory,'/SolutionData/',filename,'.xml']),solDataFolderPath)
    except: print("{} nothing to do".format(directory))
            
