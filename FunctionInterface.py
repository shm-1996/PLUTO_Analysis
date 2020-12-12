from header import *

def Compute_min_max(tstart,tend,directory='./',N=200,field='rho'):
    """
    Routine to compute min/max of data
        tstart: integer 
            directory in which files reside
        tend : integer 
            timestep of file
        directory : string
            directory where files reside
        N : integer
            Resolution
        field: string
            Data field to read in
        
    """
    directory = os.path.abspath(directory)
    if(tend<tstart):
        raise ValueError("tend has to be greater than tstart")
    
    if(tend>tstart):
    #Compute min and max density
        min_data = 1.e50
        max_data = 1.e-50
        for time in range(tstart,tend+1):
            data = np.abs(read.readsinglefile(directory,time,N,field))
            min_data = min(min_data,np.min(data))
            max_data = max(max_data,np.max(data))

    else:
        time = tstart
        data = np.abs(read.readsinglefile(directory,time,N,field))
        min_data = np.min(data)
        max_data = np.max(data)

    return min_data,max_data



