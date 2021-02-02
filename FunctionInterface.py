from header import *

def Compute_min_max(tstart,tend,directory='./',N=200,field='rho',slice_index=None):
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
        slice_index : integer
            Optional argument to compute min max in a slice of the box
        
    """
    directory = os.path.abspath(directory)
    if(tend<tstart):
        raise ValueError("tend has to be greater than tstart")

    if(slice_index == -1):
        slice_index = int(N/2.)+1

    
    if(tend>tstart):
    #Compute min and max density
        min_data = 1.e50
        max_data = 1.e-50
        for time in range(tstart,tend+1):
            data = read.readsinglefile(directory,time,N,field)
            data = data.reshape(N,N,N)
            if(slice_index):
                data = data[slice_index]            
            min_data = min(min_data,np.min(data))
            max_data = max(max_data,np.max(data))

    else:
        time = tstart
        data = read.readsinglefile(directory,time,N,field)
        data = data.reshape(N,N,N)
        if(slice_index):
            data = data[slice_index] 
        min_data = np.min(data)
        max_data = np.max(data)

    return min_data,max_data


def return_Volavg(directory,tstart,tend,field='rho',N=200,mass_weighted=False):
    directory = os.path.abspath(directory)
    time = tstart

    weighted_sum = 0.0
    sum_weights = 0.0
    average = np.zeros(tend-tstart+1)
    dVol = (4.0*constant.Parsec/(N))**3

    for time in range(tstart,tend+1):
        data = read.readsinglefile(directory,time,N,field).reshape(N,N,N)
        dV = np.full_like(data,dVol)
        data = data
        if(mass_weighted):
            rho = read.readsinglefile(directory,time,N,'rho').reshape(N,N,N)
            weighted_sum = np.sum(data*dV*rho)
            sum_weights = np.sum(dV*rho)
        else:
            weighted_sum = np.sum(data*dV)
            sum_weights = np.sum(dV)

        average[time-tstart] = weighted_sum/sum_weights

    return average


def Make_Movie(directory,basefile='CombinedColumn',tstart=0,tend=100,convert_pdfs=True):
    """
    Function to create movie in mp4 format from set of plots 
    Parameters
        directory : string
            directory where PDF plot outputs present
        basefile: string 
            Base file name of plots
        tend : integer
            last timestep upto which to convert
        convert_pdfs : Boolean
            assumes plots are in pdf. This flag converts them to png
            as ffmpeg does not accept png
    Returns
        None

    """
    # Convert PDFs to PNG
    directory = os.path.abspath(directory) + '/'
    if(convert_pdfs is True):
        print("Converting pdfs to png in directory; {}".format(directory))
        i = 0
        for i in tqdm.trange(tstart,tend+1) : 
            filename = directory + basefile+ "_{:04d}".format(i)
            filename_png = directory + basefile+ "_{:04d}".format(i-tstart)
            os.system("convert -density 400 -background white -alpha remove {}.pdf".format(filename)+" {}.png".format(filename_png))
    #Creating Movie
    print("Creating Movie")
    os.system("ffmpeg -r 10 -i {}{}_%04d.png ".format(directory,basefile)+ 
        "-s:v 2560x1440 -vcodec libx264 -y -pix_fmt yuv420p -loglevel error "+
        "{}animation.mp4".format(directory))
    print("Deleting png files")
    os.system("rm {}*.png".format(directory))



