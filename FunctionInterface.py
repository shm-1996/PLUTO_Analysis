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
            os.system("convert -density 400 {}.pdf".format(filename)+" {}.png".format(filename_png))
    #Creating Movie
    print("Creating Movie")
    os.system("ffmpeg -r 10 -i {}{}_%04d.png ".format(directory,basefile)+ 
        "-s:v 2560x1440 -vcodec libx264 -y -pix_fmt yuv420p -loglevel error "+
        "{}animation.mp4".format(directory))
    print("Deleting png files")
    os.system("rm {}*.png".format(directory))



