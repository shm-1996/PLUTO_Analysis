///
///  Fourier spectra MPI version (single-precision) for FlashUG
///
///  written by Christoph Federrath, 2012-2016
///  Edited for PLUTO Output(Double-precision) by Shyam Harimohan Menon (May 2018) 

///  USAGE : ./powspectrum Outputfilesdirectory Timestepno NX NY NZ
///  Example : ./powspectrum /home/Desktop/outputs 12 200 200 200 will run it for a 200x200x200 simulation for the 12th timestep ie .0012.dbl . 

#include "mpi.h" /// MPI lib
#include <iostream>
#include <iomanip> /// for io manipulations
#include <sstream> /// stringstream
#include <fstream> /// for filestream operation
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fftw3-mpi.h> /// Fast Fourier Transforms MPI version

// constants
int NDIM = 3;
using namespace std;
enum {X, Y, Z};
static const bool Debug = false;
static const int FAILURE = 0;
static const int MAX_NUM_BINS = 10048;
static const double pi = 3.14159265358979323846;
static const bool density_threshold = false ;
// MPI stuff
int MyPE = 0, NPE = 1;

// for FFTW
fftw_complex *fft_data_x, *fft_data_y, *fft_data_z;
fftw_plan fft_plan_x, fft_plan_y, fft_plan_z;

// for output
vector<string> OutputFileHeader;
vector< vector<double> > WriteOutTable;

/// forward function
//vector<int> GetMetaData(const string inputfile, const string datasetname);
//vector<int> GetDimensions(const string inputfile, const string datasetname);
vector<int> InitFFTW(const vector<int> nCells);

//void SwapMemOrder(double * const data, const vector<int> N);
vector<int> InitFFTW(const vector<int> nCells);
void AssignDataToFFTWContainer(const string type, const vector<int> N, const long ntot_local, double * const dens, 
	const double * const velx, const double * const vely, const double * const velz);
void ComputeSpectrum(const vector<int> Dim, const vector<int> MyInds, const bool decomposition);
void WriteOutAnalysedData(const string OutputFilename);
void Normalize(double * const data_array, const long n, const double norm);
double Mean(const double * const data, const long size);
void SwapMemOrder(double * const data, const vector<int> N);
void Window(double* const data_array, const int nx, const int ny, const int nz);
void ReadData(const string inputfile,  double * data_ptr, vector<int> N);

/// --------
///   MAIN
/// --------
int main(int argc, char * argv[])
{

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NPE);
	MPI_Comm_rank(MPI_COMM_WORLD, &MyPE);
	if (MyPE==0) cout<<"=== powspectrum === MPI num procs: "<<NPE<<endl;
	if ((argc < 6)||(argc > 6)) {
	  if (MyPE==0) cout << "Usage: powspectrum INPUTFILEDIRECTORY TIMESTEPNO NX NY NZ" << endl;
	  MPI_Finalize(); exit(FAILURE);
	}
	string inputfiledirectory = argv[1];
	int no = atoi(argv[2]);
	//For converting fileno to desired format say 0056.dbl
	ostringstream fileno;
	fileno<<setfill('0')<<setw(4)<<no;
	 // Whether to use density threshold for calculating the power spectra
	//For time
	long starttime = time(NULL);
	
	string inputfile; // The file name to pass

	//Getting the dataset dimensions TODO: Get them from grid.out file	
	static const int arr[] = {atoi(argv[3]),atoi(argv[4]),atoi(argv[5])};
	
	vector<int> N(arr,arr+sizeof(arr)/sizeof(arr[0]));
	///////////////////////////////////////////////////////////////////

	
	

	cout<<"Looking for files in directory:"<<inputfile<<"With dimensions :  "<<N[X]<<", "<<N[Y]<<", "<<N[Z];
	/// get the dataset dimensions
	//vector<int> N = GetDimensions(inputfile, "velx");
	
	// signal dimensionality of the imulation
	if (N[Z]==1) NDIM = 2;

	/// allocate FFTW containers and create FTTW plan
	vector<int> MyInds = InitFFTW(N);
	long ntot_local = MyInds[1]*N[Y]*N[Z];
	if (Debug) cout<<"["<<MyPE<<"] MyInds: "<<MyInds[0]<<" "<<MyInds[1]<<" ntot_local="<<ntot_local<<endl;

    /// parallelisation / decomposition check
    int wrong_decompostion = 0, wrong_decompostion_red = 0;
    if (MyInds[1] != 0) {
        if (N[X] % MyInds[1] != 0) wrong_decompostion = 1;
    } else {
        wrong_decompostion = 1;
    }
    MPI_Allreduce(&wrong_decompostion, &wrong_decompostion_red, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (wrong_decompostion_red > 0) {
        if (MyPE==0) cout<<"Error: Number of cores is not multiple of N[X]."<<endl;
        MPI_Finalize();
        return 0;
    }

	/// allocate arrays
	double *dens = new double[ntot_local];
	double *velx = new double[ntot_local];
	double *vely = new double[ntot_local];
	double *velz = 0;
    if (NDIM==3) {
		velz = new double[ntot_local];
	}

	/// read data
	if (MyPE==0) cout<<"start reading data from disk..."<<endl;
	inputfile = inputfiledirectory+"rho."+fileno.str()+".dbl";
	ReadData(inputfile, dens, N); if (MyPE==0) cout<<"dens read."<<endl;
	cout<<dens[0]<<dens[1]<<endl;
	inputfile = inputfiledirectory+"vx1."+fileno.str()+".dbl";
	ReadData(inputfile,velx, N); if (MyPE==0) cout<<"velx read."<<endl;
	inputfile = inputfiledirectory+"vx2."+fileno.str()+".dbl";
	ReadData(inputfile,vely, N); if (MyPE==0) cout<<"vely read."<<endl;
    if (NDIM==3) {
		inputfile = inputfiledirectory+"vx3."+fileno.str()+".dbl";
		ReadData(inputfile,velz, N); if (MyPE==0) cout<<"velz read."<<endl;
	}

	/// normalization
	double norm_dens = 1.0;
	double norm_vels = 1.0;
	double norm_mags = 1.0;
	if (MyPE==0) {
		cout << "main:  norm_dens = " << norm_dens << endl;
		cout << "main:  norm_vels = " << norm_vels << endl;
		cout << "main:  norm_mags = " << norm_mags << endl;
	}
	if (false) {
		Normalize(dens, ntot_local, norm_dens);
		Normalize(velx, ntot_local, norm_vels);
		Normalize(vely, ntot_local, norm_vels);
    	if (NDIM==3) {
			Normalize(velz, ntot_local, norm_vels);
		}
	}

	long endtime = time(NULL);
	int duration = endtime-starttime, duration_red = 0;
	if (Debug) cout << "["<<MyPE<<"] ****************** Local time for reading data = "<<duration<<"s ******************" << endl;
	MPI_Allreduce(&duration, &duration_red, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (MyPE==0) cout << "****************** Global time for reading data = "<<duration_red<<"s ******************" << endl;	


	string outfilename = "";

	/// vels spectrum
	AssignDataToFFTWContainer("vels", N, ntot_local, dens, velx, vely, velz);	
	ComputeSpectrum(N, MyInds, true); // decomposition
	if (MyPE==0) {
		outfilename = inputfiledirectory+"powspectrum/_spect_vels"+fileno.str()+".dat";
		WriteOutAnalysedData(outfilename);
	}
	/*
	/// rho3 spectrum
	AssignDataToFFTWContainer("rho3", N, ntot_local, dens, velx, vely, velz);	
	ComputeSpectrum(N, MyInds, true); // decomposition
	if (MyPE==0) {
		outfilename = inputfiledirectory+"powspectrum/_spect_rho3"+fileno.str()+".dat";
		WriteOutAnalysedData(outfilename);
	}

	/// rhov spectrum
	AssignDataToFFTWContainer("rhov", N, ntot_local, dens, velx, vely, velz);	
	ComputeSpectrum(N, MyInds, true); // decomposition
	if (MyPE==0) {
		outfilename = inputfiledirectory+"powspectrum/_spect_rhov"+fileno.str()+".dat";
		WriteOutAnalysedData(outfilename);
	}

	/// varrho spectrum
	AssignDataToFFTWContainer("varrho", N, ntot_local, dens, velx, vely, velz);	
	ComputeSpectrum(N, MyInds, false); // no decomposition
	if (MyPE==0) {
		outfilename = inputfiledirectory+"powspectrum/_spect_varrho"+fileno.str()+".dat";
		WriteOutAnalysedData(outfilename);
	}

	/// varlnrho spectrum
	AssignDataToFFTWContainer("varlnrho", N, ntot_local, dens, velx, vely, velz);	
	ComputeSpectrum(N, MyInds, false); // no decomposition
	if (MyPE==0) {
		outfilename = inputfiledirectory+"powspectrum/_spect_varlnrho"+fileno.str()+".dat";
		WriteOutAnalysedData(outfilename);
	}


	/// rho spectrum
	AssignDataToFFTWContainer("rho", N, ntot_local, dens, velx, vely, velz);	
	ComputeSpectrum(N, MyInds, false); // no decomposition
	if (MyPE==0) {
		outfilename = inputfiledirectory+"powspectrum/_spect_rho"+fileno.str()+".dat";
		WriteOutAnalysedData(outfilename);
	}

	/// lnrho spectrum
	AssignDataToFFTWContainer("lnrho", N, ntot_local, dens, velx, vely, velz);	
	ComputeSpectrum(N, MyInds, false); // no decomposition
	if (MyPE==0) {
		outfilename = inputfiledirectory+"powspectrum/_spect_lnrho"+fileno.str()+".dat";
		WriteOutAnalysedData(outfilename);
	}

	*/

	/// rhoweightedv spectrum: Calculates spectrum of (rho/rhoav)(v-vmean) : Use to compare with Gritschneder et al 2009
	AssignDataToFFTWContainer("rhoweightedv", N, ntot_local, dens, velx, vely, velz);	
	ComputeSpectrum(N, MyInds, true); // no decomposition
	if (MyPE==0) {
		outfilename = inputfiledirectory+"powspectrum/_spect_rhoweightedv"+fileno.str()+".dat";
		WriteOutAnalysedData(outfilename);
	}

	
	/// deallocate and clean
	delete [] dens; ///// DEALLOCATE
	delete [] velx; ///// DEALLOCATE
	delete [] vely; ///// DEALLOCATE
    if (NDIM==3) {
		delete [] velz; ///// DEALLOCATE
	}
	
	fftw_free(fft_data_x);
	fftw_free(fft_data_y);
    if (NDIM==3) {
		fftw_free(fft_data_z);
	}

	fftw_destroy_plan(fft_plan_x);
	fftw_destroy_plan(fft_plan_y);
    if (NDIM==3) {
		fftw_destroy_plan(fft_plan_z);
	}
	fftw_mpi_cleanup();

	endtime = time(NULL);
	duration = endtime-starttime; duration_red = 0;
	if (Debug) cout << "["<<MyPE<<"] ****************** Local time to finish = "<<duration<<"s ******************" << endl;
	MPI_Allreduce(&duration, &duration_red, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (MyPE==0) cout << "****************** Global time to finish = "<<duration_red<<"s ******************" << endl;	

	MPI_Finalize();
  	return 0;

} // end main ==========================================================





/// InitFFTW ==========================================================
vector<int> InitFFTW(const vector<int> nCells)
{
	const bool Debug = false;

	const ptrdiff_t N[3] = {nCells[X], nCells[Y], nCells[Z]};
	ptrdiff_t alloc_local = 0, local_n0 = 0, local_0_start = 0;

	fftw_mpi_init();

	// get local data size and allocate
    if (NDIM==3) {
		alloc_local = fftw_mpi_local_size_3d(N[X], N[Y], N[Z], MPI_COMM_WORLD, &local_n0, &local_0_start);
	}
    if (NDIM==2) {
		alloc_local = fftw_mpi_local_size_2d(N[X], N[Y], MPI_COMM_WORLD, &local_n0, &local_0_start);
	}
	/// ALLOCATE
	if (Debug) cout<<"["<<MyPE<<"] Allocating fft_data_x..."<<endl;
	fft_data_x = fftw_alloc_complex(alloc_local);
	if (Debug) cout<<"["<<MyPE<<"] Allocating fft_data_y..."<<endl;
	fft_data_y = fftw_alloc_complex(alloc_local);
    if (NDIM==3) {
		if (Debug) cout<<"["<<MyPE<<"] Allocating fft_data_z..."<<endl;
		fft_data_z = fftw_alloc_complex(alloc_local);
	}
	if (Debug) cout<<"["<<MyPE<<"] ...alloc done."<<endl;

	/// PLAN
    if (NDIM==3) {
		if (Debug) cout<<"["<<MyPE<<"] fft_plan_x..."<<endl;
		fft_plan_x = fftw_mpi_plan_dft_3d(N[X], N[Y], N[Z], fft_data_x, fft_data_x, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
		if (Debug) cout<<"["<<MyPE<<"] fft_plan_y..."<<endl;
		fft_plan_y = fftw_mpi_plan_dft_3d(N[X], N[Y], N[Z], fft_data_y, fft_data_y, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
		if (Debug) cout<<"["<<MyPE<<"] fft_plan_z..."<<endl;
		fft_plan_z = fftw_mpi_plan_dft_3d(N[X], N[Y], N[Z], fft_data_z, fft_data_z, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
	}
    if (NDIM==2) {
		if (Debug) cout<<"["<<MyPE<<"] fft_plan_x..."<<endl;
		fft_plan_x = fftw_mpi_plan_dft_2d(N[X], N[Y], fft_data_x, fft_data_x, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
		if (Debug) cout<<"["<<MyPE<<"] fft_plan_y..."<<endl;
		fft_plan_y = fftw_mpi_plan_dft_2d(N[X], N[Y], fft_data_y, fft_data_y, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
	}
	if (Debug) cout<<"["<<MyPE<<"] ...plans done."<<endl;

	vector<int> ReturnVector(2);
	ReturnVector[0] = local_0_start;
	ReturnVector[1] = local_n0;
	
	return ReturnVector;

} /// ================================================================



void AssignDataToFFTWContainer(const string type, const vector<int> N, const long ntot_local, double * const dens, 
								const double * const velx, const double * const vely, const double * const velz)
{
	const double onethird = 1./3.;
/*	
	if (type=="rho") {
		for (long n=0; n<ntot_local; n++) {
			fft_data_x[n][0] = sqrt(dens[n]); /// Real part
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = 0.; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = 0.; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}
	if (type=="lnrho") {
		double mean_dens = Mean(dens, ntot_local);
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean dens = "<<mean_dens<<endl;
		for (long n=0; n<ntot_local; n++) {
			fft_data_x[n][0] = log(dens[n]/mean_dens); /// Real part
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = 0.; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = 0.; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}
	if (type=="varrho") {
		double mean_dens = Mean(dens, ntot_local);
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean dens = "<<mean_dens<<endl;
		for (long n=0; n<ntot_local; n++) {
			fft_data_x[n][0] = dens[n]-mean_dens; /// Real part
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = 0.; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = 0.; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}
	if (type=="varlnrho") {
		for (long n=0; n<ntot_local; n++) dens[n] = log(dens[n]);
		double mean_lndens = Mean(dens, ntot_local);
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean log(dens) = "<<mean_lndens<<endl;
		for (long n=0; n<ntot_local; n++) dens[n] = exp(dens[n]);
		for (long n=0; n<ntot_local; n++) {
			fft_data_x[n][0] = log(dens[n])-mean_lndens; /// Real part
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = 0.; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = 0.; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}
	 
	
	if (type=="sqrtrho") {
		for (long n=0; n<ntot_local; n++) {
			fft_data_x[n][0] = sqrt(dens[n])*velx[n]; /// Real part
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = sqrt(dens[n])*vely[n]; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = sqrt(dens[n])*velz[n]; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}
	if (type=="rho3") {
		for (long n=0; n<ntot_local; n++) {
			fft_data_x[n][0] = pow(dens[n],onethird)*velx[n]; /// Real part
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = pow(dens[n],onethird)*vely[n]; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = pow(dens[n],onethird)*velz[n]; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}
	if (type=="rhov") {
		for (long n=0; n<ntot_local; n++) {
			fft_data_x[n][0] = dens[n]*velx[n]; /// Real part
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = dens[n]*vely[n]; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = dens[n]*velz[n]; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}

	*/
	if (type=="rhoweightedv") {
		
		
		double mean_velx = Mean(velx, ntot_local);
		double mean_vely = Mean(vely, ntot_local);
		double mean_velz = Mean(velz, ntot_local);

		if(density_threshold==true)
		{
			//TO take into consideration range of gas density
			for(long n=0;n<ntot_local;n++) 
			{
				if(dens[n]<0.016714)    // Taking into account gas of no densities only between 1e2 and 1e4 cm^-3 
					dens[n] = 0.0 ;
				else
					dens[n] = dens[n] ;
				
			}
		}


// To calculate weighted mean of velcity
		double *weightedvelx = new double[ntot_local];
		double *weightedvely = new double[ntot_local];
		double *weightedvelz = new double[ntot_local];
  		double mean_dens = Mean(dens, ntot_local);
		for(long n=0;n<ntot_local;n++) 
		{

			weightedvelx[n] = pow((dens[n]/mean_dens),0.5)*velx[n];
			weightedvely[n] = pow((dens[n]/mean_dens),0.5)*vely[n];
			weightedvelz[n] = pow((dens[n]/mean_dens),0.5)*velz[n];
		}


		double mean_weightedvelx = Mean(weightedvelx, ntot_local);
		double mean_weightedvely = Mean(weightedvely, ntot_local);
		double mean_weightedvelz = Mean(weightedvelz, ntot_local);
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean dens = "<<mean_dens<<endl;
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean velx = "<<mean_velx<<endl;
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean vely = "<<mean_vely<<endl;
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean velz = "<<mean_velz<<endl ;
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean Weightedvelx = "<<mean_weightedvelx<<endl ;
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean Weightedvely = "<<mean_weightedvely<<endl ;
		if (MyPE==0) cout<<"AssignDataToFFTWContainer: mean Weightedvelz = "<<mean_weightedvelz<<endl ;	
		int count =0;
		for (long n=0; n<ntot_local; n++) { 
			if(dens[n]>0.0)fft_data_x[n][0] = weightedvelx[n] - mean_weightedvelx ; /// Real part
			else	        fft_data_x[n][0] = 0.0;
			//cout<<"Point No:"<<n<<"\t Velocity = "<<weightedvelx[n]<<endl; //TODO:Debugging comment out or remove if I have forgoten to do so 
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = weightedvely[n] - mean_weightedvely; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = weightedvelz[n] - mean_weightedvelz; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}


	if (type=="vels") {
		double mean_velx = Mean(velx,ntot_local);
		double mean_vely = Mean(vely,ntot_local);
		double mean_velz = Mean(velz,ntot_local);		
		for (long n=0; n<ntot_local; n++) {
			fft_data_x[n][0] = velx[n]-mean_velx; /// Real part
			fft_data_x[n][1] = 0.; /// Imaginary part
			fft_data_y[n][0] = vely[n]-mean_vely; /// Real part
			fft_data_y[n][1] = 0.; /// Imaginary part
			if (NDIM==3) {
				fft_data_z[n][0] = velz[n]-mean_velz; /// Real part
				fft_data_z[n][1] = 0.; /// Imaginary part
			}
		}
	}


	
	/// error check
	double sum = 0.0, sum_red = 0.0;
	double ntot = (double)(N[X])*(double)(N[Y])*(double)(N[Z]);
	if (NDIM==3)
		for (long n = 0; n < ntot_local; n++)
			sum += fft_data_x[n][0]*fft_data_x[n][0]+fft_data_y[n][0]*fft_data_y[n][0]+fft_data_z[n][0]*fft_data_z[n][0];
	if (NDIM==2)
		for (long n = 0; n < ntot_local; n++)
			sum += fft_data_x[n][0]*fft_data_x[n][0]+fft_data_y[n][0]*fft_data_y[n][0];
	if (Debug) cout << "["<<MyPE<<"] Local sum in physical space = " << sum/ntot << endl;
	MPI_Allreduce(&sum, &sum_red, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if (MyPE==0) cout << "Global sum in physical space ("<<type<<") = " << sum_red/ntot << endl;
	
} /// ==============================================================


/** --------------------- ComputeSpectrum ----------------------------
 **  computes total, transversal, longitudinal spectrum functions
 ** ------------------------------------------------------------------ */
void ComputeSpectrum(const vector<int> Dim, const vector<int> MyInds, const bool decomposition)
{
	const bool Debug = false;
	long starttime = time(NULL);
	
	if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: entering." << endl;
	
	if (NDIM==3) {
		if ((Dim[X]!=Dim[Y])||(Dim[X]!=Dim[Z])||(Dim[Y]!=Dim[Z])) {
	    	cout << "Spectra can only be obtained from cubic datasets (Nx=Ny=Nz)." << endl;
	    	exit(FAILURE);
		}
	}
	if (NDIM==2) {
		if (Dim[X]!=Dim[Y]) {
	    	cout << "Spectra can only be obtained from quadratic datasets (Nx=Ny)." << endl;
	    	exit(FAILURE);
		}
	}

	/////////// EXECUTE PLAN
	if (decomposition) {
		fftw_execute(fft_plan_x);
		fftw_execute(fft_plan_y);
		if (NDIM==3) {
			fftw_execute(fft_plan_z);
		}
	} else {
		fftw_execute(fft_plan_x);		
	}
	
	/// general constants
	const int N = Dim[X]; /// assume a cubic (square in 2D) box !
	long LocalNumberOfDataPoints = 0;
	long TotalNumberOfDataPoints = 0;
	if (NDIM==3) {
		LocalNumberOfDataPoints = N*N*MyInds[1];	
		TotalNumberOfDataPoints = N*N*N;
	}
	if (NDIM==2) {
		LocalNumberOfDataPoints = N*MyInds[1];	
		TotalNumberOfDataPoints = N*N;
	}
	const double TotalNumberOfDataPointsDouble = (double)(TotalNumberOfDataPoints);
	const double sqrt_TotalNumberOfDataPoints = sqrt((double)(TotalNumberOfDataPoints));

	/// allocate containers
	if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: Allocating energy_spect..." << endl;
	double       *energy_spect     = new double[LocalNumberOfDataPoints];
	if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: Allocating energy_lt_spect..." << endl;
	double       *energy_lt_spect  = new double[LocalNumberOfDataPoints];
	if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: ...allocating done." << endl;	
	for (long n = 0; n < LocalNumberOfDataPoints; n++) {
		energy_spect[n]    = 0.0;
		energy_lt_spect[n] = 0.0;
	}

	/// FFTW normalization
	if (decomposition) 
	{
		if (NDIM==3) {
		  for (long n = 0; n < LocalNumberOfDataPoints; n++)
		  {
			fft_data_x[n][0] /= sqrt_TotalNumberOfDataPoints;
			fft_data_x[n][1] /= sqrt_TotalNumberOfDataPoints;
			fft_data_y[n][0] /= sqrt_TotalNumberOfDataPoints;
			fft_data_y[n][1] /= sqrt_TotalNumberOfDataPoints;
			fft_data_z[n][0] /= sqrt_TotalNumberOfDataPoints;
			fft_data_z[n][1] /= sqrt_TotalNumberOfDataPoints;
			energy_spect[n] += ( fft_data_x[n][0]*fft_data_x[n][0]+fft_data_x[n][1]*fft_data_x[n][1] +
		                         fft_data_y[n][0]*fft_data_y[n][0]+fft_data_y[n][1]*fft_data_y[n][1] +
		                         fft_data_z[n][0]*fft_data_z[n][0]+fft_data_z[n][1]*fft_data_z[n][1]  ) / TotalNumberOfDataPointsDouble;
		  }
		}
		if (NDIM==2) {
		  for (long n = 0; n < LocalNumberOfDataPoints; n++)
		  {
			fft_data_x[n][0] /= sqrt_TotalNumberOfDataPoints;
			fft_data_x[n][1] /= sqrt_TotalNumberOfDataPoints;
			fft_data_y[n][0] /= sqrt_TotalNumberOfDataPoints;
			fft_data_y[n][1] /= sqrt_TotalNumberOfDataPoints;
			energy_spect[n] += ( fft_data_x[n][0]*fft_data_x[n][0]+fft_data_x[n][1]*fft_data_x[n][1] +
		                         fft_data_y[n][0]*fft_data_y[n][0]+fft_data_y[n][1]*fft_data_y[n][1] ) / TotalNumberOfDataPointsDouble;
		  }
		}
	}
	else
	{
		for (long n = 0; n < LocalNumberOfDataPoints; n++)
		{
		    fft_data_x[n][0] /= sqrt_TotalNumberOfDataPoints;
		    fft_data_x[n][1] /= sqrt_TotalNumberOfDataPoints;
		    energy_spect[n] += (fft_data_x[n][0]*fft_data_x[n][0]+fft_data_x[n][1]*fft_data_x[n][1]) / TotalNumberOfDataPointsDouble;
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	
	double tot_energy_spect = 0.0, tot_energy_spect_red = 0.0;
	for (long n = 0; n < LocalNumberOfDataPoints; n++)
	    tot_energy_spect += energy_spect[n];
	if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: Local sum in spectral space = " << tot_energy_spect << endl;
	MPI_Allreduce(&tot_energy_spect, &tot_energy_spect_red, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tot_energy_spect = tot_energy_spect_red;
	if (MyPE==0) cout << "ComputeSpectrum: Global sum in spectral space = " << tot_energy_spect << endl;
        
        




       //TODO : Check this part of code. Issue seems most likely here.
	/// compute longitudinal spectrum (remember how FFTW sorts the k-values)
	double tot_energy_lt_spect = 0.0, tot_energy_lt_spect_red = 0.0;
	double dec_lt_0 = 0.0; double dec_lt_1 = 0.0;
	if (decomposition)
	{
	    int k1 = 0; int k2 = 0; int k3 = 0;
        for (int j = MyInds[0]; j < MyInds[0]+MyInds[1]; j++) // parallelized bit (only loop local part)
        {
          if (j <= Dim[X]/2.) k1 = j; else k1 = j-Dim[X];
	      for (int l = 0; l < Dim[Y]; l++)
	      {
	        if (l <= Dim[Y]/2.) k2 = l; else k2 = l-Dim[Y];
		    for (int m = 0; m < Dim[Z]; m++)
		    {
		      if (m <= Dim[Z]/2.) k3 = m; else k3 = m-Dim[Z];
	
              double k_sqr_index = k1*k1 + k2*k2 + k3*k3;
			  long index = (j-MyInds[0])*Dim[Y]*Dim[Z] + l*Dim[Z] + m; // row-major
			  if (NDIM==3) {
				dec_lt_0 = k1*fft_data_x[index][0] + k2*fft_data_y[index][0] + k3*fft_data_z[index][0];
				dec_lt_1 = k1*fft_data_x[index][1] + k2*fft_data_y[index][1] + k3*fft_data_z[index][1];
			  }
			  if (NDIM==2) {
				dec_lt_0 = k1*fft_data_x[index][0] + k2*fft_data_y[index][0];
				dec_lt_1 = k1*fft_data_x[index][1] + k2*fft_data_y[index][1];
			  }
		  	  if (k_sqr_index > 0)
		        energy_lt_spect[index] = ((dec_lt_0*dec_lt_0+dec_lt_1*dec_lt_1)/k_sqr_index)/TotalNumberOfDataPointsDouble;
	        }
	      }
	    }

		for (long n = 0; n < LocalNumberOfDataPoints; n++) tot_energy_lt_spect += energy_lt_spect[n];
		if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: Local sum of longitudinal part in spectral space = " << tot_energy_lt_spect << endl;
		MPI_Allreduce(&tot_energy_lt_spect, &tot_energy_lt_spect_red, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tot_energy_lt_spect = tot_energy_lt_spect_red;
		if (MyPE==0) cout << "ComputeSpectrum: Global sum of longitudinal part in spectral space = " << tot_energy_lt_spect << endl;

	} // decomposition





	/// compute the maximum k and construct the spect_grid i.e. the k_axis as well as a staggered grid
	const int    k_cut = N/2;
	const double log_increment_k = 0.01;
	const double increase_k_factor = pow(10.0, 2.0*log_increment_k);
	double       spect_grid_stag[MAX_NUM_BINS]; spect_grid_stag[0] = 0.0;
	double       spect_grid     [MAX_NUM_BINS]; spect_grid[0]      = 0.0;
	double       k_sqr = 1.0;
	int          bin_index = 0;
	while (k_sqr <= (k_cut+1.0)*(k_cut+1.0))
	{
	    bin_index++;
	    if (k_sqr <= (k_cut+1.0)*(k_cut+1.0))
		  k_sqr = bin_index*bin_index;
	    else
		  k_sqr = k_sqr * increase_k_factor;
	    spect_grid_stag[bin_index] = k_sqr;
	    if (bin_index >= MAX_NUM_BINS)
	    {
		  cout << "["<<MyPE<<"] ComputeSpectrum: ERROR. Number of spectral bins exceeds maximum number." << endl;
		  exit(FAILURE);
	    }
	}
	const int numbins = bin_index;

	/// construct the spectral grid
	for (int bin = 1; bin <= numbins; bin++)
	    spect_grid[bin-1] = pow(sqrt(spect_grid_stag[bin-1])+(sqrt(spect_grid_stag[bin])-sqrt(spect_grid_stag[bin-1]))/2.0, 2.0);

	/// calculate spectral densities
	if (MyPE==0 && Debug) cout << "ComputeSpectrum: Calculating spectral densities ..." << endl;
	double spect_binsum            [3][numbins]; /// contains longit., transv., total spectral densities
	double spect_binsum_sqr        [3][numbins]; /// this is used to compute the RMS and finally the sigma
	double sigma_spect_binsum      [3][numbins]; /// contains sigmas
	double spect_funct             [3][numbins]; /// contains longit., transv., total spectrum functions
	double sigma_spect_funct       [3][numbins]; /// contains sigmas
	double comp_lt_spect_funct        [numbins]; /// Kolmogorov compensated spectrum function
	double sigma_comp_lt_spect_funct  [numbins]; /// contains sigmas
	double comp_trsv_spect_funct      [numbins]; /// Kolmogorov compensated spectrum function
	double sigma_comp_trsv_spect_funct[numbins]; /// contains sigmas
	double diss_spect_funct           [numbins]; /// dissipative spectrum function
	double sigma_diss_spect_funct     [numbins]; /// contains sigmas
	double spect_binsum_lin           [numbins]; /// contains spectral densities of the non-decomposed dataset
	double spect_binsum_lin_sqr       [numbins]; /// this is used to compute the RMS and finally the sigma
	double sigma_spect_binsum_lin     [numbins]; /// contains sigmas
	double spect_funct_lin            [numbins]; /// contains lin spectrum function
	double sigma_spect_funct_lin      [numbins]; /// contains sigmas
	long   n_cells                    [numbins]; /// the number of cells inside a spherical shell in k-space
	long   n_count = 0;

	for (int bin = 0; bin < numbins; bin++) /// set containers to zero
	{
	  if (decomposition)
	  {
	    for (int type = 0; type < 3; type++) // type means long, trans, total
	    {
		  spect_binsum      [type][bin] = 0.0;
		  spect_binsum_sqr  [type][bin] = 0.0;
		  sigma_spect_binsum[type][bin] = 0.0;
		  spect_funct       [type][bin] = 0.0;
		  sigma_spect_funct [type][bin] = 0.0;
	    }
	    comp_lt_spect_funct        [bin] = 0.0;
	    sigma_comp_lt_spect_funct  [bin] = 0.0;
	    comp_trsv_spect_funct      [bin] = 0.0;
	    sigma_comp_trsv_spect_funct[bin] = 0.0;
	    diss_spect_funct           [bin] = 0.0;
	    sigma_diss_spect_funct     [bin] = 0.0;
	  }
	  else // no decomposition
	  {
		spect_binsum_lin      [bin] = 0.0;
		spect_binsum_lin_sqr  [bin] = 0.0;
		sigma_spect_binsum_lin[bin] = 0.0;
		spect_funct_lin       [bin] = 0.0;
		sigma_spect_funct_lin [bin] = 0.0;
	  }
	  n_cells[bin] = 0;
	}

	int k1 = 0; int k2 = 0; int k3 = 0; // these are the time consuming loops (start optimization here)
    for (int j = MyInds[0]; j < MyInds[0]+MyInds[1]; j++) // the parallel bit
    {
  	  if (j <= Dim[X]/2.) k1 = j; else k1 = j-Dim[X];
	  for (int l = 0; l < Dim[Y]; l++)
	  {
	    if (l <= Dim[Y]/2.) k2 = l; else k2 = l-Dim[Y];
	    for (int m = 0; m < Dim[Z]; m++)
	    {
	      if (m <= Dim[Z]/2.) k3 = m; else k3 = m-Dim[Z];
	
		  long k_sqr_index = k1*k1 + k2*k2 + k3*k3;
		  int interval_l = 0; int interval_r = numbins-1; int bin_id = 0;
		  while ((interval_r - interval_l) > 1) /// nested intervals
		  {
		    bin_id = interval_l + (interval_r - interval_l)/2;
		    if (spect_grid[bin_id] > k_sqr_index) interval_r = bin_id;
		    else                                  interval_l = bin_id;
		  }
		  bin_id = interval_r;
		  if ((bin_id <= 0) || (bin_id > numbins-1))
		  {
		    cout << "["<<MyPE<<"] ComputeSpectrum: ERROR. illegal bin index." << endl;
		    exit(FAILURE);
		  }
		  long index = (j-MyInds[0])*Dim[Y]*Dim[Z] + l*Dim[Z] + m; // row-major
		  {
		    if (decomposition)
		    {
			  double energy_trsv_spect = energy_spect[index] - energy_lt_spect[index];
			  spect_binsum    [0][bin_id] += energy_lt_spect[index];
			  spect_binsum    [1][bin_id] += energy_trsv_spect;
			  spect_binsum    [2][bin_id] += energy_spect[index];
			  spect_binsum_sqr[0][bin_id] += energy_lt_spect[index]*energy_lt_spect[index];
			  spect_binsum_sqr[1][bin_id] += energy_trsv_spect*energy_trsv_spect;
			  spect_binsum_sqr[2][bin_id] += energy_spect[index]*energy_spect[index];
		    }
		    else // no decomposition
		    {
			  spect_binsum_lin    [bin_id] += energy_spect[index];
			  spect_binsum_lin_sqr[bin_id] += energy_spect[index]*energy_spect[index];
		    }
		    n_cells[bin_id]++;
		    n_count++;
		  }
	    } // j
	  } // l
	} // m

	/// resum the number of cells and total energy in k-space for error checking
	long n_cells_tot = 0, n_cells_tot_red = 0;
	tot_energy_spect = 0.0;
	for (int bin = 0; bin < numbins; bin++)
	{
	    n_cells_tot += n_cells[bin];
	    if (decomposition)  tot_energy_spect += spect_binsum[2][bin];
	    if (!decomposition) tot_energy_spect += spect_binsum_lin[bin];
	}
	if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: Local ReSummed total number of cells   = " << n_cells_tot << endl;
	MPI_Allreduce(&n_cells_tot, &n_cells_tot_red, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	n_cells_tot = n_cells_tot_red;
	if (MyPE==0) cout << "ComputeSpectrum: Global ReSummed total number of cells = " << n_cells_tot << endl;
	if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: Local ReSummed total in spectral space = " << tot_energy_spect << endl;
	MPI_Allreduce(&tot_energy_spect, &tot_energy_spect_red, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tot_energy_spect = tot_energy_spect_red;
	if (MyPE==0) cout << "ComputeSpectrum: Global ReSummed energy in spectral space = " << tot_energy_spect << endl;

	/// MPI Allreduce of the bin containers
	int *tmp_i_red = 0; double *tmp_d_red = 0;
	int *tmp_i = 0; double *tmp_d = 0;
	
	tmp_i_red = new int[numbins]; tmp_i = new int[numbins];
	for (int n=0; n<numbins; n++) tmp_i[n] = n_cells[n];
	MPI_Allreduce(tmp_i, tmp_i_red, numbins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	for (int n=0; n<numbins; n++) n_cells[n]=tmp_i_red[n];
	delete [] tmp_i; delete [] tmp_i_red;
 
	tmp_d_red = new double[3*numbins]; tmp_d = new double[3*numbins];
   	for (int n=0; n<numbins; n++) for (int dir=0; dir<3; dir++) tmp_d[3*n+dir] = spect_binsum[dir][n];
	MPI_Allreduce(tmp_d, tmp_d_red, 3*numbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for (int n=0; n<numbins; n++) for (int dir=0; dir<3; dir++) spect_binsum[dir][n]=tmp_d_red[3*n+dir];
   	for (int n=0; n<numbins; n++) for (int dir=0; dir<3; dir++) tmp_d[3*n+dir] = spect_binsum_sqr[dir][n];
	MPI_Allreduce(tmp_d, tmp_d_red, 3*numbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for (int n=0; n<numbins; n++) for (int dir=0; dir<3; dir++) spect_binsum_sqr[dir][n]=tmp_d_red[3*n+dir];
	delete [] tmp_d; delete [] tmp_d_red;

	tmp_d_red = new double[numbins]; tmp_d = new double[numbins];
   	for (int n=0; n<numbins; n++) tmp_d[n] = spect_binsum_lin[n];
	MPI_Allreduce(tmp_d, tmp_d_red, numbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for (int n=0; n<numbins; n++) spect_binsum_lin[n]=tmp_d_red[n];
   	for (int n=0; n<numbins; n++) tmp_d[n] = spect_binsum_lin_sqr[n];
	MPI_Allreduce(tmp_d, tmp_d_red, numbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for (int n=0; n<numbins; n++) spect_binsum_lin_sqr[n]=tmp_d_red[n];
	delete [] tmp_d; delete [] tmp_d_red;
		
	/// write out (MASTER CPU only)
	if (MyPE==0)
	{
		/// calculate spectral densities and functions (normalization)
		for (int bin = 0; bin < numbins; bin++)
		{
	  		if (n_cells[bin] > 0)
	  		{
	    		if (decomposition)
	    		{
	      			for (int dir = 0; dir < 3; dir++) /// long., transv., total
	      			{
						spect_binsum      [dir][bin] /= static_cast<double>(n_cells[bin]);
						spect_binsum_sqr  [dir][bin] /= static_cast<double>(n_cells[bin]);
						sigma_spect_binsum[dir][bin]  = sqrt(spect_binsum_sqr[dir][bin] - spect_binsum[dir][bin]*spect_binsum[dir][bin]);
						if (NDIM==3) {
						  spect_funct       [dir][bin]  = 4*pi*spect_grid_stag[bin]*spect_binsum      [dir][bin];
						  sigma_spect_funct [dir][bin]  = 4*pi*spect_grid_stag[bin]*sigma_spect_binsum[dir][bin];
						}
						if (NDIM==2) {
						  spect_funct       [dir][bin]  = 2*pi*sqrt(spect_grid_stag[bin])*spect_binsum      [dir][bin];
						  sigma_spect_funct [dir][bin]  = 2*pi*sqrt(spect_grid_stag[bin])*sigma_spect_binsum[dir][bin];
						}
	      			}
	      			comp_lt_spect_funct        [bin] = pow(spect_grid_stag[bin], 2.0/2.0) * spect_funct      [0][bin];
	      			sigma_comp_lt_spect_funct  [bin] = pow(spect_grid_stag[bin], 2.0/2.0) * sigma_spect_funct[0][bin];
	      			comp_trsv_spect_funct      [bin] = pow(spect_grid_stag[bin], 5.0/6.0) * spect_funct      [1][bin];
	      			sigma_comp_trsv_spect_funct[bin] = pow(spect_grid_stag[bin], 5.0/6.0) * sigma_spect_funct[1][bin];
	      			diss_spect_funct           [bin] = spect_grid_stag[bin] * spect_funct      [2][bin];
	      			sigma_diss_spect_funct     [bin] = spect_grid_stag[bin] * sigma_spect_funct[2][bin];
	    		}
	    		else // no decomposition
	    		{
					spect_binsum_lin      [bin] /= static_cast<double>(n_cells[bin]);
					spect_binsum_lin_sqr  [bin] /= static_cast<double>(n_cells[bin]);
					sigma_spect_binsum_lin[bin]  = sqrt(spect_binsum_lin_sqr[bin] - spect_binsum_lin[bin]*spect_binsum_lin[bin]);
					if (NDIM==3) {
					  spect_funct_lin       [bin]  = 4*pi*spect_grid_stag[bin]*spect_binsum_lin[bin];
					  sigma_spect_funct_lin [bin]  = 4*pi*spect_grid_stag[bin]*sigma_spect_binsum_lin[bin];
					}
					if (NDIM==2) {
					  spect_funct_lin       [bin]  = 2*pi*sqrt(spect_grid_stag[bin])*spect_binsum_lin[bin];
					  sigma_spect_funct_lin [bin]  = 2*pi*sqrt(spect_grid_stag[bin])*sigma_spect_binsum_lin[bin];
					}
	    		}
	  		}
		}

		/// prepare OutputFileHeader
		OutputFileHeader.resize(0);
		stringstream dummystream;
		if (decomposition)
		{
	    	dummystream.precision(8);
	    	dummystream << "E_tot = " << endl;
	    	dummystream << scientific << tot_energy_spect << endl;
	    	dummystream << "E_lgt = " << endl;
	    	dummystream << scientific << tot_energy_lt_spect << endl << endl;
	    	dummystream << setw(30) << left << "#00_BinIndex";
	    	dummystream << setw(30) << left << "#01_KStag"              << setw(30) << left << "#02_K";
	    	dummystream << setw(30) << left << "#03_DK"                 << setw(30) << left << "#04_NCells";
	    	dummystream << setw(30) << left << "#05_SpectDensLgt"       << setw(30) << left << "#06_SpectDensLgtSigma";
	    	dummystream << setw(30) << left << "#07_SpectDensTrv"       << setw(30) << left << "#08_SpectDensTrvSigma";
	    	dummystream << setw(30) << left << "#09_SpectDensTot"       << setw(30) << left << "#10_SpectDensTotSigma";
	    	dummystream << setw(30) << left << "#11_SpectFunctLgt"      << setw(30) << left << "#12_SpectFunctLgtSigma";
	    	dummystream << setw(30) << left << "#13_SpectFunctTrv"      << setw(30) << left << "#14_SpectFunctTrvSigma";
	    	dummystream << setw(30) << left << "#15_SpectFunctTot"      << setw(30) << left << "#16_SpectFunctTotSigma";
	    	dummystream << setw(30) << left << "#17_CompSpectFunctLgt"  << setw(30) << left << "#18_CompSpectFunctLgtSigma";
	    	dummystream << setw(30) << left << "#19_CompSpectFunctTrv"  << setw(30) << left << "#20_CompSpectFunctTrvSigma";
	    	dummystream << setw(30) << left << "#21_DissSpectFunct"     << setw(30) << left << "#22_DissSpectFunctSigma";
	    	OutputFileHeader.push_back(dummystream.str()); dummystream.clear(); dummystream.str("");
		}
		else // no decomposition
		{
	    	dummystream << setw(30) << left << "#00_BinIndex";
	    	dummystream << setw(30) << left << "#01_KStag"           << setw(30) << left << "#02_K";
	    	dummystream << setw(30) << left << "#03_DK"              << setw(30) << left << "#04_NCells";
	    	dummystream << setw(30) << left << "#05_SpectDens"       << setw(30) << left << "#06_SpectDensSigma";
	    	dummystream << setw(30) << left << "#07_SpectFunct"      << setw(30) << left << "#08_SpectFunctSigma";
	    	OutputFileHeader.push_back(dummystream.str()); dummystream.clear(); dummystream.str("");
		}

		if (decomposition)
		{
	  		/// resize and fill WriteOutTable
	  		WriteOutTable.resize(numbins-2); /// spectrum output has numbins-2 lines
	  		for (unsigned int i = 0; i < WriteOutTable.size(); i++)
				WriteOutTable[i].resize(23); /// dec energy spectrum output has 23 columns
	  		for (int bin = 1; bin < numbins-1; bin++)
	  		{
	      		int wob = bin-1;
	      		WriteOutTable[wob][ 0] = bin;
	      		WriteOutTable[wob][ 1] = sqrt(spect_grid_stag[bin]); /// k (staggered)
	      		WriteOutTable[wob][ 2] = sqrt(spect_grid     [bin]); /// k
	      		WriteOutTable[wob][ 3] = sqrt(spect_grid[bin])-sqrt(spect_grid[bin-1]); /// delta k
	      		WriteOutTable[wob][ 4] = n_cells                    [bin]; /// the number of cells in bin
	      		WriteOutTable[wob][ 5] = spect_binsum           [0] [bin]; /// longitudinal spectral density
	      		WriteOutTable[wob][ 6] = sigma_spect_binsum     [0] [bin]; /// sigma
	      		WriteOutTable[wob][ 7] = spect_binsum           [1] [bin]; /// transversal spectral density
	      		WriteOutTable[wob][ 8] = sigma_spect_binsum     [1] [bin]; /// sigma
	      		WriteOutTable[wob][ 9] = spect_binsum           [2] [bin]; /// total spectral density
	      		WriteOutTable[wob][10] = sigma_spect_binsum     [2] [bin]; /// sigma
	      		WriteOutTable[wob][11] = spect_funct            [0] [bin]; /// longitudinal spectrum function
	      		WriteOutTable[wob][12] = sigma_spect_funct      [0] [bin]; /// sigma
	      		WriteOutTable[wob][13] = spect_funct            [1] [bin]; /// transversal spectrum function
	      		WriteOutTable[wob][14] = sigma_spect_funct      [1] [bin]; /// sigma
	      		WriteOutTable[wob][15] = spect_funct            [2] [bin]; /// total spectrum function
	      		WriteOutTable[wob][16] = sigma_spect_funct      [2] [bin]; /// sigma
	      		WriteOutTable[wob][17] = comp_lt_spect_funct        [bin]; /// compensated longitudinal spectrum function
	      		WriteOutTable[wob][18] = sigma_comp_lt_spect_funct  [bin]; /// sigma
	      		WriteOutTable[wob][19] = comp_trsv_spect_funct      [bin]; /// compensated tranversal spectrum function
	      		WriteOutTable[wob][20] = sigma_comp_trsv_spect_funct[bin]; /// sigma
	      		WriteOutTable[wob][21] = diss_spect_funct           [bin]; /// dissipative spectrum function
	      		WriteOutTable[wob][22] = sigma_diss_spect_funct     [bin]; /// sigma
	  		}
		}
		else // no decomposition
		{
	  		/// resize and fill WriteOutTable
	  		WriteOutTable.resize(numbins-2); /// spectrum output has numbins-2 lines
	  		for (unsigned int i = 0; i < WriteOutTable.size(); i++)
				WriteOutTable[i].resize(9); /// density spectrum output has 9 columns
	  		for (int bin = 1; bin < numbins-1; bin++)
	  		{
	      		int wob = bin-1;
	      		WriteOutTable[wob][0] = bin;
	      		WriteOutTable[wob][1] = sqrt(spect_grid_stag[bin]); /// k (staggered)
	      		WriteOutTable[wob][2] = sqrt(spect_grid     [bin]); /// k
	      		WriteOutTable[wob][3] = sqrt(spect_grid[bin])-sqrt(spect_grid[bin-1]); /// delta k
	      		WriteOutTable[wob][4] = n_cells                   [bin]; /// the number of cells in bin
	      		WriteOutTable[wob][5] = spect_binsum_lin          [bin]; /// spectral density of non-decomposed dataset
	      		WriteOutTable[wob][6] = sigma_spect_binsum_lin    [bin]; /// sigma
	      		WriteOutTable[wob][7] = spect_funct_lin           [bin]; /// spectrum function of non-decomposed dataset
	      		WriteOutTable[wob][8] = sigma_spect_funct_lin     [bin]; /// sigma
	  		}
		}
		
	} // MyPE==0

	/// clean up
	delete [] energy_spect;
	delete [] energy_lt_spect;

	long endtime = time(NULL);
	int duration = endtime-starttime, duration_red = 0;
	if (Debug) cout << "["<<MyPE<<"] ****************** Local elapsed time for spectrum function computation = "<<duration<<"s ******************" << endl;
	MPI_Allreduce(&duration, &duration_red, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (MyPE==0) cout << "****************** Global elapsed time for spectrum function computation = "<<duration_red<<"s ******************" << endl;
	if (Debug) cout << "["<<MyPE<<"] ComputeSpectrum: exiting." << endl;
} /// =======================================================================


/** --------------------- Normalize -----------------------------------------------
 ** divide by norm
 ** ------------------------------------------------------------------------------- */
void Normalize(double * const data_array, const long n, const double norm)
{
	for (long i = 0; i < n; i++) data_array[i] /= norm;
} /// =======================================================================


/** ----------------------------- Mean -------------------------------
 **  computes the mean of a pointer-array
 ** ------------------------------------------------------------------ */
double Mean(const double * const data, const long size)
{
	long local_size = size;
	long global_size = 0;
	double value = 0.0, value_red = 0.0;
	for (long n = 0; n < local_size; n++)
		value += data[n];
	MPI_Allreduce(&value, &value_red, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_size, &global_size, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	value_red /= static_cast<double>(global_size);
	return value_red;
} /// =======================================================================


/** ------------------------ Window -----------------------------------------------
 ** apply window function to data ("Hann" or "Hanning" window in 3D)
 ** ------------------------------------------------------------------------------- */
void Window(double* const data_array, const int nx, const int ny, const int nz)
{
	if ((nx!=ny)||(nx!=nz))
	{
		cout << "Window: only works for nx=ny=nz." << endl;
		exit(FAILURE);
	}
	const long nxy = nx*ny;
	const double twopi = 2.*pi;
	const double L = (double)(nx);

	for (int k = 0; k < nz; k++)
	{
	  for (int j = 0; j < ny; j++)
	  {
	    for (int i = 0; i < nx; i++)
	    {
		double dx = (double)(i)+0.5-((double)(nx)/2.);
		double dy = (double)(j)+0.5-((double)(ny)/2.);
		double dz = (double)(k)+0.5-((double)(nz)/2.);
		double r = sqrt(dx*dx+dy*dy+dz*dz);
		long index = k*nxy+j*nx+i;
		if (r < L/2.)
			data_array[index] *= 0.5*(1.0+cos((twopi*r)/L));
		else
			data_array[index] = 0.;
	    } //i
	  } //j
	} //k
} /// =======================================================================


/** -------------------- WriteOutAnalysedData ---------------------------------
 **  Writes out a variable table of data and a FileHeader to a specified file
 ** --------------------------------------------------------------------------- */
void WriteOutAnalysedData(const string OutputFilename)
{
	/// open output file
	ofstream Outputfile(OutputFilename.c_str());
	
	/// check for file
	if (!Outputfile)
	{
	    cout << "WriteOutAnalysedData:  File system error. Could not create '" << OutputFilename.c_str() << "'."<< endl;
	    exit (FAILURE);
	}
	/// write data to output file
	else
	{
	    cout << "WriteOutAnalysedData:  Writing output file '" << OutputFilename.c_str() << "' ..." << endl;

	    for (unsigned int row = 0; row < OutputFileHeader.size(); row++)
	    {
		Outputfile << setw(61) << left << OutputFileHeader[row] << endl;      /// header
		if (false && Debug) cout << setw(61) << left << OutputFileHeader[row] << endl;
	    }
	    for (unsigned int row = 0; row < WriteOutTable.size(); row++)                  /// data
	    {
		for (unsigned int col = 0; col < WriteOutTable[row].size(); col++)
		{
		    Outputfile << scientific << setw(30) << left << setprecision(8) << WriteOutTable[row][col];
		    if (false && Debug) cout << scientific << setw(30) << left << setprecision(8) << WriteOutTable[row][col];
		}
		Outputfile << endl; if (false && Debug) cout << endl;
	    }

	    Outputfile.close();
	    Outputfile.clear();

	    cout << "WriteOutAnalysedData:  done!" << endl;
	}
} /// =======================================================================


//TODO: Write this function for arbitary dataset
/**--------------------------------------ReadData -----------------------------------------------------
 ** Reads Data from PLUTO dbl file and initializes it to pointer data_ptr
 ** ---------------------------------------------------------------------------------------------------- */
// Remember N[0] : X direction, N[1]: Y direction, N[2]: Z direction
void ReadData(const string inputfile,  double * data_ptr, vector<int> N)   
{

	///open input file
	FILE *FP ;
	int i,j,k;
	int IBEG = 0, JBEG=0, KBEG=0; //Hopefully true
        int IEND= N[X]-1 ,JEND=N[Y]-1,KEND=N[Z]-1 ; // Again hopefully true
	char * Vc; 
	size_t dsize = sizeof(double); // TODO: Change this if working with single precision doubles
	FP = fopen(inputfile.c_str(),"rb");
	if (FP == NULL){
        	cout<<"! OpenBinaryFile: file"<<inputfile<<"does not exist"<<endl ; 
        exit(EXIT_FAILURE);
    }

        ///opened binary file
	
	Vc = (char *) data_ptr;    
	for (k = KBEG; k <= KEND; k++) {
         	for (j = JBEG; j <= JEND; j++) {
           		i = IBEG + N[X]*(j + N[Y]*k);
           		fread (Vc + i*dsize, dsize, N[X] , FP);
   }}

/// SwapMemOrder
	
	SwapMemOrder(data_ptr,N);

	if (Debug) {
		double sum = 0.0, sum_red = 0.0;
		for (long n = 0; n < N[X]*N[Y]*N[Z]; n++) sum += data_ptr[n];
	    MPI_Allreduce(&sum, &sum_red, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if (MyPE==0) cout<<"after SwapMemOrder: sum="<<sum_red<<endl;
	}

	
}
/// =======================================================================

/// SwapMemOrder =====================================================
void SwapMemOrder(double * data, vector<int> N)
{
	const long ntot = N[0]*N[1]*N[2];
	float *tmp = new float[ntot];
	for (int i=0; i<N[0]; i++) for (int j=0; j<N[1]; j++) for (int k=0; k<N[2]; k++) {
		long ind1 = i*N[1]*N[2] + j*N[2] + k;
		long ind2 = k*N[1]*N[0] + j*N[0] + i;
		tmp[ind1] = data[ind2];
	}
	for (long i=0; i<ntot; i++) data[i] = tmp[i];
	delete [] tmp;
} /// ==============================================================


