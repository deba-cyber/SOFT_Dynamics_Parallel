# if ! defined (_DATA_PROCESS_H)
# define _DATA_PROCESS_H

# include <iosfwd>
# include <vector>
# include <fftw3.h>


/***----------------------------------------------------------***/
/***----------------------------------------------------------***/

/**-------------------------------------------**/
/* Class for stride array and corresponding 
 * member functions .. will be used for 
 * extracting relevant elements at queried 
 * multidimensional grid points ..........
**---------------------------------------------**/

class STRIDE_ARR_BACK_N_FORTH_DP	{
	public:
		std::vector<int> basis_size_vect;
		STRIDE_ARR_BACK_N_FORTH_DP(const std::vector<int>& basis_size_vect);
		std::vector<int> stride_arr();
		int multidim_index_dpb(const std::vector<int>& onedim_index_vect);
		std::vector<int> onedim_indices(const int& multidim_index);
};

/***----------------------------------------------------------***/

/** Routine for calculating overlap of wave packet (important for 
 *  normalization constant) ..
 *  overlap can be calculated from this same routine for wave packet 
 *  in momentum space **/ 

double get_normalizn_fact(const std::vector<fftw_complex>& init_wp,
		STRIDE_ARR_BACK_N_FORTH_DP& stridemultobj,const std::vector<std::vector<double>>& all_onedimgrids,
		const std::vector<double>& stepsizevect);


/***----------------------------------------------------------***/

/** Function for calculating expectation position/momentum value for all coordinates of a multidimensional wave function 
 *  wave function is passed in argument(data type = fftw_complex) .. multidim function as 1D vector **/ 

std::vector<double> get_expect_pos_or_momentum(const std::vector<fftw_complex>& wp,STRIDE_ARR_BACK_N_FORTH_DP& stridemultdimobj,
		const std::vector<std::vector<double>>& all_onedimgrids,const std::vector<double>& stepsizevect,
		double init_state_overlap);	

/***----------------------------------------------------------***/
/** Function for calculating <x**2> for all coordinates .. to 
 * be used for dispersion of wave packet in coordinate space **/

std::vector<double> get_expect_pos_sqr(const std::vector<fftw_complex>& wp,STRIDE_ARR_BACK_N_FORTH_DP& stridemultdimobj,
		const std::vector<std::vector<double>>& all_onedimgrids,const std::vector<double>& stepsizevect,
		double init_state_overlap);


/***----------------------------------------------------------***/

/** Function for numerical integration using 1D Simpson 
 *  passed argument has odd number of points **/

double onedim_simpson(const std::vector<fftw_complex>& vect,double stepsize);


/** Function template for 1D simpson integration with usual data types **/

template<typename T>
double onedim_simpson_template(const std::vector<T>& vect,double stepsize)	{
	double integral = pow(vect[0],2.)+pow(vect[vect.size()-1],2.);
	for (unsigned int i=1; i<vect.size(); i++)	{
		if (i%2 != 0)
			integral += 4.*(pow(vect[i],2.));
		else if (i%2 == 0)
			integral += 2.*(pow(vect[i],2.));
	}
	integral*= (stepsize/3);
	return integral;
}

/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/

/** Routine to calculate overlap of wave packet at different times with 
 * prepared direct product states ...
 * prepared direct product states are real 
 * overlap calculation has both real and imaginary parts 
 * both real and imaginary parts will be returned ..
 * data type fftw_complex **/


std::vector<fftw_complex> get_overlap_w_dp(const std::vector<fftw_complex>& wp_cur,const std::vector<double>& Q1_anh_quanta_val,
		const std::vector<double>& Q10_anh_quanta_val,const std::vector<double>& Q29_anh_quanta_val,
		STRIDE_ARR_BACK_N_FORTH_DP& stridemultdimobj,const std::vector<double>& stepsizevect);


/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/

/**** Routine for calculation of overlap of wave function with 
 * itself when the first coordinate range is halved ***/

double get_tun_prob(const std::vector<fftw_complex>& wp,STRIDE_ARR_BACK_N_FORTH_DP& stridefullgridobj,
		STRIDE_ARR_BACK_N_FORTH_DP& stridehalfgridobj,int DPG_IO_size,int shiftval,
		const std::vector<double>& stepsizevect_IO);




/*-------------------------------------------------------------*/


/** Routine to calculate reduced probability density by integrating
 *  out multiple coordinates  ...
 *  Full Dim (X1 X2 X3 X4 ... XN)
 *  Arguments mean ...
 *  Full wave function is given as 1D vector(std::vector<fftw_complex>)
 *  std::vector<int> query_multgrid_id has size equal to full dimension
 *  of the dynamics ... 
 *  coordinates not integrated out are set at const 
 *  values where reduced prob. density is desired ..  
 *  coordinates to be integrated out are set to zero (as those will change)
 *  std::vector<int> serialvect_IO .. have indices in serial order for 
 *  the coordinates to be integrated out (0-based indexing) 
 *  stride_IO_obj corresponds to stride object involving the coordinates 
 *  to be integrated out ...
 *  stride_fulldim_obj corresponds to stride object involving all coordinates
 *  in the dynamics ....
 *  DPG_IO_size corresponds to direct product grid size for the coordinates 
 *  to be integrated out ...
 *  std::vector<double> stepsizevect_IO .. grid spacing for the coordinates 
 *  to be integrated out (needed to calculate the weight factor)
 **/

double get_reduced_prob_density_IO_multcoord(const std::vector<fftw_complex>& wp,const std::vector<int>& query_multgrid_id,
		const std::vector<int>& serialvect_IO,STRIDE_ARR_BACK_N_FORTH_DP& stride_IO_obj,
		STRIDE_ARR_BACK_N_FORTH_DP& stride_fulldim_obj,int DPG_IO_size,const std::vector<double>& stepsizevect_IO);



/***----------------------------------------------------------***/
/***----------------------------------------------------------***/

/** Routine to calculate probability density by integrating out one coordinate ..
 *	Full Dim (X1 X2 X3 X4 ... XN)
 *	Arguments mean ...
 *	Full wave function is given as 1D vector (std::vector<fftw_complex>)
 *	std::vector<int> query_multgrid_id has size equal to full dimension of the dynamics
 *	Coordinates not integrated out are set constant values where reduced prob. density 
 *	is desired ...
 *	Coordinates to be integrated out are set to zero (as those will change)  
 *  int id_coord_IO is the index of the coordinate to be integrated out.(0-based indexing) 
 *  stride_fulldim_obj corresponds to stride object involving all coordinates in the dynamics ....
 *  IO_coord_size .. no of grid points for the relevant coordinate to be integrated out (for the 1D integration)
 *  double stepsize ... grid spacing for the coordinate to be integrated out
 *  **/


double get_reduced_prob_density_IO_coord(const std::vector<fftw_complex>& wp,const std::vector<int>& query_multgrid_id,
		int id_coord_IO,STRIDE_ARR_BACK_N_FORTH_DP& stride_fulldim_obj,const int IO_coord_size,const double stepsize);


/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/




















# endif
