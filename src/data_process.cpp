# include <SOFT/data_process.hpp>
# include <vector>
# include <cmath>
# include <algorithm>


// constructor for STRIDE_ARR_BACK_N_FORTH_DP class for data process //


STRIDE_ARR_BACK_N_FORTH_DP::STRIDE_ARR_BACK_N_FORTH_DP(const std::vector<int>& basis_size_vect)	{
	this->basis_size_vect = basis_size_vect;
}

std::vector<int> STRIDE_ARR_BACK_N_FORTH_DP::stride_arr()	{
	/** preparing stride array **/
	int cur_product;
	int cur_index_pdt = 1;
	for (size_t j=1; j<basis_size_vect.size(); j++)	{
		cur_index_pdt*= basis_size_vect[j];
	}
	/** Initialising stride array with first element **/
	std::vector<int>stride_arr_init;
	stride_arr_init.emplace_back(cur_index_pdt);
	/** other elements of stride array will be prepared from the first element of the stride array **/
	int cur_index_init = 1;		// initialising current index for generating other elements of stride array 
	while (true)	{
		if (cur_index_init == basis_size_vect.size()-1)	{
			break;
	}
		else 
		{
			cur_product = int(cur_index_pdt/basis_size_vect[cur_index_init]);
			cur_index_init +=1;
			stride_arr_init.emplace_back(cur_product);
			cur_index_pdt = cur_product;
		}
	}
	return stride_arr_init;
}


/**********************************************************
 ********* member function to generate multidimensional 
			index from onedimensional index array  *******
 **********************************************************/

int STRIDE_ARR_BACK_N_FORTH_DP::multidim_index_dpb(const std::vector<int> &onedim_index_vect)	{
	/** Given one dimensional indices , will return multidimensional index 
	 ** array/vector containing one dimensional indices are zero-based indices **/
	std::vector<int>stride_arr = this->stride_arr();	// calling stride array
	int multidim_basis_index = onedim_index_vect.back();
	multidim_basis_index++;		// adding one for 1-based indexing 
	for (size_t i=0; i<stride_arr.size(); i++)	{
		multidim_basis_index += stride_arr[i]*onedim_index_vect[i];
	}
	return multidim_basis_index;
}


/********************************************************************************************
 * ********* member function to one dimensional index array (returned as vector)	*********
 *	*******		from multidimensional index (1-based indexing)
 *******************************************************************************************/
std::vector<int> STRIDE_ARR_BACK_N_FORTH_DP::onedim_indices(const int & multidim_index)	{
	/** Given multidimensional index for direct product basis, returns
	 ** one dimensional indices ( 0- based indexing ); multidim_index has 1-based indexing ***/
	std::vector<int>stride_arr = this->stride_arr();	// calling stride array 
	int multidim_index_4_caln;
	multidim_index_4_caln = multidim_index -1;
	std::vector<int>onedim_index_vect;
    /** multidim_index will change for finding each of the 1-d indices in the loop .. 
	 * here it is first initialized **/
	int cur_onedim_index;
	for (size_t i=0; i<stride_arr.size(); i++)	{
		cur_onedim_index = int(multidim_index_4_caln/stride_arr[i]);
		onedim_index_vect.emplace_back(cur_onedim_index);
		multidim_index_4_caln -= cur_onedim_index*stride_arr[i];
	}
	onedim_index_vect.emplace_back(multidim_index_4_caln);
	/** returns 1 dimensional index array .. zero based indexing **/
	return onedim_index_vect;
}


/***----------------------------------------------------------***/

/** Routine for calculating overlap of wave packet (important for 
 *  normalization constant) ..
 *  overlap can be calculated from this same routine for wave packet 
 *  in momentum space **/ 

double get_normalizn_fact(const std::vector<fftw_complex>& init_wp,
		STRIDE_ARR_BACK_N_FORTH_DP& stridemultobj,const std::vector<std::vector<double>>& all_onedimgrids,
		const std::vector<double>& stepsizevect)	{
	double integral_val = 0.;
	for (unsigned int i=0; i<init_wp.size(); i++)	{
		std::vector<int> cur_id_vect = stridemultobj.onedim_indices(i+1);
		double cur_wt = 1.;
		for (unsigned int j=0; j<cur_id_vect.size(); j++)	{
			if (cur_id_vect[j] == 0 || cur_id_vect[j] == stridemultobj.basis_size_vect[j]-1)
				cur_wt *= 1.;
			else if (cur_id_vect[j] %2 == 0)	
				cur_wt *= 2.;
			else if (cur_id_vect[j] %2 != 0)	
				cur_wt *= 4.;
		}
		integral_val += cur_wt*(pow(init_wp[i][0],2.) + pow(init_wp[i][1],2.));
	}
	for (unsigned int i=0; i<stepsizevect.size(); i++)
		integral_val *= (stepsizevect[i]/3);
	return integral_val;
}

/***----------------------------------------------------------***/

/** Function for calculating <x**2> for all coordinates .. to 
 * be used for dispersion of wave packet in coordinate space **/

std::vector<double> get_expect_pos_sqr(const std::vector<fftw_complex>& wp,STRIDE_ARR_BACK_N_FORTH_DP& stridemultdimobj,
		const std::vector<std::vector<double>>& all_onedimgrids,const std::vector<double>& stepsizevect,
		double init_state_overlap)	{
	std::vector<double> all_expect_vect(all_onedimgrids.size());
	for (unsigned int i=0; i<wp.size(); i++)	{
		std::vector<int> cur_id_vect = stridemultdimobj.onedim_indices(i+1);
		double cur_wt = 1.;
		for (unsigned int j=0; j<cur_id_vect.size(); j++)	{
			if (cur_id_vect[j] == 0 || cur_id_vect[j] == stridemultdimobj.basis_size_vect[j]-1)	
				cur_wt *= 1.;
			else if (cur_id_vect[j] %2 == 0)	
				cur_wt *= 2.;
			else if (cur_id_vect[j] %2 != 0)	
				cur_wt *= 4.;
		}
		all_expect_vect[0] += cur_wt*(pow(wp[i][0],2.)+pow(wp[i][1],2.))*(pow(all_onedimgrids[0][cur_id_vect[0]],2.));
		all_expect_vect[1] += cur_wt*(pow(wp[i][0],2.)+pow(wp[i][1],2.))*(pow(all_onedimgrids[1][cur_id_vect[1]],2.));
		all_expect_vect[2] += cur_wt*(pow(wp[i][0],2.)+pow(wp[i][1],2.))*(pow(all_onedimgrids[2][cur_id_vect[2]],2.));
	}
	for (unsigned int i=0; i<stepsizevect.size(); i++)	{
		all_expect_vect[0] *= (stepsizevect[i]/3);
		all_expect_vect[1] *= (stepsizevect[i]/3);
		all_expect_vect[2] *= (stepsizevect[i]/3);
	}
	for_each(all_expect_vect.begin(),all_expect_vect.end(),[&](double& elem) { elem *= (1/init_state_overlap);});
	return all_expect_vect;
}

/***----------------------------------------------------------***/

/** Function for calculating expectation position/momentum value for all coordinates of a multidimensional wave function 
 *  wave function is passed in argument(data type = fftw_complex) .. multidim function as 1D vector **/ 

std::vector<double> get_expect_pos_or_momentum(const std::vector<fftw_complex>& wp,STRIDE_ARR_BACK_N_FORTH_DP& stridemultdimobj,
		const std::vector<std::vector<double>>& all_onedimgrids,const std::vector<double>& stepsizevect,
		double init_state_overlap)	{
	std::vector<double> all_expect_vect(all_onedimgrids.size());
	for (unsigned int i=0; i<wp.size(); i++)	{
		std::vector<int> cur_id_vect = stridemultdimobj.onedim_indices(i+1);	
		double cur_wt = 1.;
		for (unsigned int j=0; j<cur_id_vect.size(); j++)	{
			if (cur_id_vect[j] == 0 || cur_id_vect[j] == stridemultdimobj.basis_size_vect[j]-1)	
				cur_wt *= 1.;
			else if (cur_id_vect[j] %2 == 0)	
				cur_wt *= 2.;
			else if (cur_id_vect[j] %2 != 0)	
				cur_wt *= 4.;
		}
		all_expect_vect[0] += cur_wt*(pow(wp[i][0],2.)+pow(wp[i][1],2.))*all_onedimgrids[0][cur_id_vect[0]]; 
		all_expect_vect[1] += cur_wt*(pow(wp[i][0],2.)+pow(wp[i][1],2.))*all_onedimgrids[1][cur_id_vect[1]]; 
		all_expect_vect[2] += cur_wt*(pow(wp[i][0],2.)+pow(wp[i][1],2.))*all_onedimgrids[2][cur_id_vect[2]];
	}	
	for (unsigned int i=0; i<stepsizevect.size(); i++)	{
		all_expect_vect[0] *= (stepsizevect[i]/3);
		all_expect_vect[1] *= (stepsizevect[i]/3);
		all_expect_vect[2] *= (stepsizevect[i]/3);
	}
	for_each(all_expect_vect.begin(),all_expect_vect.end(),[&](double& element) {element *= (1/init_state_overlap);});	
	return all_expect_vect;
}

/***----------------------------------------------------------***/
/***----------------------------------------------------------***/

/** Function for numerical integration using 1D Simpson 
 *  passed argument has odd number of points **/

double onedim_simpson(const std::vector<fftw_complex>& vect,double stepsize)	{
	double integral = pow(vect[0][0],2.)+pow(vect[0][1],2.)+pow(vect[vect.size()-1][0],2.)+pow(vect[vect.size()-1][1],2.);
	for (unsigned int i=1; i<vect.size()-1; i++)	{
		if (i%2 != 0)
			integral += 4.*(pow(vect[i][0],2.)+pow(vect[i][1],2.));		
		else if (i%2 == 0)
			integral += 2.*(pow(vect[i][0],2.)+pow(vect[i][1],2.));
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
		STRIDE_ARR_BACK_N_FORTH_DP& stridemultdimobj,const std::vector<double>& stepsizevect)	{
	std::vector<fftw_complex> overlap_fftw(1);
	double overlap_re = 0.;
	double overlap_im = 0.;
	for (unsigned int i=0; i<wp_cur.size(); i++)	{
		std::vector<int> cur_id_vect = stridemultdimobj.onedim_indices(i+1);
		double DP_Q1_Q10_Q29 = Q1_anh_quanta_val[cur_id_vect[0]]*Q10_anh_quanta_val[cur_id_vect[1]]*
			Q29_anh_quanta_val[cur_id_vect[2]];
		double cur_wt = 1.;
		for (unsigned int j=0; j<cur_id_vect.size(); j++)	{
			if (cur_id_vect[j] == 0 || cur_id_vect[j] == stridemultdimobj.basis_size_vect[j]-1)	
				cur_wt *= 1.;
			else if (cur_id_vect[j] %2 == 0)	
				cur_wt *= 2.;
			else if (cur_id_vect[j] %2 != 0)	
				cur_wt *= 4.;
		}
		overlap_re += cur_wt*DP_Q1_Q10_Q29*wp_cur[i][0];
		overlap_im += cur_wt*DP_Q1_Q10_Q29*wp_cur[i][1];
	}
	for (unsigned int i=0; i<stepsizevect.size(); i++)	{
		overlap_re *= (stepsizevect[i]/3);
		overlap_im *= (stepsizevect[i]/3);
	}
	overlap_fftw[0][0] = overlap_re;
	overlap_fftw[0][1] = overlap_im;
	return overlap_fftw;
}


/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/

/**** Routine for calculation of overlap of wave function with 
 * itself when the first coordinate range is halved ***/

double get_tun_prob(const std::vector<fftw_complex>& wp,STRIDE_ARR_BACK_N_FORTH_DP& stridefullgridobj,
		STRIDE_ARR_BACK_N_FORTH_DP& stridehalfgridobj,int DPG_IO_size,int shiftval,
		const std::vector<double>& stepsizevect_IO)	{
	double integral = 0.;
	std::vector<int> cur_id_vect_full(stepsizevect_IO.size());
	for (unsigned int i=0; i<DPG_IO_size; i++)	{
		std::vector<int> cur_id_vect_IO = stridehalfgridobj.onedim_indices(i+1);
		int j=0.;
		for_each(cur_id_vect_full.begin(),cur_id_vect_full.end(),[&cur_id_vect_IO,&j](int& elem){
				elem = cur_id_vect_IO[j];
				j++;});
		cur_id_vect_full[0] += shiftval;
		int cur_fulldim_index = stridefullgridobj.multidim_index_dpb(cur_id_vect_full);
		double cur_wt = 1.;
		for (unsigned int k=0; k<cur_id_vect_IO.size(); k++)	{
			if (cur_id_vect_IO[k] == 0 || cur_id_vect_IO[k] == stridehalfgridobj.basis_size_vect[j]-1)
				cur_wt *= 1.;
			else if (cur_id_vect_IO[k] %2 == 0)
				cur_wt *= 2.;
			else if (cur_id_vect_IO[k] %2 != 0)
				cur_wt *= 4.;
		}
		integral += cur_wt*(pow(wp[cur_fulldim_index-1][0],2.)+pow(wp[cur_fulldim_index-1][1],2.)); 
	}
	for (unsigned int i=0; i<stepsizevect_IO.size(); i++)	
		integral *= (stepsizevect_IO[i]/3);
	return integral;
}

/***----------------------------------------------------------***/
/***----------------------------------------------------------***/


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
		STRIDE_ARR_BACK_N_FORTH_DP& stride_fulldim_obj,int DPG_IO_size,const std::vector<double>& stepsizevect_IO)	{
	double red_prob_density = 0.;
	std::vector<int> cur_id_vect_fulldim = query_multgrid_id;
	for (unsigned int i=1; i<DPG_IO_size+1; i++)	{
		std::vector<int> cur_id_vect_IO = stride_IO_obj.onedim_indices(i);
		int j=0;
		for_each(cur_id_vect_IO.begin(),cur_id_vect_IO.end(),[&cur_id_vect_fulldim,&serialvect_IO,&j](int& elem){cur_id_vect_fulldim[serialvect_IO[j]] = elem;
				j++;});
		int cur_fulldim_index = stride_fulldim_obj.multidim_index_dpb(cur_id_vect_fulldim); // 1-based indexing 
		double cur_wt = 1.;
		for (unsigned int k=0; k<cur_id_vect_IO.size(); k++)	{
			if (cur_id_vect_IO[k] == 0 || cur_id_vect_IO[k] == stride_IO_obj.basis_size_vect[j]-1)
				cur_wt *= 1.;
			else if (cur_id_vect_IO[k] %2 == 0)
				cur_wt *= 2.;
			else if (cur_id_vect_IO[k] %2 != 0)
				cur_wt *= 4.;
		}
		red_prob_density += cur_wt*(pow(wp[cur_fulldim_index-1][0],2.)+pow(wp[cur_fulldim_index-1][1],2.));
	}
	for (unsigned int i=0; i<stepsizevect_IO.size(); i++)	
		red_prob_density *= (stepsizevect_IO[i]/3);
	return red_prob_density;
}



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
		int id_coord_IO,STRIDE_ARR_BACK_N_FORTH_DP& stride_fulldim_obj,const int IO_coord_size,const double stepsize)	{
	std::vector<int> cur_id_vect_fulldim = query_multgrid_id;
	std::vector<fftw_complex> wp_desired_loc(IO_coord_size);
	for (unsigned int i=0; i<IO_coord_size; i++)	{
		cur_id_vect_fulldim[id_coord_IO] = i;
		int cur_fulldim_index = stride_fulldim_obj.multidim_index_dpb(cur_id_vect_fulldim);
		wp_desired_loc[i][0] = wp[cur_fulldim_index-1][0];
		wp_desired_loc[i][1] = wp[cur_fulldim_index-1][1];
	}
	double red_prob_density = onedim_simpson(wp_desired_loc,stepsize);
	return red_prob_density;
}


/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/









