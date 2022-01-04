

# include <SOFT/tunneling_init_state.hpp>
# include <SOFT/fileops.hpp>
# include <SOFT/constants.hpp>
# include <iostream>
# include <cmath>
# include <string>
# include <vector>
# include <numeric>
# include <algorithm>
# include <boost/math/special_functions/hermite.hpp>


/*******************************************************
 * class constructor for class STRIDE_ARR_BACK_N_FORTH *
 * *****************************************************/

STRIDE_ARR_BACK_N_FORTH::STRIDE_ARR_BACK_N_FORTH(int* basis_size_arr,int size)	{
	this->basis_size_arr = new int[size];
	this->size = size;
	for (int i=0; i<size; i++)	{
		this->basis_size_arr[i] = basis_size_arr[i];
	}
}

/**********************************************************
 ********* member function to generate stride array *******
 **********************************************************/

std::vector<int> STRIDE_ARR_BACK_N_FORTH::stride_arr()	{
	/** preparing stride array **/
	int cur_product;
	int cur_index_pdt = 1;
	for (int j=1; j<size; j++)	{
		cur_index_pdt*= basis_size_arr[j];
	}
	/** Initialising stride array with first element **/
	std::vector<int>stride_arr_init;
	stride_arr_init.emplace_back(cur_index_pdt);
	/** other elements of stride array will be prepared from the first element of the stride array **/
	int cur_index_init = 1;		// initialising current index for generating other elements of stride array 
	while (true)	{
		if (cur_index_init == size-1)	{
			break;
	}
		else 
		{
			cur_product = int(cur_index_pdt/basis_size_arr[cur_index_init]);
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

int STRIDE_ARR_BACK_N_FORTH::multidim_index_dpb(const std::vector<int> &onedim_index_vect)	{
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
std::vector<int> STRIDE_ARR_BACK_N_FORTH::onedim_indices(const int & multidim_index)	{
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

/*******************************************************
 * class destructor for class STRIDE_ARR_BACK_N_FORTH *
 * *****************************************************/

STRIDE_ARR_BACK_N_FORTH::~STRIDE_ARR_BACK_N_FORTH()	{
	delete[] basis_size_arr;
}

/**----------------------------------------------------------------**/
/**----------------------------------------------------------------**/
/**----------------------------------------------------------------**/
/**----------------------------------------------------------------**/

/***************************************************
 * class constructor for derived class init_wp_val *
 * *************************************************/

init_wp_val::init_wp_val(int* basis_size_arr, int size, int DPB_size,
int* mode_id_arr, const std::vector<std::string> &dvr_type_vect, double* Q_shift_arr, int tot_prim_dvr_pts,
int* n_row_so_2_anh_arr, int* n_col_so_2_anh_arr, int* n_row_anh_2_po_arr, int* n_col_anh_2_po_arr,
int tot_size_so_2_anh_all, int tot_size_anh_2_po_all,
double* so_2_anh_all_arr, double* anh_2_po_all_arr, double* sodvr_pts_arr, 
const std::string & dyn_type)	:STRIDE_ARR_BACK_N_FORTH(basis_size_arr, size)	{
	this-> DPB_size = DPB_size;
	this-> tot_prim_dvr_pts = tot_prim_dvr_pts;
	this-> mode_id_arr = new int[size];
	this-> Q_shift_arr = new double[size];
	this-> n_row_so_2_anh_arr = new int[size];
	this-> n_col_so_2_anh_arr = new int[size];
	this-> n_row_anh_2_po_arr = new int[size];
	this-> n_col_anh_2_po_arr = new int[size];
	this-> sodvr_pts_arr = new double[this-> tot_prim_dvr_pts];
	this-> tot_size_so_2_anh_all = tot_size_so_2_anh_all;
	this-> tot_size_anh_2_po_all = tot_size_anh_2_po_all;
	for (int i=0; i<size; i++)	{
		this-> mode_id_arr[i] = mode_id_arr[i];
		this-> Q_shift_arr[i] = Q_shift_arr[i];
		this-> n_row_so_2_anh_arr[i] = n_row_so_2_anh_arr[i];
		this-> n_col_so_2_anh_arr[i] = n_col_so_2_anh_arr[i];
		this-> n_row_anh_2_po_arr[i] = n_row_anh_2_po_arr[i];
		this-> n_col_anh_2_po_arr[i] = n_col_anh_2_po_arr[i];
	}
	this-> so_2_anh_all_arr = new double [tot_size_so_2_anh_all];
	this-> anh_2_po_all_arr = new double [tot_size_anh_2_po_all];
	for (int i=0; i<(this->tot_prim_dvr_pts); i++)	{
		this-> sodvr_pts_arr[i] = sodvr_pts_arr[i];
	}
	for (int i=0; i<tot_size_so_2_anh_all; i++)	{
		this-> so_2_anh_all_arr[i] = so_2_anh_all_arr[i];
	}
	for (int i=0; i<tot_size_anh_2_po_all; i++)	{
		this-> anh_2_po_all_arr[i] = anh_2_po_all_arr[i];
	}
	this-> dvr_type_vect = dvr_type_vect;
	this-> dyn_type = dyn_type;
}

/**********************************************************************
 * Member function to get one dimensional index arrays for all values 
 * up to DPB size	..	by calling function onedim_indices 
 * of STRIDE_ARR_BACK_N_FORTH 
**********************************************************************/

std::vector<int> init_wp_val::get_1_dim_index_vect()	{
  /** all one dimensional indices for all DPB saved in one dimensional vector 
   *  Indexing for DPB starts from 1 ........................................
   *  Returning as 1-dimensional vector	***/
	std::vector<int> onedim_id_full_vect;	
	std::vector<int> cur_onedim_index_vect;
	for (int i=1; i< (this->DPB_size +1); i++)	{
		cur_onedim_index_vect = STRIDE_ARR_BACK_N_FORTH::onedim_indices(i);
		onedim_id_full_vect.insert(std::end(onedim_id_full_vect),std::begin(cur_onedim_index_vect),std::end(cur_onedim_index_vect));
	}
	return onedim_id_full_vect;
}

/***********************************************************************
 * Member function to get value of the HO function  at the queried
 * grid point	........................................................
 * *********************************************************************/


std::vector<double> init_wp_val::prim_basis_val_vect(int N_HO_basis, double Q_cur_val)	{
/** N_HO_basis is no of HO functions for which values 
   * will be calculated at the current grid value **/
	double HO_cur_func_at_grid;
	std::vector<double> HO_func_all_at_grid;
	for (int i=0; i<N_HO_basis; i++)	{
		HO_cur_func_at_grid = (1/sqrt((pow(2,i))*(std::tgamma(i+1))))*
			(1/(pow(M_PI,0.25)))*(boost::math::hermite(i,Q_cur_val))*(exp(-0.5*(pow(Q_cur_val,2))));
		HO_func_all_at_grid.emplace_back(HO_cur_func_at_grid);
	}
/** returning values of N_HO_basis at the queried grid point **/
	return HO_func_all_at_grid;
}



/**************************************************************************
 * Member function for returning values of all primitive DVR functions  
 * (N_prim for the concerned coordinate at the current coordinate value
 * If sinc-dvr is used for a particular coordinate, shift value is kept 
 * at zero	............................................................
 * ************************************************************************/



std::vector<double> init_wp_val::dvr_basis_val_vect(int mode_id, double cur_coord_val)	{
  // getting the location of the queried mode , 0-based indexing // 
	int mode_index_id = std::distance(this-> mode_id_arr, std::find(this-> mode_id_arr, this-> mode_id_arr+size, mode_id));
  // Getting DVR type (HO/SINC) for the current mode //
	std::string cur_dvr_type = this->dvr_type_vect[mode_index_id];
  /* Getting 1-dim dvr points for the relevant coordinate
   * no of points for an individual coordinate is equal to 
   * elements in column n_col_so_2_anh_arr ***********/ 
	std::vector<double> dvr_pts_cur; 
	int start_id = 0;
	int end_id = 0;
	for (int i=0; i<mode_index_id; i++)	{
		start_id += this-> n_col_so_2_anh_arr[i];
		end_id += this-> n_col_so_2_anh_arr[i];
	}
	end_id += this-> n_col_so_2_anh_arr[mode_index_id];
	int N_prim = end_id - start_id;		// No of primitive basis used
	//////////////////////////////////
	for (int i=start_id; i<end_id; i++)	{
		dvr_pts_cur.emplace_back(this->sodvr_pts_arr[i]);
	}
	///////////////////////////////////
	double delta = (dvr_pts_cur.back()-dvr_pts_cur.front())/(dvr_pts_cur.size()-1);
	double arg1;	// argument of sin theta in sinc-dvr 
	if (cur_dvr_type == "sinc")	{
		std::vector<double>sinc_dvr_all_at_grid;
		for (int i=0; i<N_prim; i++)	{
			if (cur_coord_val != dvr_pts_cur[i])	{
				arg1 = (M_PI/delta)*(cur_coord_val - dvr_pts_cur[i]); 
				sinc_dvr_all_at_grid.emplace_back((1/sqrt(delta))*(sin(arg1))/(arg1));
			}
			else 
				sinc_dvr_all_at_grid.emplace_back(1.);
		}
		return sinc_dvr_all_at_grid;
  }
  else if (cur_dvr_type == "HO")	{
	  // constructing FBR ==> DVR transformation matrix .. again as 1-d vector//
	  std::vector<double> FBR_2_DVR_vect;	// total no of elements N_prim*N_prim
	  double fbr_2_dvr_cur_val;		// placeholder for current fbr_2_dvr matrix element
	  /** Hermite polynomial nodes and weights corresponding to N_prim size **/
	  std::vector<double> herm_nodes = get_specific_line(N_prim,"../input_data/Gauss_Hermite_nodes") ;
	  std::vector<double> herm_wts = get_specific_line(N_prim,"../input_data/Gauss_Hermite_weights");
	  for (int i=0; i<N_prim; i++)	{
		  for (int j=0; j<N_prim; j++)	{
			  fbr_2_dvr_cur_val = (sqrt(herm_wts[i]))*((1/sqrt((pow(2,j))*(std::tgamma(j+1))))
					  *(1/pow(M_PI,0.25)))*(boost::math::hermite(j,herm_nodes[i])); 
			  FBR_2_DVR_vect.emplace_back(fbr_2_dvr_cur_val);
		  }
	  }
	  ///////////////////////////////////////////////////
	  double Q_shift_cur = this-> Q_shift_arr[mode_index_id];	// extracting the shift for DVR calculation for the current coordinate
	  int N_HO = this-> n_col_so_2_anh_arr[mode_index_id];		// extracting the number of HO basis functions //
	  /** value of the wavefunction to be evaluated when there's no shift **/
	  double cur_coord_wo_shift = cur_coord_val - Q_shift_cur;
	  /** values of all HO functions as a vector when there's no shift **/
	  std::vector<double> HO_func_all_vect = this-> prim_basis_val_vect(N_HO,cur_coord_wo_shift);
	  // storing values of all HO-DVR functions at the current grid //
	  std::vector<double> HO_dvr_all_at_cur_grid = get_mat_vec_pdt(FBR_2_DVR_vect,HO_func_all_vect,N_HO,N_HO);
	  return HO_dvr_all_at_cur_grid;
  }
}


/**********************************************************************************
 * Member function for returning value of an one dimensional podvr basis function
 * at the corresponding coordinate value of the multidimensional grid	...
 * cur_coord_po_or_anh_basis_id (zero-based indexing) gives 1-d po/anh function 
 * index value .. based on dynamics type (IVR/TUNNELING) anh/po will be used
 * of which is to be evaluated at the current grid	.............................
 * It'll need transformation matrices SODVR=>ANH_EIG=>PODVR .....................
 * ********************************************************************************/


double init_wp_val::pobasis_val(int mode_id, double cur_coord_val, int cur_coord_po_or_anh_basis_id)	{
  // getting the location of the queried mode , 0-based indexing // 
	int mode_index_id = std::distance(this-> mode_id_arr, std::find(this-> mode_id_arr, this-> mode_id_arr+size, mode_id));
  // relevant transformation matrix for DVR => ANH as 1-dimensional vector //
	std::vector<double> dvr_2_anh_cur_vect;
  // relevant transformation matrix for ANH ==> PO as 1-dimensional vector //
	std::vector<double> anh_2_po_cur_vect;
	//////////////////////////////////////
	int start_id1 = 0;
	int end_id1 = 0;
	int start_id2 = 0;
	int end_id2 = 0;
	for (int i=0; i<mode_index_id; i++)	{
		start_id1 += (this-> n_row_so_2_anh_arr[i])*(this-> n_col_so_2_anh_arr[i]); 
		end_id1 += (this-> n_row_so_2_anh_arr[i])*(this-> n_col_so_2_anh_arr[i]);
		start_id2 += (this-> n_row_anh_2_po_arr[i])*(this-> n_col_anh_2_po_arr[i]);
		end_id2 += (this-> n_row_anh_2_po_arr[i])*(this-> n_col_anh_2_po_arr[i]);
	}
	end_id1 += (this-> n_row_so_2_anh_arr[mode_index_id])*(this-> n_col_so_2_anh_arr[mode_index_id]);
	end_id2 += (this-> n_row_anh_2_po_arr[mode_index_id])*(this-> n_col_anh_2_po_arr[mode_index_id]);
	//////////////////////////////////////////////////////////////
	for (int i=start_id1; i<end_id1; i++)	{
		dvr_2_anh_cur_vect.emplace_back(this-> so_2_anh_all_arr[i]); 
	}
	for (int i=start_id2; i<end_id2; i++)	{
		anh_2_po_cur_vect.emplace_back(this-> anh_2_po_all_arr[i]);
	}
	//////////////////////////////////////////////////////////////
	if (this-> dyn_type == "IVR")	{
	/** extracting the particular row (anharm state) from the dvr_2_anh_cur_vect **/
		std::vector<double> dvr_2_anh_cur_anh_vect;
	/** extracting the particular row (anharm state) from the dvr_2_anh_cur_vect	**/
		dvr_2_anh_cur_anh_vect = get_specific_row_or_col_elements(dvr_2_anh_cur_vect,
				cur_coord_po_or_anh_basis_id,this->n_row_so_2_anh_arr[mode_index_id],this->n_col_so_2_anh_arr[mode_index_id],"ROW");
		////////////////////////////////////////////////
	/** getting all primitive-dvr basis values at the queried grid **/
		std::vector<double> prim_dvr_func_all_at_grid = this-> dvr_basis_val_vect(mode_id,cur_coord_val);  
		///////////////////////////////////////////////
		/** getting the corresponding anharm state value **/
		double anh_cur_val_at_cur_coord = 0.;
		for (size_t i=0; i<dvr_2_anh_cur_anh_vect.size(); i++)	{
			anh_cur_val_at_cur_coord += dvr_2_anh_cur_anh_vect[i]*prim_dvr_func_all_at_grid[i];
		}
		return anh_cur_val_at_cur_coord;
	}
	else if (this-> dyn_type == "TUNNELING")	{
  // extracting the relevant array for the pobasis index .. resultant again 1-dim vector //
		std::vector<double> anh_2_po_cur_po_vect;
	/** extracting the particular row (pobasis) from the anh_2_po_cur_vect **/
		anh_2_po_cur_po_vect = get_specific_row_or_col_elements(anh_2_po_cur_vect,
				cur_coord_po_or_anh_basis_id,this->n_row_anh_2_po_arr[mode_index_id],this->n_col_anh_2_po_arr[mode_index_id],"ROW");
	/////////////////////////////////////////////////////////////
	/** getting all primitive-dvr basis values at the queried grid **/
		std::vector<double> prim_dvr_func_all_at_grid = this-> dvr_basis_val_vect(mode_id,cur_coord_val);  
	/** getting all converged anharm_eigvector values for the current coordinate **/
		std::vector<double> anharm_all_at_cur_grid = 
			get_mat_vec_pdt(dvr_2_anh_cur_vect,prim_dvr_func_all_at_grid,
					this->n_row_so_2_anh_arr[mode_index_id],this->n_col_so_2_anh_arr[mode_index_id]);
	////////////////////////////////////////////////////////////
	/** getting the corresponding podvr basis value **/
		double pobasis_cur_val_at_cur_coord = 0.;
		for (size_t i=0; i<anh_2_po_cur_po_vect.size(); i++)	{
			pobasis_cur_val_at_cur_coord += anh_2_po_cur_po_vect[i]*anharm_all_at_cur_grid[i];
		}
		return pobasis_cur_val_at_cur_coord;
	}
}


/**************************************************************************
 * Member function to generated meshgrid in two dimensions for fixed ******
 * values of other coordinates (zero/nonzero) ...
 * plt_coord_grid_1 & plt_coord_grid_2 are the one dimensional grids 
 * to be used for grid generation .. 
 * other_coord_val_vect gives values of other coordinates (zero/nonzero)
 * for those constant values 2-dimensional meshgrid in plt_coord_grid_1
 * & plt_coord_grid_2 will be generated	................................
 * keeping the values zero at locations of the plotting coordinates ..
 * plt_coord_1_id & plt_coord_2_id are mode numbers (say 1,5, etc) 
 * *************************************************************************/


std::vector<double> init_wp_val::generate_multidim_meshgrid(const std::vector<double> & plt_coord_grid_1, 
		const std::vector<double> & plt_coord_grid_2,const std::vector<double> & other_coord_val_vect, 
		int plt_coord_1_id, int plt_coord_2_id)	{
  /**************************************************
   * plt_coord_grid_1 has the slow index .. plt_coord_grid_2 has the fast index 
   * plt_coord_1_id & plt_coord_2_id are integers giving locations (0-based) 
   * of plotting coordinates in multidimensional calculation
   * ***********************************************/
	int plt_coord_1_index_id = std::distance(this-> mode_id_arr, std::find(this-> mode_id_arr, this-> mode_id_arr+size , plt_coord_1_id)); 
	int plt_coord_2_index_id = std::distance(this-> mode_id_arr, std::find(this-> mode_id_arr, this-> mode_id_arr+size , plt_coord_2_id)); 
	std::string str;	// checking if any other coordinate values are nonzero 
	for (size_t i=0; i<other_coord_val_vect.size(); i++)	{
		if (other_coord_val_vect[i] != 0.)	{
			str = "false";
			break;
		}
	}
	/*** vector for returning multidimensional grid as one dimensional vector ***/
	std::vector<double>multidim_meshgrid_vect;
	/////////////////////////////
	/** if other coordinates are all zero, multidimensional meshgrid is just 2-dimensional 
	 * direct product grid **/
	if (other_coord_val_vect.size()	!= INIT_STATE_CONST::NDIM_4_WP_INIT)	{
		throw "Size of other coordinate vector should be equal to dynamics dimension";
	}
	else {
		if (str == "false")	{
			std::vector<std::vector<double>> simple_2dim_grid;
			std::vector<double> tmp_2dim_grid;
			for (size_t i=0; i<plt_coord_grid_1.size(); i++)	{
				for (size_t j=0; j<plt_coord_grid_2.size(); j++)	{
					tmp_2dim_grid.emplace_back(plt_coord_grid_1[i]);
					tmp_2dim_grid.emplace_back(plt_coord_grid_2[j]);
					simple_2dim_grid.push_back(tmp_2dim_grid);
					tmp_2dim_grid.clear();
				}
			}
			/*** temporary placeholder for current multigrid vector to be inserted to the final vector ***/
			std::vector<double> tmp_cur_meshvect;
			int grid_size = (plt_coord_grid_1.size())*(plt_coord_grid_2.size());
			for (int i=0; i<grid_size; i++)	{
				tmp_cur_meshvect.insert(std::end(tmp_cur_meshvect),std::begin(other_coord_val_vect),std::end(other_coord_val_vect));
				tmp_cur_meshvect[plt_coord_1_index_id] = simple_2dim_grid[i][0];
				tmp_cur_meshvect[plt_coord_2_index_id] = simple_2dim_grid[i][1];
				multidim_meshgrid_vect.insert(std::end(multidim_meshgrid_vect),std::begin(tmp_cur_meshvect),std::end(tmp_cur_meshvect));
				tmp_cur_meshvect.clear();
			}	
		}	
		else 
		{
			for (size_t i=0; i<plt_coord_grid_1.size(); i++)	{
				for (size_t j=0; j<plt_coord_grid_2.size(); j++)	{
					multidim_meshgrid_vect.emplace_back(plt_coord_grid_1[i]);
					multidim_meshgrid_vect.emplace_back(plt_coord_grid_2[j]);
				}
			}
		}
		return multidim_meshgrid_vect;
	}
}


/************************************************************
 * Member function to generate meshgrid in desired dimension
 * Desired dimension to be obtained from size of input 
 * argument low_range_vect or up_range_vect or no_of_pts_vect
 * size of other_coord_val_vect is equal to the size of 
 * dynamics, if not throw error ........................
 * meshgrid_dyn_mode_id_vect gives mode numbers to be 
 * included in meshgrid preparation ....................
 * elements in vectors low_range_vect,up_range_vect,
 * no_of_pts_vect & meshgrid_dyn_mode_id_vect has 
 * one to one correspondence ...........................
 * slowest index to fastest as usual for multidimensional
 * meshgrid preparation .................................
 * **********************************************************/

std::vector<double> init_wp_val::generate_meshgrid_desired_dim(const std::vector<double> & low_range_vect,
		const std::vector<double> & up_range_vect, const std::vector<int> & no_of_pts_vect,
		const std::vector<double> & other_coord_val_vect,const std::vector<int> & meshgrid_dyn_mode_id_vect)	{
  	/******************************************************/
	/*** multidimensional vector for storing each one dimensional grid vectors 
	 *	out of which multidimensional meshgrid will be prepared **/
		std::vector<std::vector<double>> meshgrid_coords_multvect;
		for (size_t i=0; i<no_of_pts_vect.size(); i++)	{
			double delta_cur = (up_range_vect[i]-low_range_vect[i])/(no_of_pts_vect[i]-1); 
			std::vector<double> cur_onedim_grid_vect;
			for (int j=0; j<no_of_pts_vect[i]; j++)	{
				cur_onedim_grid_vect.emplace_back(low_range_vect[i]+j*delta_cur);
			}
			meshgrid_coords_multvect.push_back(cur_onedim_grid_vect);
			cur_onedim_grid_vect.clear();
		}
	/************************************************/
	std::vector<int> meshgrid_coord_index_id_vect;	
	// vector above will store locations of the modes (0-based)	//
	for (size_t i=0; i<meshgrid_dyn_mode_id_vect.size(); i++)	{
		int plt_coord_cur_id = std::distance(this->mode_id_arr, std::find(this-> mode_id_arr, this-> mode_id_arr+size, meshgrid_dyn_mode_id_vect[i]));
		meshgrid_coord_index_id_vect.emplace_back(plt_coord_cur_id);
	}
	// meshgrid_coord_index_id_vect is vector containing mode indices(0-based) involved in meshgrid //
	int simple_meshgrid_size  = std::accumulate(std::begin(no_of_pts_vect),std::end(no_of_pts_vect),1,std::multiplies<double>());
	int N_mesh = no_of_pts_vect.size();		// dimension of the meshgrid(direct product)
	/*** vector for returning multidimensional grid as one dimensional vector ***/
	std::vector<double>multidim_meshgrid_vect;
	/////////////////////////////
	/** if other coordinates are all zero, multidimensional meshgrid is just N-dimensional 
	 * direct product grid where N is size of no_of_pts_vect
	 * **/
	if (other_coord_val_vect.size()	!= INIT_STATE_CONST::NDIM_4_WP_INIT)	{
		throw "Size of other coordinate vector should be equal to dynamics dimension";
	}	
	else {
		std::vector<std::vector<double>> simple_meshgrid;
		std::vector<double> tmp_meshgrid;
		// calling stride array generator from fileop.hpp //
		std::vector<int> stride_vect_4_simple_mesh = GET_STRIDE_ARR_4_ANY(no_of_pts_vect); 
		for (int i=0; i<simple_meshgrid_size; i++)	{
			std::vector<int> tmp_stride_ind = get_specific_row_or_col_elements(stride_vect_4_simple_mesh,i,simple_meshgrid_size,N_mesh,"ROW"); 
			for (size_t j=0; j<tmp_stride_ind.size(); j++)	{
				tmp_meshgrid.emplace_back(meshgrid_coords_multvect[j][tmp_stride_ind[j]]);
			}
			simple_meshgrid.push_back(tmp_meshgrid);
			tmp_meshgrid.clear();
		}
		/////////////////////////
		std::vector<double> tmp_cur_meshvect;
		for (int i=0; i<simple_meshgrid_size; i++)	{
			tmp_cur_meshvect.insert(std::end(tmp_cur_meshvect),std::begin(other_coord_val_vect),std::end(other_coord_val_vect));
			for (int j=0; j<N_mesh; j++)	{
				tmp_cur_meshvect[meshgrid_coord_index_id_vect[j]] = simple_meshgrid[i][j];
			}
			multidim_meshgrid_vect.insert(std::end(multidim_meshgrid_vect),std::begin(tmp_cur_meshvect),std::end(tmp_cur_meshvect));
			tmp_cur_meshvect.clear();
		}
		return multidim_meshgrid_vect;
	}
}



/***********************************************************
 * Function for getting wave packet value at any queried 
 * location(multidimensional) ..........................
 * all_eigvect is vector storing all the relevant eigenvector
 * coefficients in serial order ........................
 * all_podvr_1dim_val_vect is a multidimensional vector
 * storing podvr function values for all relevant modes
 * (in serial order i.e. first podvr functions of first 
 * element in mode_id_arr is listed and so on	.......
 * each vector in this multidimensional vector stores 
 * all podvr basis function values for the mode at the 
 * queried grid ......................................
 * multidim_meshgrid_vect is multidimensional grid ..
 * length equal to INIT_STATE_CONST::NDIM_4_WP_INIT	.................
 * coeff_4_wp_tun_vect gives the contribution of 
 * different multidimensional eigenstates in the 
 * initial wave packet ................................
 * *******************************************************/

double init_wp_val::WP_init_val(std::vector<double> & all_eigvect, const std::vector<std::vector<double>> & all_podvr_1dim_val_vect,
		const std::vector<double> & multidim_meshgrid_vect, const std::vector<double> & coeff_4_wp_tun_vect)	{
/**************************************************************
 * multigrid_meshgrid_vect gives a particular multidimensional
 * grid out of the full meshgrid .. size of multigrid_meshgrid_vect
 * is equal to INIT_STATE_CONST::NDIM_4_WP_INIT	..
 * all_eigvect has all the coefficients for all the relevant 
 * eigenvectors for the wave packet ..
 * all_podvr_1dim_val_vect is a multidimensional vector 
 * dimensional equal to INIT_STATE_CONST::NDIM_4_WP_INIT	....
 * size of each element of all_podvr_1dim_val_vect is equal 
 * to no of podvr basis functions ... has one-to-one 
 * correspondence with mode_id_arr i.e. first Q1 basis ,
 * then Q5 and so on if mode_id_arr is =[1,5,...].......
 * *****************************************************/
	if (coeff_4_wp_tun_vect.size() > INIT_STATE_CONST::coeff_size_max)	{
		throw "Work with less number of basis functions for preparing initial state";
	}
	else	{
		/** vector for storing DPB values that is common for all the relevant eigenvectors **/
		std::vector<double> DPB_val_vect(this->DPB_size);
		/** vector for all onedim indices by caliing member function get_1_dim_index_vect **/
		std::vector<int> all_onedim_ids = this-> get_1_dim_index_vect();
		//////////////////////////////
		for (int i=0; i<this->DPB_size; i++)	{
			std::vector<int> cur_onedim_ids = get_specific_row_or_col_elements(all_onedim_ids,
					i,this->DPB_size,INIT_STATE_CONST::NDIM_4_WP_INIT,"ROW");
			double cur_dpb_val = 1;
			for (size_t j=0; j<cur_onedim_ids.size(); j++)	{
				cur_dpb_val *= all_podvr_1dim_val_vect[j][cur_onedim_ids[j]]; 
			}
			DPB_val_vect[i] = cur_dpb_val;
		}
		/////////////////////////////////
		/** vector for storing values of all relevant eigenstate values at the queried grid **/
		std::vector<double> eigvect_all_val_vect(coeff_4_wp_tun_vect.size());
		for (size_t i=0; i<coeff_4_wp_tun_vect.size(); i++)	{
			std::vector<double> cur_eigenvect = get_specific_row_or_col_elements(all_eigvect,i,coeff_4_wp_tun_vect.size(),this->DPB_size,"ROW");
			double cur_eigvect_val = 0.;
			for (int j=0; j<this->DPB_size; j++)	{
				cur_eigvect_val += cur_eigenvect[j]*DPB_val_vect[j];
			}
			eigvect_all_val_vect[i] = cur_eigvect_val;
		}
		///////////////////////////////
		/** variable for returning wave packet value at the queried grid **/
		double wp_cur_val = 0.;
		for (size_t i=0; i<coeff_4_wp_tun_vect.size(); i++)	{
			wp_cur_val += coeff_4_wp_tun_vect[i]*eigvect_all_val_vect[i];
		}
		return wp_cur_val;
	}
}

/**************************************************
 * class destructor for derived class init_wp_val *
 * ************************************************/

init_wp_val::~init_wp_val()	{
	delete[] this-> mode_id_arr;
	delete[] this-> Q_shift_arr;
	delete[] this-> n_row_so_2_anh_arr;
	delete[] this-> n_col_so_2_anh_arr;
	delete[] this-> n_row_anh_2_po_arr;
	delete[] this-> n_col_anh_2_po_arr;
	delete[] this-> sodvr_pts_arr;
	delete[] this-> so_2_anh_all_arr;
	delete[] this-> anh_2_po_all_arr;
}

