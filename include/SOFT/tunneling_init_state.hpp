# if !defined (TUNN_INIT)
# define TUNN_INIT

# include <iosfwd>
# include <string>
# include <vector>


////////////////////////////////////////////////////////////////////////////////////////////

/***** stride array class for calculating stride array
 * *** generating one dimensional index array from multidimensional index
 * *** generating multidimensional index for direct product basis from one 
 * *** dimensional index array *****/


class STRIDE_ARR_BACK_N_FORTH	{
	protected:
		int* basis_size_arr;
		int size;
	public:
		//======== Constructor ==========//
		STRIDE_ARR_BACK_N_FORTH(int* basis_size_arr,int size);	
		//=========== Member Functions ===========//
		std::vector<int> stride_arr();
		int multidim_index_dpb(const std::vector<int> &onedim_index_vect);
		std::vector<int> onedim_indices(const int & multidim_index);
		//========	Destructor	===========//
		~STRIDE_ARR_BACK_N_FORTH();
};

/*------------------------------------------------------------*/

/*** Following class will generate value of the wavefunction when queried at a multidimensional 
 * *	grid point.. Inherits attributes from STRIDE_ARR_BACK_N_FORTH	****/


/* ONLY ONE DIMENSIONAL ARRAY WILL BE USED */

class init_wp_val : public STRIDE_ARR_BACK_N_FORTH {
	private:
		int DPB_size;
		int* mode_id_arr;
		double* Q_shift_arr;
		int tot_prim_dvr_pts;						// integer .. total no of primitive dvr points combining all coordinates
		std::vector<std::string> dvr_type_vect;
		int* n_row_so_2_anh_arr;					// array containing no of rows for each so_2_anh transformation matrix for each coordinate
		int* n_col_so_2_anh_arr;					// array containing no of columns for each so_2_anh transformation matrix for each coordinate
		int* n_row_anh_2_po_arr;					// array containing no of rows for each anh_2_po transformation matrix for each coordinate
		int* n_col_anh_2_po_arr;					// array containing no of columns for each anh_2_po transformation matrix for each coordinate
		int tot_size_so_2_anh_all;					// total no of elements including all so_2_anh transformation matrices in one dimensional array
		int tot_size_anh_2_po_all;					// total no of elements including all anh_2_po trandformation matrices in one dimensional array
		double* so_2_anh_all_arr;					// 3-dimensional array as one dimensional array 	
		double* anh_2_po_all_arr;					// 3-dimensional array as one dimensional array
		double* sodvr_pts_arr;						// sodvr_pts_arr will be needed for sinc-dvr values, SODVR points for all coordinates in one array 
		std::string dyn_type; 
	public:
		//==== constructor ====//
		init_wp_val(int* basis_size_arr, int size, int DPB_size,
		int* mode_id_arr, const std::vector<std::string>& dvr_type_vect, double* Q_shift_arr, int tot_prim_dvr_pts,
		int* n_row_so_2_anh_arr, int* n_col_so_2_anh_arr, int* n_row_anh_2_po_arr, int* n_col_anh_2_po_arr,
		int tot_size_so_2_anh_all, int tot_size_anh_2_po_all,
		double* so_2_anh_all_arr, double* anh_2_po_all_arr, double* sodvr_pts_arr, 
		const std::string& dyn_type);	
	//========= Member Functions ========//
		std::vector<int>get_1_dim_index_vect();
		std::vector<double> prim_basis_val_vect(int N_HO_basis, double Q_cur_val);
		std::vector<double>dvr_basis_val_vect(int mode_id, double cur_coord_val);
		double pobasis_val(int mode_id, double cur_coord_val, int cur_coord_pobasis_id);
		std::vector<double> generate_multidim_meshgrid(const std::vector<double>& plt_coord_grid_1, const std::vector<double>& plt_coord_grid_2,
				const std::vector<double>& other_coord_val_vect, int plt_coord_1, int plt_coord_2);
		std::vector<double> generate_meshgrid_desired_dim(const std::vector<double>& low_range_vect,
				const std::vector<double>& up_range_vect, const std::vector<int>& no_of_pts_vect,
				const std::vector<double>& other_coord_val_vect,const std::vector<int>& meshgrid_dyn_mode_id_vect);
		double WP_init_val(std::vector<double>& all_eigvect, const std::vector<std::vector<double>>& all_podvr_1dim_val_vect,
				const std::vector<double>& multidim_meshgrid_vect, const std::vector<double>& coeff_4_wp_tun_vect);
	//======= destructor =============//
		~init_wp_val();
};


# endif
