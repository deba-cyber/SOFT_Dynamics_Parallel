# include <SOFT/tunneling_init_state.hpp>
# include <SOFT/constants.hpp>
# include <SOFT/fileops.hpp>
# include <iostream>
# include <vector>
# include <algorithm>
# include <string>
# include <numeric>
# include <chrono>
# include <omp.h>



int main()	{
	/*** Necessary parameters ****/
	int size = INIT_STATE_CONST::NDIM_4_WP_INIT;
	int basis_size_arr[size];
	int mode_id_arr[size];
	double Q_shift_arr[size];
	int n_row_so_2_anh_arr[size];
	int n_col_so_2_anh_arr[size];
	int n_row_anh_2_po_arr[size];
	int n_col_anh_2_po_arr[size];
	for (int i=0; i<size; i++)	{
		basis_size_arr[i] = INIT_STATE_CONST::basis_size_vect[i];
		mode_id_arr[i] = INIT_STATE_CONST::mode_id_vect[i];
		Q_shift_arr[i] = INIT_STATE_CONST::Q_shift_vect[i];
		n_row_so_2_anh_arr[i] = INIT_STATE_CONST::n_row_so_2_anh_vect[i];
		n_col_so_2_anh_arr[i] = INIT_STATE_CONST::n_col_so_2_anh_vect[i];
		n_row_anh_2_po_arr[i] = INIT_STATE_CONST::n_row_anh_2_po_vect[i];
		n_col_anh_2_po_arr[i] = INIT_STATE_CONST::n_col_anh_2_po_vect[i];	
	}
	int DPB_size = std::accumulate(basis_size_arr,basis_size_arr+size, 1, std::multiplies<int>());
	int tot_prim_dvr_pts = INIT_STATE_CONST::tot_prim_dvr_pts ; 
	/* total no of elements including all transformation matrices so=>anh */
	int tot_size_so_2_anh_all = 0;
	/* total no of elements including all transformation matrices anh=>po */
	int tot_size_anh_2_po_all = 0;
	for (int i=0; i<size; i++)	{
		tot_size_so_2_anh_all += n_row_so_2_anh_arr[i]*n_col_so_2_anh_arr[i];
		tot_size_anh_2_po_all += n_row_anh_2_po_arr[i]*n_col_anh_2_po_arr[i];
	}
	/**** collecting data from input files *****/
	std::vector<double> sodvr_pts_all = vect_accumulate_all_col(mode_id_arr, INIT_STATE_CONST::sodvr_pts_file, size); 
	std::vector<double> sodvr_2_anh_all = vect_accumulate_all_col(mode_id_arr, INIT_STATE_CONST::sodvr_2_anharm_file, size);
	std::vector<double> anh_2_po_all = vect_accumulate_all_col(mode_id_arr, INIT_STATE_CONST::anharm_2_podvr_file, size);
	double sodvr_pts_arr[sodvr_pts_all.size()];
	double so_2_anh_all_arr[tot_size_so_2_anh_all];
	double anh_2_po_all_arr[tot_size_anh_2_po_all];
	std::copy(sodvr_pts_all.begin(), sodvr_pts_all.end(), sodvr_pts_arr);
	std::copy(sodvr_2_anh_all.begin(), sodvr_2_anh_all.end(), so_2_anh_all_arr);
	std::copy(anh_2_po_all.begin(), anh_2_po_all.end(), anh_2_po_all_arr);
	empty_swap(sodvr_pts_all);
	empty_swap(sodvr_2_anh_all);
	empty_swap(anh_2_po_all);
	/*---------------------------------------------*/
	std::string dyn_type = INIT_STATE_CONST::dyn_type;
	/*---------------------------------------------*/
	//============ creating init_wp_val object  ==============// 
	init_wp_val INIT_WP(basis_size_arr,size,DPB_size,mode_id_arr,INIT_STATE_CONST::dvr_type_vect,Q_shift_arr,
			tot_prim_dvr_pts,n_row_so_2_anh_arr,n_col_so_2_anh_arr,n_row_anh_2_po_arr,n_col_anh_2_po_arr,
			tot_size_so_2_anh_all,tot_size_anh_2_po_all,so_2_anh_all_arr,anh_2_po_all_arr,sodvr_pts_arr,
			INIT_STATE_CONST::dyn_type);
	/**----------------------------------------------*/
	/*---- getting wave packet values ----*/
	//------------ preparing meshgrid --------------//
	std::vector<double> gh1 = INIT_WP.generate_meshgrid_desired_dim(DYN_CONST::grid_low_lim_vect,DYN_CONST::grid_up_lim_vect,
	DYN_CONST::grid_pt_vect,INIT_STATE_CONST::other_coord_val_vect,INIT_STATE_CONST::meshgrid_4_dyn_mode);
	/**----------------------------------------------*/
	/**----------------------------------------------*/
	/**----------------------------------------------*/
	/* Preparing 2-dimensional vector containing one dimensional grid vectors for all 
	 * relevant coordinates in the calculation 
	 * */
	const int DPG_size = std::accumulate(std::begin(DYN_CONST::grid_pt_vect),
			std::end(DYN_CONST::grid_pt_vect),1,std::multiplies<int>());
	std::vector<std::vector<double>> all_onedim_grid_vect;
	if (DYN_CONST::grid_pt_vect.size() != INIT_STATE_CONST::NDIM_4_WP_INIT)	{
		int counter = 0;
		for (int i=0; i<INIT_STATE_CONST::NDIM_4_WP_INIT; i++)	{
			int location = std::find(INIT_STATE_CONST::meshgrid_4_dyn_mode.begin(),
					INIT_STATE_CONST::meshgrid_4_dyn_mode.end(),INIT_STATE_CONST::mode_id_vect[i])!=INIT_STATE_CONST::meshgrid_4_dyn_mode.end();
			if (location == 0)	{
				all_onedim_grid_vect.push_back({INIT_STATE_CONST::other_coord_val_vect[i]});
			}
			else if (location != 0)	{
				std::vector<double> cur_1d_tmp_vect(DYN_CONST::grid_pt_vect[counter]);
				double delta_cur = (DYN_CONST::grid_up_lim_vect[counter]-DYN_CONST::grid_low_lim_vect[counter])/(
						DYN_CONST::grid_pt_vect[counter]-1);
				for (size_t j=0; j<cur_1d_tmp_vect.size(); j++)	{
					cur_1d_tmp_vect[j] = (DYN_CONST::grid_low_lim_vect[counter]+j*delta_cur);
				}
				all_onedim_grid_vect.push_back(cur_1d_tmp_vect);
				counter++;
			}
		}
	}
	else if (DYN_CONST::grid_pt_vect.size() == INIT_STATE_CONST::NDIM_4_WP_INIT)	{
		for (int i=0; i<INIT_STATE_CONST::NDIM_4_WP_INIT; i++)	{
			std::vector<double> cur_1d_tmp_vect(DYN_CONST::grid_pt_vect[i]);
			double delta_cur = (DYN_CONST::grid_up_lim_vect[i]-DYN_CONST::grid_low_lim_vect[i])/(
					DYN_CONST::grid_pt_vect[i]-1); 
			for (size_t j=0; j<cur_1d_tmp_vect.size(); j++)	{
				cur_1d_tmp_vect[j] = (DYN_CONST::grid_low_lim_vect[i]+j*delta_cur);
			}
			all_onedim_grid_vect.push_back(cur_1d_tmp_vect);
			cur_1d_tmp_vect.clear();
		}
	}
	/**----------------------------------------------*/
	/**----------------------------------------------*/
	/**----------------------------------------------*/
	if (dyn_type == "TUNNELING")	{
		std::vector<double> all_eigvec_coeffs = get_all_sp_lines(INIT_STATE_CONST::eigstate_srno_wp_tun_vect,
				INIT_STATE_CONST::multidim_eigvect_filestring);
		std::vector<std::vector<double>> all_pobasis_4_tun_vect;
		for (int i=0; i<INIT_STATE_CONST::NDIM_4_WP_INIT; i++)	{
			std::vector<double> tmp_cur_coord_pobasis_vect;
			for (int j=0; j<INIT_STATE_CONST::basis_size_vect[i]; j++)	{
				for (size_t k=0; k<all_onedim_grid_vect[i].size(); k++)	{	
					double pobasis_val_cur = INIT_WP.pobasis_val(INIT_STATE_CONST::mode_id_vect[i],
							all_onedim_grid_vect[i][k],j);
					tmp_cur_coord_pobasis_vect.emplace_back(pobasis_val_cur);
				}
			}
			all_pobasis_4_tun_vect.push_back(tmp_cur_coord_pobasis_vect);
			tmp_cur_coord_pobasis_vect.clear();
		}
		/**----------------------------**/
		std::vector<std::vector<double>> cur_mesh_all;
		std::vector<std::vector<std::vector<double>>> all_podvr_1dim_val_vect_tot;
		for (int i=0; i<DPG_size; i++)	{
			/** current multidimensional grid **/
			std::vector<double> cur_mesh = get_specific_row_or_col_elements(gh1,i,DPG_size,INIT_STATE_CONST::NDIM_4_WP_INIT,"ROW");
			/** multidimensional vector for storing all pobasis values for all coords at the given grid **/
			std::vector<std::vector<double>> all_podvr_1dim_val_vect;
			for (size_t j=0; j<cur_mesh.size(); j++)	{
				auto it = std::find(all_onedim_grid_vect[j].begin(),all_onedim_grid_vect[j].end(),cur_mesh[j]);
				int mesh_id_4_cur_coord = std::distance(all_onedim_grid_vect[j].begin(),it);
				std::vector<double> podvr_cur_coord_vect = get_specific_row_or_col_elements(all_pobasis_4_tun_vect[j],
						mesh_id_4_cur_coord,INIT_STATE_CONST::basis_size_vect[j],all_onedim_grid_vect[j].size(),"COLUMN");
				all_podvr_1dim_val_vect.push_back(podvr_cur_coord_vect);
			}
			cur_mesh_all.push_back(cur_mesh);
			all_podvr_1dim_val_vect_tot.push_back(all_podvr_1dim_val_vect);
		}
		/**----------------------------**/
		std::vector<double> all_wp_init_val_vect(DPG_size);	// vector for storing wp initial values
		auto t1 = std::chrono::high_resolution_clock::now();
		#pragma omp parallel for default(none) shared(all_wp_init_val_vect,all_eigvec_coeffs,\
				all_podvr_1dim_val_vect_tot,cur_mesh_all,INIT_WP) schedule(static)
		for (int i=0; i<DPG_size; i++)	{
			double wp_cur_val = INIT_WP.WP_init_val(all_eigvec_coeffs,all_podvr_1dim_val_vect_tot[i],
					cur_mesh_all[i],INIT_STATE_CONST::coeff_4_wp_tun_vect);
			all_wp_init_val_vect[i] = wp_cur_val;
		}
		auto t2 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count(); 
		std::cout	<< "Time taken : "	<< duration	<< " microseconds "	<< "\n";
		save_binary(all_wp_init_val_vect,DYN_CONST::initfile,"append");
	}
	/**----------------------------------------------*/
	/**----------------------------------------------*/
	/**----------------------------------------------*/
	else if (dyn_type == "IVR")	{
		std::vector<std::vector<double>> all_anhbasis_4_ivr_vect;
		std::vector<std::vector<double>> anhbasis_4_Q1_pair_vect;
		for (int i=0; i<INIT_STATE_CONST::NDIM_4_WP_INIT; i++)	{
			std::vector<double> tmp_cur_coord_anh_vect;
			if (INIT_STATE_CONST::mode_id_vect[i] == 1)	{
				std::vector<double> tmp_cur_coord_anh_pair_vect;
				for (size_t j=0; j<all_onedim_grid_vect[i].size(); j++)	{
					double anh_val_cur = INIT_WP.pobasis_val(INIT_STATE_CONST::mode_id_vect[i],
							all_onedim_grid_vect[i][j],INIT_STATE_CONST::IVR_init_quanta_vect[i]);
					double anh_val_cur_pair = INIT_WP.pobasis_val(INIT_STATE_CONST::mode_id_vect[i],
							all_onedim_grid_vect[i][j],INIT_STATE_CONST::IVR_init_quanta_vect[i]+1);
					tmp_cur_coord_anh_vect.emplace_back(anh_val_cur);
					tmp_cur_coord_anh_pair_vect.emplace_back(anh_val_cur_pair);
				}
				all_anhbasis_4_ivr_vect.push_back(tmp_cur_coord_anh_vect);
				anhbasis_4_Q1_pair_vect.push_back(tmp_cur_coord_anh_pair_vect);
			}
			else {
				for (size_t j=0; j<all_onedim_grid_vect[i].size(); j++)	{
					double anh_val_cur = INIT_WP.pobasis_val(INIT_STATE_CONST::mode_id_vect[i],
							all_onedim_grid_vect[i][j],INIT_STATE_CONST::IVR_init_quanta_vect[i]);
					tmp_cur_coord_anh_vect.emplace_back(anh_val_cur);
				}
				all_anhbasis_4_ivr_vect.push_back(tmp_cur_coord_anh_vect);
			}
		}
		/**--------------------------------------**/
		std::vector<double> wp_init_ivr_vect(DPG_size);
		#pragma omp parallel for schedule(static) default(none) shared(all_onedim_grid_vect,\
				all_anhbasis_4_ivr_vect,anhbasis_4_Q1_pair_vect,gh1,wp_init_ivr_vect)
		for (int i=0; i<DPG_size; i++)	{
			std::vector<double> cur_mesh = get_specific_row_or_col_elements(gh1,i,DPG_size,
					INIT_STATE_CONST::NDIM_4_WP_INIT,"ROW"); 
			double wp_init_ivr_cur = 1.;
			double wp_init_ivr_cur_pair = 1.;
			auto it1 = std::find(all_onedim_grid_vect[0].begin(),all_onedim_grid_vect[0].end(),cur_mesh[0]);
			int Q1_mesh_id = std::distance(all_onedim_grid_vect[0].begin(),it1);
			wp_init_ivr_cur *= all_anhbasis_4_ivr_vect[0][Q1_mesh_id];
			wp_init_ivr_cur_pair *= anhbasis_4_Q1_pair_vect[0][Q1_mesh_id];
			int j;
			for (j=1; j<cur_mesh.size(); j++)	{
				auto it = std::find(all_onedim_grid_vect[j].begin(),all_onedim_grid_vect[j].end(),cur_mesh[j]);
				int mesh_id_4_cur_coord = std::distance(all_onedim_grid_vect[j].begin(),it); 
				wp_init_ivr_cur *= all_anhbasis_4_ivr_vect[j][mesh_id_4_cur_coord];
				wp_init_ivr_cur_pair *= all_anhbasis_4_ivr_vect[j][mesh_id_4_cur_coord];
			}
			wp_init_ivr_vect[i] = ((1./sqrt(2)*(wp_init_ivr_cur+wp_init_ivr_cur_pair)));
		}
		save_binary(wp_init_ivr_vect,DYN_CONST::initfile,"append");
	}
	return EXIT_SUCCESS;
}


