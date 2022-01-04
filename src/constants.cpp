# include <SOFT/constants.hpp>
# include <iosfwd>
# include <string>
# include <cmath>
# include <vector>

namespace INIT_STATE_CONST	{
	extern const std::string dyn_type = "IVR";
	extern const std::vector<int> eigstate_srno_wp_tun_vect = {8,11};
	extern const std::vector<double> coeff_4_wp_tun_vect = {1/sqrt(2.),1/sqrt(2.)};
	extern const std::vector<int> basis_size_vect = {18,17,7};
	extern const std::vector<int> mode_id_vect = {1,10,29};			// vector containing mode indices 
	extern const std::vector<double> Q_shift_vect = {0.,1.680548,0.};		
	extern const std::vector<std::string> dvr_type_vect = {"sinc","HO","sinc"};
	extern const std::vector<int> n_row_so_2_anh_vect = {18,17,7};
	extern const std::vector<int> n_col_so_2_anh_vect = {54,35,76};
	extern const std::vector<int> n_row_anh_2_po_vect = {18,17,7};
	extern const std::vector<int> n_col_anh_2_po_vect = {18,17,7};
	extern const std::string sodvr_pts_file = "../input_data/sodvr_pts_";
	extern const std::string sodvr_2_anharm_file = "../input_data/sodvr_2_anharm_"; 
	extern const std::string anharm_2_podvr_file = "../input_data/anharm_2_podvr_";
	extern const std::string multidim_eigvect_filestring = "../input_data/Eigenvectors_1_10_29_3d_18_17_7_size";
	extern const std::vector<int> IVR_init_quanta_vect = {2,0,0};
	extern const std::vector<double> other_coord_val_vect = {0.,0.,0.};
	extern const std::vector<int> meshgrid_4_dyn_mode = {1,10,29};
}

namespace DYN_CONST {
	extern const std::string eigvalfile = "../input_data/Eigvals_1_10_twodim_18_17_size";
	extern const std::string eigstate_4_anal_dyn_file = "../input_data/eigstate_4_anal_dyn_";
	extern const std::vector<double> stepno_vect_4_anal = {50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};	
	extern const std::vector<int> grid_pt_vect = {87,121,81};
	extern const std::vector<double> grid_low_lim_vect = {-4.3,-4.,-4.};
	extern const std::vector<double> grid_up_lim_vect = {4.3,8.,4.};
	extern const std::vector<double> freq_vect = {1215.7032,675.7,2058.2};
	extern const std::string initfiletype = "bin";	// whether initial potential and wp files are in binary/dat format
	extern const std::string potfile = "../input_data/Potential_3D_extended_bound_chk";		// file containing potential data
	extern const std::string initfile = "../input_data/Q1_Q10_Q29_init_case14";			// file containing initial wave packet data
	extern const std::string quanta_4_overlap_file = "../input_data/3D_quanta_arr_dyn";
	extern const std::string Q1_anhbasis_file = "../input_data/Q1_anhbasis_val_3D_dyn";
	extern const std::string Q10_anhbasis_file = "../input_data/Q10_anhbasis_val_3D_dyn";
	extern const std::string Q29_anhbasis_file = "../input_data/Q29_anhbasis_val_3D_dyn";
	extern const std::string outfile_red_density_Q1 = "../output_data/Q1_Q10_Q29_red_density_Q1_case14";
	extern const std::string outfile_red_density_Q10 = "../output_data/Q1_Q10_Q29_red_density_Q10_case14";
	extern const std::string outfile_red_density_Q29 = "../output_data/Q1_Q10_Q29_red_density_Q29_case14";
	extern const std::string outfile_red_density_Q1_Q10 = "../output_data/Q1_Q10_Q29_red_density_Q1_Q10_case14";
	extern const std::string outfile_red_density_Q1_Q29 = "../output_data/Q1_Q10_Q29_red_density_Q1_Q29_case14";
	extern const std::string outfile_red_density_Q10_Q29 = "../output_data/Q1_Q10_Q29_red_density_Q10_Q29_case14";
	extern const std::string outfile_expect = "../output_data/Q1_Q10_Q29_expect_all_case14";
	extern const std::string outfile_overlap = "../output_data/Q1_Q10_Q29_overlap_case14";
	extern const std::string outfile_closure = "../output_data/Q1_Q10_Q29_closure_case14";
	extern const std::string outfile_tun_prob = "../output_data/Q1_Q10_Q29_tun_prob_case14";
	extern const std::vector<double> stepsizevect = {0.1,0.1,0.1};
}
