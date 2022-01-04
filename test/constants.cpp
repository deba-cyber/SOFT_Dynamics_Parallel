# include <SOFT/constants.hpp>
# include <iosfwd>
# include <string>
# include <cmath>
# include <vector>

namespace INIT_STATE_CONST	{
	extern const std::string dyn_type = "TUNNELING";
	extern const std::vector<int> eigstate_srno_wp_tun_vect = {1,2};
	extern const std::vector<double> coeff_4_wp_tun_vect = {1/sqrt(2.),1/sqrt(2.)};
	extern const std::vector<int> basis_size_vect = {18,17};
	extern const std::vector<int> mode_id_vect = {1,10};			// vector containing mode indices 
	extern const std::vector<double> Q_shift_vect = {0.,1.680548};		
	extern const std::vector<std::string> dvr_type_vect = {"sinc","HO"};
	extern const std::vector<int> n_row_so_2_anh_vect = {18,17};
	extern const std::vector<int> n_col_so_2_anh_vect = {54,35};
	extern const std::vector<int> n_row_anh_2_po_vect = {18,17};
	extern const std::vector<int> n_col_anh_2_po_vect = {18,17};
	extern const std::string sodvr_pts_file = "../input_data/sodvr_pts_";
	extern const std::string sodvr_2_anharm_file = "../input_data/sodvr_2_anharm_"; 
	extern const std::string anharm_2_podvr_file = "../input_data/anharm_2_podvr_";
	extern const std::string multidim_eigvect_filestring = "../input_data/Eigenvectors_1_10_2d_18_17_size";
	extern const std::vector<int> IVR_init_quanta_vect = {0,0};
	extern const std::vector<double> other_coord_val_vect = {0.,0.};
	extern const std::vector<int> meshgrid_4_dyn_mode = {1,10};
}

namespace DYN_CONST {
	extern const std::string eigvalfile = "../input_data/Eigvals_1_10_twodim_18_17_size";
	extern const std::string eigstate_4_anal_dyn_file = "../input_data/eigstate_4_anal_dyn_";
	extern const std::vector<double> stepno_vect_4_anal = {50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};	
	extern const std::vector<int> grid_pt_vect = {101,101};
	extern const std::vector<double> grid_low_lim_vect = {-4.,-4.};
	extern const std::vector<double> grid_up_lim_vect = {4.,4.};
	extern const std::vector<double> freq_vect = {1215.7032,675.7};
	extern const std::string initfiletype = "bin";	// whether initial potential and wp files are in binary/dat format
	extern const std::string potfile = "../input_data/Potential_1_10_for_wp_dyn";		// file containing potential data
	extern const std::string initfile = "../input_data/Q1_Q10_1_2_wp_sym";			// file containing initial wave packet data
	extern const std::string restart_outbinfile = "../output_data/psi_dyn_restart_Q1_Q10_1_2_tun";
	extern const std::string save_outbinfile = "../output_data/psi_dyn_Q1_Q10_1_2_tun";		// binary file name for saving 
	extern const std::string anal_outbinfile_abs = "../output_data/psi_dyn_anal_abs";	// binary file name for saving analytical results
	extern const std::string anal_outbinfile_real = "../output_data/psi_dyn_anal_real";	// binary file name for saving analytical results
	extern const std::string anal_outbinfile_imag = "../output_data/psi_dyn_anal_imag";	// binary file name for saving analytical results
}
