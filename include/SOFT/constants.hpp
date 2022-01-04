# if ! defined(CONST_H)
# define CONST_H 

# include <iosfwd>
# include <vector>
# include <string>
# include <complex>


/* List for fundamental constants .. 
 * source	........................ 
 * codata recommended values of the fundamental constants 2014 */ 

namespace FUNDAMENTAL_CONST	{
	inline constexpr double c = 29979245800;			// speed of light cm s^-1 
	inline constexpr double h_J = 6.626070040e-34;	// Planck constant in J-s
	inline constexpr double h_ev = 4.135667662e-15;	// Planck constant in ev-s
	inline constexpr double e  = 1.6021766208e-19;	// elementary charge in Coulomb
	inline constexpr double a0 = 0.52917721067e-10;  // Bohr radius in m 
	inline constexpr double H2J = 4.359744650e-18;   // Hartree to Joule conversion
	inline constexpr double H2ev = 27.21138602;      // hartree to ev 
	inline constexpr double m_e = 9.10938356e-31;	// electron mass in Kg 
	inline constexpr double k = 1.38064852e-23;      // Boltzmann constant in J K^-1
	inline constexpr double amu = 1.660539040e-27;   // atomic mass constant in Kg
	inline constexpr double N_A = 6.022140857e+23;   // Avogadro constant 
	inline constexpr double R = 8.3144598;           // Molar gas constant J K^-1 mol^-1
	inline constexpr double wn2ev = 1.239848738e-04;	// conversion factor for Wavenumber to eV
}


/* relevant parameters for preparing initial state */

namespace INIT_STATE_CONST	{
	inline constexpr int coeff_size_max = 6;		// no of eigenstates in the linear combination for generating the initial state
	inline constexpr int NDIM_4_WP_INIT = 3;		// No of coordnates involved to generate the initial state(some of the coordinate may have constant value in the entire grid)
	inline constexpr int NDIM_4_DYNAMICS = 3;		// dimension of the dynamics
	extern const std::string dyn_type; 
	extern const std::vector<int> eigstate_srno_wp_tun_vect;
	extern const std::vector<double> coeff_4_wp_tun_vect;	// coefficient vector for initial state generation
	extern const std::vector<int> basis_size_vect;		// vector containing basis size for different coordinates in the dynamics
	extern const std::vector<int> mode_id_vect;			// vector containing mode indices 
	extern const std::vector<double> Q_shift_vect;		
	extern const std::vector<std::string> dvr_type_vect;
	extern const std::vector<int> n_row_so_2_anh_vect;
	extern const std::vector<int> n_col_so_2_anh_vect;
	extern const std::vector<int> n_row_anh_2_po_vect;
	extern const std::vector<int> n_col_anh_2_po_vect;
	inline constexpr int tot_prim_dvr_pts = 165; 
	extern const std::string sodvr_pts_file;
	extern const std::string sodvr_2_anharm_file; 
	extern const std::string anharm_2_podvr_file;
	extern const std::string multidim_eigvect_filestring;
	extern const std::vector<int> IVR_init_quanta_vect;
	extern const std::vector<double> other_coord_val_vect;
	extern const std::vector<int> meshgrid_4_dyn_mode;
}

/* relevant parameters for carrying out the dynamics */

namespace DYN_CONST	{
	inline constexpr double timestep = 0.02;	// timestep for dynamics in femtosecond
	extern const std::vector<double> stepno_vect_4_anal;	
	inline constexpr int max_step = 2.5e+04;	// total number of steps for doing the dynamics
	inline constexpr double twopic = 2.0*M_PI*(FUNDAMENTAL_CONST::c)*1e-15;
	inline constexpr std::complex<double> EYE(0.,1.);
	inline constexpr int step_4_file = 25;	// in every step_4_file wave function will be saved to file(binary)  
	inline constexpr int N_save = 1001;		// no of times data to be saved 
	extern const std::string eigvalfile;		// file containing eigenvalues of the relevant subspace 
	extern const std::string eigstate_4_anal_dyn_file;	// file name without the state number having values in the grid
	extern const std::vector<int> grid_pt_vect;
	extern const std::vector<double> grid_low_lim_vect;	// contains lower limit of all the relevant coordinates to the dynamics
	extern const std::vector<double> grid_up_lim_vect;	// contains upper limit of all the relevant coordinates to the dynamics
	/*containing only relevant frequencies for the dynamics(magnitude of the imaginary frequency)*/
	extern const std::vector<double> freq_vect;
	extern const std::string initfiletype;
	extern const std::string potfile;
	extern const std::string initfile;
	extern const std::string quanta_4_overlap_file;
	inline constexpr int N_quanta = 560;
	inline constexpr double shiftval_tun_prob = 43;
	extern const std::string Q1_anhbasis_file;
	extern const std::string Q10_anhbasis_file;
	extern const std::string Q29_anhbasis_file;
	extern const std::string outfile_red_density_Q1;
	extern const std::string outfile_red_density_Q10;
	extern const std::string outfile_red_density_Q29;
	extern const std::string outfile_red_density_Q1_Q10;
	extern const std::string outfile_red_density_Q1_Q29;
	extern const std::string outfile_red_density_Q10_Q29;
	extern const std::string outfile_expect;
	extern const std::string outfile_overlap;
	extern const std::string outfile_closure;
	extern const std::string outfile_tun_prob; 
	extern const std::vector<double> stepsizevect;
}


# endif
