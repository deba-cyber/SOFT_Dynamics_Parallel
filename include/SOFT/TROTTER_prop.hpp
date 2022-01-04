# if !defined (TROTTER)
# define TROTTER

# include <SOFT/constants.hpp>
# include <iosfwd>
# include <vector>
# include <string>
# include <fftw3.h> 

/* function to generate potential part of Trotter product 
 * of type fftw_complex */

std::vector<fftw_complex> get_pot_trotter_part(const std::vector<double>& V_vect); 

/*------------------------------------------------------------*/

/* function to generate ke part of Trotter product of type 
 * fftw_complex */

std::vector<fftw_complex> get_ke_trotter_part(std::vector<std::vector<double>>& p_grid_vect,
		const std::vector<double>& freq_vect);


/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/
/** Throughout the propagation, wave function will be modified 
 *	in place ..
 *	In each step of propagation ..
 *	first real space wave function will be acted on by potential 
 *	part of the Trotter product .. in place 
 *	then it would be Fourier Transformed .. in place 
 *  momentum space wave function will be acted on by momentum part 
 *  of the Trotter product .. in place 
 *  then it would be back Fourier Transformed to real space .. in place **/
/*-------------------------------------------------------------*/


/* function for multiplying two fftw_complex vectors 
 * psi_cur_vect  is to be modified in place 
 * after multiplying with the first vector */

void generate_fftw_complex_product(std::vector<fftw_complex>& trotter_act,
		std::vector<fftw_complex>& psi_cur_vect);

/*-------------------------------------------------------------*/

/* function for generating momentum grid from real space grid 
 * size & real space grid spacing **/

void onedim_momentum_grid(std::vector<double>& p_grid,int N_pts,
		double delta);

/*-------------------------------------------------------------*/

/* function for shifting DC component to the center for plotting
 * purpose in Fourier domain ..**/

void fftwshift_4_plot(std::vector<fftw_complex>& inplotvect,
		std::vector<fftw_complex>& invect);

/*-------------------------------------------------------------*/

/** function for shifting DC component of momentum grid to 
 * the center for plotting purpose **/

void fftwshift_momentum_grid_4_plot(std::vector<double>& p_grid,
		std::vector<double>& p_grid_4_plot);

/*-------------------------------------------------------------*/

/* complex multidimensional transform using basic interface */

/***************** plan forward ***********************/

void plan_Multidim_ForwardFFT(fftw_plan& plan, int N_pts, 
		const int* sizevect,int rank, const std::string& type);

/***************** plan backward ***********************/

void plan_Multidim_BackwardFFT(fftw_plan& plan, int N_pts, 
		const int* sizevect,int rank, const std::string& type);


/*-------------------------------------------------------------*/

/* complex one dimensional FFT using basic interface */

/***************** plan forward ***********************/

void plan_ForwardFFT_1d(fftw_plan& plan,int N_pts,
		const std::string& type);

/***************** plan backward ***********************/

void plan_BackwardFFT_1d(fftw_plan& plan,int N_pts,
		const std::string& type);

/****************** Execute FFT **********************/

void execute_inplaceFFT(std::vector<fftw_complex>& in,
		const std::string& direction,fftw_plan plan);

void execute_outofplaceFFT(std::vector<fftw_complex>& in,
		std::vector<fftw_complex>& out,
		const std::string& direction,fftw_plan plan);

/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/

/* function to return real,imaginary parts & absolute value of 
 * the wave function from <fftw_complex> type */
		
void rtrn_real_imag_norm(std::vector<double>& rtrn_real_imag_nrm_vect,
		const std::vector<fftw_complex>& psi_cur_vect);

/*--------------------------------------------------------------*/

# endif
