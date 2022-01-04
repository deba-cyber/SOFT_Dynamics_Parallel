
# include <SOFT/TROTTER_prop.hpp>
# include <SOFT/constants.hpp>
# include <SOFT/fileops.hpp>
# include <iostream>
# include <cmath>
# include <numeric>
# include <vector>

/*-------------------------------------------------------------*/
/* function to generate potential part of Trotter product 
 * of type fftw_complex */

std::vector<fftw_complex> get_pot_trotter_part(const std::vector<double>& V_vect)	{ 
	std::vector<fftw_complex> pot_trotter_part(V_vect.size());
	for (size_t i=0; i<V_vect.size(); i++)	{
		pot_trotter_part[i][0] = exp(-0.5*DYN_CONST::EYE*V_vect[i]*DYN_CONST::twopic*DYN_CONST::timestep).real(); 
		pot_trotter_part[i][1] = exp(-0.5*DYN_CONST::EYE*V_vect[i]*DYN_CONST::twopic*DYN_CONST::timestep).imag();
	}
	return pot_trotter_part;
}

/*------------------------------------------------------------*/

/* function to generate ke part of Trotter product of type 
 * fftw_complex */

std::vector<fftw_complex> get_ke_trotter_part(std::vector<std::vector<double>>& p_grid_vect,
		const std::vector<double>& freq_vect)	{
	int size = std::accumulate(std::begin(DYN_CONST::grid_pt_vect),
			std::end(DYN_CONST::grid_pt_vect),1,std::multiplies<int>());
	std::vector<fftw_complex> ke_trotter_part(size); 
	std::vector<int> stride_vect_4_dp_grid = GET_STRIDE_ARR_4_ANY(DYN_CONST::grid_pt_vect);
	for (int i=0; i<size; i++)	{
		std::vector<int> tmp_stride_cur_id_vect = get_specific_row_or_col_elements(stride_vect_4_dp_grid,i,
				size,freq_vect.size(),"ROW");
		double cur_exponent_cnt = 0.;
		for (size_t j=0; j<tmp_stride_cur_id_vect.size(); j++)	{
			cur_exponent_cnt += freq_vect[j]*pow(p_grid_vect[j][tmp_stride_cur_id_vect[j]],2.);
		}
		ke_trotter_part[i][0] = exp(-0.5*DYN_CONST::EYE*cur_exponent_cnt*
				DYN_CONST::twopic*DYN_CONST::timestep).real();
		ke_trotter_part[i][1] = exp(-0.5*DYN_CONST::EYE*cur_exponent_cnt*
				DYN_CONST::twopic*DYN_CONST::timestep).imag();
	}
	return ke_trotter_part;
}


/*-------------------------------------------------------------*/

/* function for multiplying two fftw_complex vectors 
 * psi_cur_vect  is to be modified in place 
 * after multiplying with the first vector */

void generate_fftw_complex_product(std::vector<fftw_complex>& trotter_act,
		std::vector<fftw_complex>& psi_cur_vect)	{
	double tmpreal,tmpimag;
	for (size_t i=0; i<psi_cur_vect.size(); i++)	{
		tmpreal = trotter_act[i][0]*psi_cur_vect[i][0] - trotter_act[i][1]*psi_cur_vect[i][1];
		tmpimag = trotter_act[i][0]*psi_cur_vect[i][1] + trotter_act[i][1]*psi_cur_vect[i][0];
		psi_cur_vect[i][0] = tmpreal;
		psi_cur_vect[i][1] = tmpimag;
	}
}

/*-------------------------------------------------------------*/

/* function for generating momentum grid from real space grid 
 * size & real space grid spacing **/

void onedim_momentum_grid(std::vector<double>& p_grid,int N_pts,
		double delta)	{
	p_grid[0] = 0.;
	if (N_pts %2 != 0)	{
		for (int i=1; i<((N_pts-1)/2+1); i++)	{
			p_grid[i] = i*2.*M_PI*(1/(N_pts*delta));
			p_grid[N_pts-i] = -1.*i*2.*M_PI*(1/(N_pts*delta));
		}
	}
	else if (N_pts%2 == 0)	{
		int halfpoint = N_pts/2;
		p_grid[halfpoint] = -1.*(N_pts/2)*2.*M_PI*(1/(N_pts*delta));
		for (int i=1; i<(N_pts/2); i++)	{
			p_grid [i] = i*2.*M_PI*(1/(N_pts*delta));
			p_grid [N_pts-i] = -1.*i*2.*M_PI*(1/(N_pts*delta));
		}
	}
}
		

/*-------------------------------------------------------------*/

/* function for shifting DC component to the center for plotting
 * purpose in Fourier domain ..**/

void fftwshift_4_plot(std::vector<fftw_complex>& inplotvect,
		std::vector<fftw_complex>& invect)	{
	const int size = invect.size();
	const int halfsize = int(size/2);
	if (size%2 != 0)	{
		inplotvect [halfsize][0] = invect [0][0];
		inplotvect [halfsize][1] = invect [0][1];
		for (int i=0; i<halfsize; i++)	{
			inplotvect [i][0] = invect[halfsize+i+1][0];
			inplotvect [i][1] = invect[halfsize+i+1][1];
			inplotvect [halfsize+i+1][0] = invect[i+1][0];
			inplotvect [halfsize+i+1][1] = invect[i+1][1];
		}
	}
	else if (size%2 == 0)	{
		for (int i=0; i<size; i++)	{
			if (i >= halfsize)	{
				inplotvect [i][0] = invect[i-halfsize][0];
				inplotvect [i][1] = invect[i-halfsize][1];
			}
			else	{
				inplotvect [i][0] = invect[i+halfsize][0];
				inplotvect [i][1] = invect[i+halfsize][1];
			}
		}
	}
}


/** function for shifting DC component of momentum grid to 
 * the center for plotting purpose **/

void fftwshift_momentum_grid_4_plot(std::vector<double>& p_grid,
		std::vector<double>& p_grid_4_plot)	{
	const int size = p_grid.size();
	const int halfsize = int(size/2);
	if (size%2 != 0)	{
		p_grid_4_plot [halfsize] = p_grid[0];
		for (int i=0; i<halfsize; i++)	{
			p_grid_4_plot [i] = p_grid[halfsize+i+1];
			p_grid_4_plot [halfsize+i+1] = p_grid[i+1];
		}
	}
	else if (size%2 == 0)	{
		for (int i=0; i<size; i++)	{
			if (i >= halfsize)	{
				p_grid_4_plot[i] = p_grid[i-halfsize];
			}
			else	{
				p_grid_4_plot[i] = p_grid[i+halfsize];
			}
		}
	}
}
			
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/

/* complex multidimensional transform using basic interface */

/***************** plan forward ***********************/

void plan_Multidim_ForwardFFT(fftw_plan& plan, int N_pts, 
		const int* sizevect,int rank, 
		const std::string& type)	{
	if (type == "in place")	{
		std::vector<fftw_complex> invect(N_pts);
		plan = fftw_plan_dft(rank,sizevect,invect.data(),
				invect.data(),FFTW_FORWARD,FFTW_MEASURE);
	}
	else if (type == "out of place")	{
		std::vector<fftw_complex> invect(N_pts);
		std::vector<fftw_complex> outvect(N_pts);
		plan = fftw_plan_dft(rank,sizevect,invect.data(),
				outvect.data(),FFTW_FORWARD,FFTW_MEASURE);
	}
}

/***************** plan backward ***********************/

void plan_Multidim_BackwardFFT(fftw_plan& plan, int N_pts, 
		const int* sizevect,int rank, 
		const std::string& type)	{
	if (type == "in place")	{
		std::vector<fftw_complex> invect(N_pts);
		plan = fftw_plan_dft(rank,sizevect,invect.data(),
				invect.data(),FFTW_BACKWARD,FFTW_MEASURE);
	}
	else if (type == "out of place")	{
		std::vector<fftw_complex> invect(N_pts);
		std::vector<fftw_complex> outvect(N_pts);
		plan = fftw_plan_dft(rank,sizevect,invect.data(),
				outvect.data(),FFTW_BACKWARD,FFTW_MEASURE);
	}
}

/*-------------------------------------------------------------*/

/* complex one dimensional FFT using basic interface */

/***************** plan forward ***********************/

void plan_ForwardFFT_1d(fftw_plan& plan,int N_pts,
		const std::string& type)	{
	if (type == "in place")	{
		std::vector<fftw_complex> invect(N_pts);
		plan = fftw_plan_dft_1d(N_pts, invect.data(), 
				invect.data(),FFTW_FORWARD,FFTW_MEASURE);
	}
	else if (type == "out of place")	{
		std::vector<fftw_complex> invect(N_pts);
		std::vector<fftw_complex> outvect(N_pts);
		plan = fftw_plan_dft_1d(N_pts, invect.data(),
				outvect.data(),FFTW_FORWARD,FFTW_MEASURE);
	}
}

/***************** plan backward ***********************/

void plan_BackwardFFT_1d(fftw_plan& plan,int N_pts,
		const std::string& type)	{
	if (type == "in place")	{
		std::vector<fftw_complex> invect(N_pts);
		plan = fftw_plan_dft_1d(N_pts, invect.data(), 
				invect.data(),FFTW_BACKWARD,FFTW_MEASURE);
	}
	else if (type == "out of place")	{
		std::vector<fftw_complex> invect(N_pts);
		std::vector<fftw_complex> outvect(N_pts);
		plan = fftw_plan_dft_1d(N_pts, invect.data(),
				outvect.data(),FFTW_BACKWARD,FFTW_MEASURE);
	}
}

/****************** Execute FFT **********************/

void execute_inplaceFFT(std::vector<fftw_complex>& in,
		const std::string& direction,fftw_plan plan)	{
	fftw_execute_dft(plan, in.data(), in.data());
	if (direction == "BACKWARD")	{
	/* normalization */
		for (size_t i=0; i<in.size(); i++)	{
			in[i][0] *= (1/static_cast<double>(in.size()));
			in[i][1] *= (1/static_cast<double>(in.size()));
		}
	}
}


void execute_outofplaceFFT(std::vector<fftw_complex>& in,
		std::vector<fftw_complex>& out,
		const std::string& direction,fftw_plan plan)	{
	fftw_execute_dft(plan, in.data(), out.data());
	if (direction == "BACKWARD")	{
	/* normalization */
		for (size_t i=0; i<out.size(); i++)	{
			out[i][0] *= (1/static_cast<double>(out.size()));
			out[i][1] *= (1/static_cast<double>(out.size()));
		}
	}
}

/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/

/****************************************
 * function for returning real,imaginary
 * parts and norm of each element of the 
 * vector passed ........................
 * returning a 1-dimensional vector .....
 * size 2 times the size of passed vector
 * in order real part, imaginary part ..
 * Norm can be calculated using real 
 * and imaginary parts ...............
 * ***************************************/

void rtrn_real_imag_norm(std::vector<double>& rtrn_real_imag_nrm_vect,
		const std::vector<fftw_complex>& psi_cur_vect)	{
	int cur_id = 0;
	for (size_t i=0; i<psi_cur_vect.size(); i++)	{
		rtrn_real_imag_nrm_vect[cur_id] = psi_cur_vect[i][0];
		rtrn_real_imag_nrm_vect[cur_id+1] = psi_cur_vect[i][1];
		cur_id += 2;
	}
}


