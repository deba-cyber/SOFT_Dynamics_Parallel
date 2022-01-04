
# include <SOFT/TROTTER_prop.hpp>
# include <SOFT/constants.hpp>
# include <SOFT/fileops.hpp>
# include <SOFT/data_process.hpp>
# include <iostream>
# include <fstream>
# include <string>
# include <cmath>
# include <vector>
# include <numeric>
# include <chrono>
# include <fftw3.h>
# include <omp.h>





int main()	{
	std::vector<int> dummyvectint;
	std::vector<double> dummyvectdouble;
	const int wp_size = std::accumulate(std::begin(
				DYN_CONST::grid_pt_vect),std::end(DYN_CONST::grid_pt_vect)
			,1,std::multiplies<int>());
	std::vector<fftw_complex> psi_x(wp_size);
	std::vector<int> line_extract_no_vect(wp_size);
	for (int i=0; i<wp_size; i++)	{
		line_extract_no_vect[i] = i+1;
	}
	/**------------------------------------**/
	/** Extracting data from potential and initial wave packet files **/
	std::vector<double> pot_vect;
	std::vector<double> psi_x_init_real;
	if (DYN_CONST::initfiletype == "dat")	{
		pot_vect = get_all_sp_lines(line_extract_no_vect,DYN_CONST::potfile);
		psi_x_init_real = get_all_sp_lines(line_extract_no_vect,DYN_CONST::initfile); 
		for (size_t i=0; i<psi_x.size(); i++)	{
			psi_x [i][0] = psi_x_init_real[i];
			psi_x [i][1] = 0.;
		}
	}
	else if (DYN_CONST::initfiletype == "bin")	{
		pot_vect = rtrn_vec_from_bin(dummyvectdouble,DYN_CONST::potfile,1,wp_size);
		psi_x_init_real = rtrn_vec_from_bin(dummyvectdouble,DYN_CONST::initfile,1,wp_size);
		for (size_t i=0; i<psi_x.size(); i++)	{
				psi_x [i][0] = psi_x_init_real[i];
				psi_x [i][1] = 0.;
		}
	}
	/**------------------------------------**/
	std::vector<std::vector<double>> p_mult_grid_vect; 
	/*----------------------------*/
	std::vector<double> delta_vect(DYN_CONST::grid_pt_vect.size());
	for (size_t i=0; i<delta_vect.size(); i++)	{
		delta_vect[i] = (DYN_CONST::grid_up_lim_vect[i]-
				DYN_CONST::grid_low_lim_vect[i])/(DYN_CONST::grid_pt_vect[i]-1);
	}
	/** product of all elements of delta_vect **/
	double multdelta = std::accumulate(delta_vect.begin(),delta_vect.end(),1.,std::multiplies<>());
	/*----------------------------*/
	for (size_t i=0; i<delta_vect.size(); i++)	{
		std::vector<double> p_grid_1d(DYN_CONST::grid_pt_vect[i]);
		onedim_momentum_grid(p_grid_1d,DYN_CONST::grid_pt_vect[i],delta_vect[i]);
		p_mult_grid_vect.push_back(p_grid_1d);
		empty_swap(p_grid_1d);
	}
	std::vector<double> stepsizevect_tilde = {p_mult_grid_vect[0][1]-p_mult_grid_vect[0][0],
		p_mult_grid_vect[1][1]-p_mult_grid_vect[1][0],p_mult_grid_vect[2][1]-p_mult_grid_vect[2][0]};
	double deltapdt_tilde = std::accumulate(stepsizevect_tilde.begin(),stepsizevect_tilde.end(),1.,std::multiplies<>());
	/*----------------------------*/
	std::vector<fftw_complex> trotter_pot = get_pot_trotter_part(pot_vect);
	std::vector<fftw_complex> trotter_ke = get_ke_trotter_part(p_mult_grid_vect,DYN_CONST::freq_vect);
	int sizevect[INIT_STATE_CONST::NDIM_4_DYNAMICS];
	int rank = INIT_STATE_CONST::NDIM_4_DYNAMICS; 
	for (int i=0; i<INIT_STATE_CONST::NDIM_4_DYNAMICS; i++)	{
		sizevect[i] = DYN_CONST::grid_pt_vect[i];
	}
	/** Parallel FFTW **/
	fftw_init_threads();
	int N_threads;
	#pragma omp parallel
	{
		#pragma omp single
		{
			N_threads = omp_get_max_threads();
		}
	}
	fftw_plan_with_nthreads(N_threads);
	fftw_plan plan_forward;
	fftw_plan plan_backward;
	fftw_plan plan_forward_out_of_place;
	plan_Multidim_ForwardFFT(plan_forward,psi_x.size(),sizevect,rank,"in place");	// in place forward FFT
	plan_Multidim_BackwardFFT(plan_backward,psi_x.size(),sizevect,rank,"in place");	// in place backward FFT
	plan_Multidim_ForwardFFT(plan_forward_out_of_place,psi_x.size(),sizevect,rank,"out of place");	// out of place forward FFT
	int dyn_counter = 0;	// counter for dynamics run
	int dyn_savecounter = 0;	// counter for dynamics data process
//	std::vector<double> psi_real_imag_abs_vect(2*wp_size);		// vector containing real,imaginary values of psi
	/**----------------------------------------**/
	/**----------------------------------------**/
	//=========================================================================//
	// relevant data for analysis that can be stored in advance //
	std::vector<int> linenovect_Q1_basis = Get_lineno_vect(INIT_STATE_CONST::basis_size_vect[0]);
	std::vector<int> linenovect_Q10_basis = Get_lineno_vect(INIT_STATE_CONST::basis_size_vect[1]);
	std::vector<int> linenovect_Q29_basis = Get_lineno_vect(INIT_STATE_CONST::basis_size_vect[2]);
	std::vector<std::vector<double>> Q1_anh_all = Get_all_lines_multvect(linenovect_Q1_basis,
			dummyvectdouble,DYN_CONST::Q1_anhbasis_file);
	std::vector<std::vector<double>> Q10_anh_all = Get_all_lines_multvect(linenovect_Q10_basis,
			dummyvectdouble,DYN_CONST::Q10_anhbasis_file);
	std::vector<std::vector<double>> Q29_anh_all = Get_all_lines_multvect(linenovect_Q29_basis,
			dummyvectdouble,DYN_CONST::Q29_anhbasis_file);
	// Grid normalization of 1D eigenstates of reduced Hamiltonians //
	for (unsigned int i=0; i<Q1_anh_all.size(); i++)	{
		double integral = onedim_simpson_template(Q1_anh_all[i],0.1);
		for_each(Q1_anh_all[i].begin(),Q1_anh_all[i].end(),[&](double& elem){elem/= sqrt(integral);});
	}
	for (unsigned int i=0; i<Q10_anh_all.size(); i++)	{
		double integral = onedim_simpson_template(Q10_anh_all[i],0.1);
		for_each(Q10_anh_all[i].begin(),Q10_anh_all[i].end(),[&](double& elem){elem/= sqrt(integral);});
	}
	for (unsigned int i=0; i<Q29_anh_all.size(); i++)	{
		double integral = onedim_simpson_template(Q29_anh_all[i],0.1);
		for_each(Q29_anh_all[i].begin(),Q29_anh_all[i].end(),[&](double& elem){elem/= sqrt(integral);});
	}
	// quanta combinations for overlap analysis //
	std::vector<int> linenovect_quanta = Get_lineno_vect(DYN_CONST::N_quanta);
	std::vector<std::vector<int>> quanta_vect = Get_all_lines_multvect(linenovect_quanta,
			dummyvectint,DYN_CONST::quanta_4_overlap_file); 
	/**
	rtrn_real_imag_norm(psi_real_imag_abs_vect,psi_x);
	save_binary(psi_real_imag_abs_vect,DYN_CONST::save_outbinfile,"append");
	psi_real_imag_abs_vect.assign(psi_real_imag_abs_vect.size(),0);
	**/
	/** Processing data from  initial state **/
	std::vector<int> grid_size_vect_Q1_Q10 = {87,121};
	std::vector<int> grid_size_vect_Q1_Q29 = {87,81};
	std::vector<int> grid_size_vect_Q10_Q29 = {121,81};
	std::vector<double> stepsizevect_IO = {0.1,0.1};
	double stepsize_IO = 0.1;
	std::vector<int> serialvect_IO_Q1_Q10 = {0,1};
	std::vector<int> serialvect_IO_Q1_Q29 = {0,2};
	std::vector<int> serialvect_IO_Q10_Q29 = {1,2};
	STRIDE_ARR_BACK_N_FORTH_DP stridedpfulldimobj(DYN_CONST::grid_pt_vect);
	STRIDE_ARR_BACK_N_FORTH_DP stridedp_Q1_Q10(grid_size_vect_Q1_Q10);
	STRIDE_ARR_BACK_N_FORTH_DP stridedp_Q1_Q29(grid_size_vect_Q1_Q29);
	STRIDE_ARR_BACK_N_FORTH_DP stridedp_Q10_Q29(grid_size_vect_Q10_Q29);
	std::vector<double> Q1_grid = generate_onedim_grid(DYN_CONST::grid_low_lim_vect[0],DYN_CONST::grid_up_lim_vect[0],
			DYN_CONST::grid_pt_vect[0]);
	std::vector<double> Q10_grid = generate_onedim_grid(DYN_CONST::grid_low_lim_vect[1],DYN_CONST::grid_up_lim_vect[1],
			DYN_CONST::grid_pt_vect[1]);
	std::vector<double> Q29_grid = generate_onedim_grid(DYN_CONST::grid_low_lim_vect[2],DYN_CONST::grid_up_lim_vect[2],
			DYN_CONST::grid_pt_vect[2]);
	std::vector<std::vector<double>> Q1_Q10_grid = generate_2d_meshgrid(Q1_grid,Q10_grid);
	std::vector<std::vector<double>> Q1_Q29_grid = generate_2d_meshgrid(Q1_grid,Q29_grid);
	std::vector<std::vector<double>> Q10_Q29_grid = generate_2d_meshgrid(Q10_grid,Q29_grid);
	std::vector<std::vector<double>> all_onedimgrids;
	all_onedimgrids.push_back(Q1_grid);
	all_onedimgrids.push_back(Q10_grid);
	all_onedimgrids.push_back(Q29_grid);
	int DPG_IO_size = 43*121*81;
	int DPG_IO_Q1_Q10 = 87*121;
	int DPG_IO_Q1_Q29 = 87*81;
	int DPG_IO_Q10_Q29 = 121*81;
	std::vector<int> grid_pt_vect_tun_prob = {43,121,81};
	STRIDE_ARR_BACK_N_FORTH_DP stridehalfgridobj(grid_pt_vect_tun_prob);
	//================================================================//
	double normalization = get_normalizn_fact(psi_x,stridedpfulldimobj,all_onedimgrids,DYN_CONST::stepsizevect);  
	std::vector<std::vector<double>> expect_all_vect;
	std::vector<std::vector<double>> overlap_vect;	
	std::vector<double> closure_vect(DYN_CONST::N_save);
	std::vector<double> tun_prob_vect(DYN_CONST::N_save);
	std::vector<std::vector<double>> Q1_red_density_vect;
	std::vector<std::vector<double>> Q10_red_density_vect;
	std::vector<std::vector<double>> Q29_red_density_vect;
	std::vector<std::vector<double>> Q1_Q10_red_density_vect;
	std::vector<std::vector<double>> Q1_Q29_red_density_vect;
	std::vector<std::vector<double>> Q10_Q29_red_density_vect;
	/// --- Tunneling probability calculation ---- ///
	tun_prob_vect[dyn_savecounter] =  get_tun_prob(psi_x,stridedpfulldimobj,stridehalfgridobj,DPG_IO_size,
			DYN_CONST::shiftval_tun_prob,DYN_CONST::stepsizevect);
	///--- Expectation calculation -----///
	std::vector<double> expect_cur_vect(12);
	std::vector<double> expect_pos_sqr = get_expect_pos_sqr(psi_x,stridedpfulldimobj,all_onedimgrids,
			DYN_CONST::stepsizevect,normalization);
	std::vector<double> expect_pos = get_expect_pos_or_momentum(psi_x,stridedpfulldimobj,all_onedimgrids,
			DYN_CONST::stepsizevect,normalization);
	expect_cur_vect[0] = expect_pos_sqr[0];
	expect_cur_vect[1] = expect_pos_sqr[1];
	expect_cur_vect[2] = expect_pos_sqr[2];
	expect_cur_vect[3] = expect_pos[0];
	expect_cur_vect[4] = expect_pos[1];
	expect_cur_vect[5] = expect_pos[2];
	std::vector<double> p_sqr_expect_vect(INIT_STATE_CONST::NDIM_4_DYNAMICS);
	std::vector<double> p_expect_vect(INIT_STATE_CONST::NDIM_4_DYNAMICS);
	std::vector<fftw_complex> wp_fftw_cur_tilde(psi_x.size());
	execute_outofplaceFFT(psi_x,wp_fftw_cur_tilde,"FORWARD",plan_forward_out_of_place);
	for_each(wp_fftw_cur_tilde.begin(),wp_fftw_cur_tilde.end(),[&multdelta](fftw_complex& elem){
			elem[0] *= multdelta;
			elem[1] *= multdelta;});
	double p_overlap = 0.;		// momemtum space overlap  
	for_each(wp_fftw_cur_tilde.begin(),wp_fftw_cur_tilde.end(),[&p_overlap,&deltapdt_tilde](fftw_complex& elem){
			p_overlap += deltapdt_tilde*(pow(elem[0],2.)+pow(elem[1],2.));});
	{
		for (unsigned int i=0; i<p_sqr_expect_vect.size(); i++)	{
			double p_sqr_expect = 0.;
			double p_expect = 0.;
			int counter_p_sqr_expect = 0;
			for_each(wp_fftw_cur_tilde.begin(),wp_fftw_cur_tilde.end(),[&](fftw_complex& elem){
					std::vector<int> cur_id_vect = stridedpfulldimobj.onedim_indices(counter_p_sqr_expect+1);
					p_sqr_expect += deltapdt_tilde*(pow(elem[0],2.)+pow(elem[1],2.))*(pow(p_mult_grid_vect[i][cur_id_vect[i]],2.)); 
					p_expect += deltapdt_tilde*(pow(elem[0],2.)+pow(elem[1],2.))*p_mult_grid_vect[i][cur_id_vect[i]];
					counter_p_sqr_expect++;});
			p_sqr_expect_vect[i] = p_sqr_expect/p_overlap;
			p_expect_vect[i] = p_expect/p_overlap;
		}	
	}
	expect_cur_vect[6] = p_sqr_expect_vect[0];
	expect_cur_vect[7] = p_sqr_expect_vect[1];
	expect_cur_vect[8] = p_sqr_expect_vect[2];
	expect_cur_vect[9] = p_expect_vect[0];
	expect_cur_vect[10] = p_expect_vect[1];
	expect_cur_vect[11] = p_expect_vect[2];
	expect_all_vect.push_back(expect_cur_vect);
	expect_cur_vect.assign(expect_cur_vect.size(),0.);
	p_sqr_expect_vect.assign(p_sqr_expect_vect.size(),0.);
	p_expect_vect.assign(p_expect_vect.size(),0.);
	for_each(wp_fftw_cur_tilde.begin(),wp_fftw_cur_tilde.end(),[&wp_fftw_cur_tilde](fftw_complex& elem){
			elem[0] = 0.;
			elem[1] = 0.;});
	//--- Reduced probability density calculation -----//
	int id_coord_IO_Q1 = 0;
	int id_coord_IO_Q10 = 1;
	int id_coord_IO_Q29 = 2;
	int counter_Q1_fixed = 0;
	int counter_Q10_fixed = 0; 
	int counter_Q29_fixed = 0;
	int counter_Q1_Q10_fixed = 0;
	int counter_Q1_Q29_fixed = 0;
	int counter_Q10_Q29_fixed = 0;
	std::vector<double> red_density_Q1_cur_vect(Q1_grid);
	std::vector<double> red_density_Q10_cur_vect(Q10_grid);
	std::vector<double> red_density_Q29_cur_vect(Q29_grid);
	std::vector<double> red_density_Q1_Q10_cur_vect(DPG_IO_Q1_Q10);
	std::vector<double> red_density_Q1_Q29_cur_vect(DPG_IO_Q1_Q29);
	std::vector<double> red_density_Q10_Q29_cur_vect(DPG_IO_Q10_Q29);
	for_each(Q1_grid.begin(),Q1_grid.end(),[&](double& elem){std::vector<int> query_grid_id = {counter_Q1_fixed,0,0};
			red_density_Q1_cur_vect[counter_Q1_fixed] = get_reduced_prob_density_IO_multcoord(psi_x,query_grid_id,serialvect_IO_Q10_Q29,
					stridedp_Q10_Q29,stridedpfulldimobj,DPG_IO_Q10_Q29,stepsizevect_IO);
			counter_Q1_fixed++;});
	Q1_red_density_vect.push_back(red_density_Q1_cur_vect);
	for_each(Q10_grid.begin(),Q10_grid.end(),[&](double& elem){std::vector<int> query_grid_id = {0,counter_Q10_fixed,0};
			red_density_Q10_cur_vect[counter_Q10_fixed] = get_reduced_prob_density_IO_multcoord(psi_x,query_grid_id,serialvect_IO_Q1_Q29,
					stridedp_Q1_Q29,stridedpfulldimobj,DPG_IO_Q1_Q29,stepsizevect_IO);
			counter_Q10_fixed++;});
	Q10_red_density_vect.push_back(red_density_Q10_cur_vect);
	for_each(Q29_grid.begin(),Q29_grid.end(),[&](double& elem){std::vector<int> query_grid_id = {0,0,counter_Q29_fixed};
			red_density_Q29_cur_vect[counter_Q29_fixed] =  get_reduced_prob_density_IO_multcoord(psi_x,query_grid_id,serialvect_IO_Q1_Q10,
					stridedp_Q1_Q10,stridedpfulldimobj,DPG_IO_Q1_Q10,stepsizevect_IO);
			counter_Q29_fixed++;});
	Q29_red_density_vect.push_back(red_density_Q29_cur_vect);
	for_each(Q1_Q10_grid.begin(),Q1_Q10_grid.end(),[&](std::vector<double>& elem){
			std::vector<int> cur_id_vect = stridedp_Q1_Q10.onedim_indices(counter_Q1_Q10_fixed+1); 
			std::vector<int> query_multgrid_id = {cur_id_vect[0],cur_id_vect[1],0};
			red_density_Q1_Q10_cur_vect[counter_Q1_Q10_fixed] = get_reduced_prob_density_IO_coord(psi_x,query_multgrid_id,id_coord_IO_Q29,
					stridedpfulldimobj,Q29_grid.size(),stepsize_IO); 
			counter_Q1_Q10_fixed++;});
	Q1_Q10_red_density_vect.push_back(red_density_Q1_Q10_cur_vect);
	for_each(Q1_Q29_grid.begin(),Q1_Q29_grid.end(),[&](std::vector<double>& elem){
			std::vector<int> cur_id_vect = stridedp_Q1_Q29.onedim_indices(counter_Q1_Q29_fixed+1);
			std::vector<int> query_multgrid_id = {cur_id_vect[0],0,cur_id_vect[1]};
			red_density_Q1_Q29_cur_vect[counter_Q1_Q29_fixed] = get_reduced_prob_density_IO_coord(psi_x,query_multgrid_id,id_coord_IO_Q10,
					stridedpfulldimobj,Q10_grid.size(),stepsize_IO);
			counter_Q1_Q29_fixed++;});
	Q1_Q29_red_density_vect.push_back(red_density_Q1_Q29_cur_vect);
	for_each(Q10_Q29_grid.begin(),Q10_Q29_grid.end(),[&](std::vector<double>& elem){
			std::vector<int> cur_id_vect = stridedp_Q10_Q29.onedim_indices(counter_Q10_Q29_fixed+1);
			std::vector<int> query_multgrid_id = {0,cur_id_vect[0],cur_id_vect[1]};
			red_density_Q10_Q29_cur_vect[counter_Q10_Q29_fixed] = get_reduced_prob_density_IO_coord(psi_x,query_multgrid_id,id_coord_IO_Q1,
					stridedpfulldimobj,Q1_grid.size(),stepsize_IO);
			counter_Q10_Q29_fixed++;});
	Q10_Q29_red_density_vect.push_back(red_density_Q10_Q29_cur_vect);
	red_density_Q1_cur_vect.assign(Q1_grid.size(),0.);
	red_density_Q10_cur_vect.assign(Q10_grid.size(),0.);
	red_density_Q29_cur_vect.assign(Q29_grid.size(),0.);
	red_density_Q1_Q10_cur_vect.assign(Q1_Q10_grid.size(),0.);
	red_density_Q1_Q29_cur_vect.assign(Q1_Q29_grid.size(),0.);
	red_density_Q10_Q29_cur_vect.assign(Q10_Q29_grid.size(),0.);
	//--- Overlap calculation -----//
	double chk_closure = 0.;
	std::vector<double> overlap_cur_vect(2*quanta_vect.size());
#pragma omp parallel for default(none) shared(psi_x,Q1_anh_all,Q10_anh_all,Q29_anh_all,\
		stridedpfulldimobj,overlap_cur_vect,quanta_vect) reduction(+: chk_closure) schedule(static)
	for (unsigned int i=0; i<quanta_vect.size(); i++)	{
		std::vector<fftw_complex> overlap_cur = get_overlap_w_dp(psi_x,Q1_anh_all[quanta_vect[i][0]-1],
				Q10_anh_all[quanta_vect[i][1]-1],Q29_anh_all[quanta_vect[i][2]-1],
				stridedpfulldimobj,DYN_CONST::stepsizevect);
		overlap_cur_vect[2*i] = overlap_cur[0][0];
		overlap_cur_vect[2*i + 1] = overlap_cur[0][1];
		chk_closure += pow(overlap_cur[0][0],2.)+pow(overlap_cur[0][1],2.);
	}
	overlap_vect.push_back(overlap_cur_vect);
	overlap_cur_vect.assign(2*quanta_vect.size(),0.);
	closure_vect[dyn_savecounter] = chk_closure;
	chk_closure = 0.;
	counter_Q1_fixed = 0;
	counter_Q10_fixed = 0; 
	counter_Q29_fixed = 0;
	counter_Q1_Q10_fixed = 0;
	counter_Q1_Q29_fixed = 0;
	counter_Q10_Q29_fixed = 0;
	/**----------------------------------------**/
	/**----------------------------------------**/
	while (dyn_counter < DYN_CONST::max_step+1)	{
		generate_fftw_complex_product(trotter_pot,psi_x);
		execute_inplaceFFT(psi_x,"FORWARD",plan_forward);		
		for_each(psi_x.begin(),psi_x.end(),[&multdelta](fftw_complex& elem){
				elem[0] *= multdelta;
				elem[1] *= multdelta;});
		generate_fftw_complex_product(trotter_ke,psi_x);
		execute_inplaceFFT(psi_x,"BACKWARD",plan_backward);
		for_each(psi_x.begin(),psi_x.end(),[&multdelta](fftw_complex& elem){
				elem[0] *= (1/multdelta);
				elem[1] *= (1/multdelta);});
		generate_fftw_complex_product(trotter_pot,psi_x);
		dyn_counter ++;
		if (dyn_counter% DYN_CONST::step_4_file == 0)	{
			dyn_savecounter++;
		    expect_pos_sqr = get_expect_pos_sqr(psi_x,stridedpfulldimobj,all_onedimgrids,
					DYN_CONST::stepsizevect,normalization);
			expect_pos = get_expect_pos_or_momentum(psi_x,stridedpfulldimobj,all_onedimgrids,
					DYN_CONST::stepsizevect,normalization);
			expect_cur_vect[0] = expect_pos_sqr[0];
			expect_cur_vect[1] = expect_pos_sqr[1];
			expect_cur_vect[2] = expect_pos_sqr[2];
			expect_cur_vect[3] = expect_pos[0];
			expect_cur_vect[4] = expect_pos[1];
			expect_cur_vect[5] = expect_pos[2];
			execute_outofplaceFFT(psi_x,wp_fftw_cur_tilde,"FORWARD",plan_forward_out_of_place);
			for_each(wp_fftw_cur_tilde.begin(),wp_fftw_cur_tilde.end(),[&multdelta](fftw_complex& elem){
					elem[0] *= multdelta;
					elem[1] *= multdelta;});
			p_overlap = 0.;		// momemtum space overlap  
			for_each(wp_fftw_cur_tilde.begin(),wp_fftw_cur_tilde.end(),[&p_overlap,&deltapdt_tilde](fftw_complex& elem){
					p_overlap += deltapdt_tilde*(pow(elem[0],2.)+pow(elem[1],2.));});
			{
				for (unsigned int i=0; i<p_sqr_expect_vect.size(); i++)	{
					double p_sqr_expect = 0.;
					double p_expect = 0.;
					int counter_p_sqr_expect = 0;
					for_each(wp_fftw_cur_tilde.begin(),wp_fftw_cur_tilde.end(),[&](fftw_complex& elem){
							std::vector<int> cur_id_vect = stridedpfulldimobj.onedim_indices(counter_p_sqr_expect+1);
							p_sqr_expect += deltapdt_tilde*(pow(elem[0],2.)+pow(elem[1],2.))*(pow(p_mult_grid_vect[i][cur_id_vect[i]],2.)); 
							p_expect += deltapdt_tilde*(pow(elem[0],2.)+pow(elem[1],2.))*p_mult_grid_vect[i][cur_id_vect[i]];
							counter_p_sqr_expect++;});
					p_sqr_expect_vect[i] = p_sqr_expect/p_overlap;
					p_expect_vect[i] = p_expect/p_overlap;
				}	
			}
			expect_cur_vect[6] = p_sqr_expect_vect[0];
			expect_cur_vect[7] = p_sqr_expect_vect[1];
			expect_cur_vect[8] = p_sqr_expect_vect[2];
			expect_cur_vect[9] = p_expect_vect[0];
			expect_cur_vect[10] = p_expect_vect[1];
			expect_cur_vect[11] = p_expect_vect[2];
			expect_all_vect.push_back(expect_cur_vect);
			tun_prob_vect[dyn_savecounter] =  get_tun_prob(psi_x,stridedpfulldimobj,stridehalfgridobj,DPG_IO_size,
					DYN_CONST::shiftval_tun_prob,DYN_CONST::stepsizevect);
			for_each(Q1_grid.begin(),Q1_grid.end(),[&](double& elem){std::vector<int> query_grid_id = {counter_Q1_fixed,0,0};
					red_density_Q1_cur_vect[counter_Q1_fixed] = get_reduced_prob_density_IO_multcoord(psi_x,query_grid_id,serialvect_IO_Q10_Q29,
							stridedp_Q10_Q29,stridedpfulldimobj,DPG_IO_Q10_Q29,stepsizevect_IO);
					counter_Q1_fixed++;});
			Q1_red_density_vect.push_back(red_density_Q1_cur_vect);
			for_each(Q10_grid.begin(),Q10_grid.end(),[&](double& elem){std::vector<int> query_grid_id = {0,counter_Q10_fixed,0};
					red_density_Q10_cur_vect[counter_Q10_fixed] = get_reduced_prob_density_IO_multcoord(psi_x,query_grid_id,serialvect_IO_Q1_Q29,
							stridedp_Q1_Q29,stridedpfulldimobj,DPG_IO_Q1_Q29,stepsizevect_IO);
					counter_Q10_fixed++;});
			Q10_red_density_vect.push_back(red_density_Q10_cur_vect);
			for_each(Q29_grid.begin(),Q29_grid.end(),[&](double& elem){std::vector<int> query_grid_id = {0,0,counter_Q29_fixed};
					red_density_Q29_cur_vect[counter_Q29_fixed] =  get_reduced_prob_density_IO_multcoord(psi_x,query_grid_id,serialvect_IO_Q1_Q10,
							stridedp_Q1_Q10,stridedpfulldimobj,DPG_IO_Q1_Q10,stepsizevect_IO);
					counter_Q29_fixed++;});
			Q29_red_density_vect.push_back(red_density_Q29_cur_vect);
			for_each(Q1_Q10_grid.begin(),Q1_Q10_grid.end(),[&](std::vector<double>& elem){
					std::vector<int> cur_id_vect = stridedp_Q1_Q10.onedim_indices(counter_Q1_Q10_fixed+1); 
					std::vector<int> query_multgrid_id = {cur_id_vect[0],cur_id_vect[1],0};
					red_density_Q1_Q10_cur_vect[counter_Q1_Q10_fixed] = get_reduced_prob_density_IO_coord(psi_x,query_multgrid_id,id_coord_IO_Q29,
							stridedpfulldimobj,Q29_grid.size(),stepsize_IO); 
					counter_Q1_Q10_fixed++;});
			Q1_Q10_red_density_vect.push_back(red_density_Q1_Q10_cur_vect);
			for_each(Q1_Q29_grid.begin(),Q1_Q29_grid.end(),[&](std::vector<double>& elem){
					std::vector<int> cur_id_vect = stridedp_Q1_Q29.onedim_indices(counter_Q1_Q29_fixed+1);
					std::vector<int> query_multgrid_id = {cur_id_vect[0],0,cur_id_vect[1]};
					red_density_Q1_Q29_cur_vect[counter_Q1_Q29_fixed] = get_reduced_prob_density_IO_coord(psi_x,query_multgrid_id,id_coord_IO_Q10,
							stridedpfulldimobj,Q10_grid.size(),stepsize_IO);
					counter_Q1_Q29_fixed++;});
			Q1_Q29_red_density_vect.push_back(red_density_Q1_Q29_cur_vect);
			for_each(Q10_Q29_grid.begin(),Q10_Q29_grid.end(),[&](std::vector<double>& elem){
					std::vector<int> cur_id_vect = stridedp_Q10_Q29.onedim_indices(counter_Q10_Q29_fixed+1);
					std::vector<int> query_multgrid_id = {0,cur_id_vect[0],cur_id_vect[1]};
					red_density_Q10_Q29_cur_vect[counter_Q10_Q29_fixed] = get_reduced_prob_density_IO_coord(psi_x,query_multgrid_id,id_coord_IO_Q1,
							stridedpfulldimobj,Q1_grid.size(),stepsize_IO);
					counter_Q10_Q29_fixed++;});
			Q10_Q29_red_density_vect.push_back(red_density_Q10_Q29_cur_vect);
			red_density_Q1_cur_vect.assign(Q1_grid.size(),0.);
			red_density_Q10_cur_vect.assign(Q10_grid.size(),0.);
			red_density_Q29_cur_vect.assign(Q29_grid.size(),0.);
			red_density_Q1_Q10_cur_vect.assign(Q1_Q10_grid.size(),0.);
			red_density_Q1_Q29_cur_vect.assign(Q1_Q29_grid.size(),0.);
			red_density_Q10_Q29_cur_vect.assign(Q10_Q29_grid.size(),0.);
			counter_Q1_fixed = 0;
			counter_Q10_fixed = 0;
			counter_Q29_fixed = 0;
			counter_Q1_Q10_fixed = 0;
			counter_Q1_Q29_fixed = 0;
			counter_Q10_Q29_fixed = 0;
			expect_cur_vect.assign(expect_cur_vect.size(),0.);
			p_sqr_expect_vect.assign(p_sqr_expect_vect.size(),0.);
			p_expect_vect.assign(p_expect_vect.size(),0.);
			for_each(wp_fftw_cur_tilde.begin(),wp_fftw_cur_tilde.end(),[&wp_fftw_cur_tilde](fftw_complex& elem){
					elem[0] = 0.;
					elem[1] = 0.;});
#pragma omp parallel for default(none) shared(psi_x,Q1_anh_all,Q10_anh_all,Q29_anh_all,\
		stridedpfulldimobj,overlap_cur_vect,quanta_vect) reduction(+: chk_closure) schedule(static)
			for (unsigned int i=0; i<quanta_vect.size(); i++)	{
				std::vector<fftw_complex> overlap_cur = get_overlap_w_dp(psi_x,Q1_anh_all[quanta_vect[i][0]-1],
						Q10_anh_all[quanta_vect[i][1]-1],Q29_anh_all[quanta_vect[i][2]-1],
						stridedpfulldimobj,DYN_CONST::stepsizevect);
				overlap_cur_vect[2*i] = overlap_cur[0][0];
				overlap_cur_vect[2*i + 1] = overlap_cur[0][1];
				chk_closure += pow(overlap_cur[0][0],2.)+pow(overlap_cur[0][1],2.);
			}
			overlap_vect.push_back(overlap_cur_vect);
			overlap_cur_vect.assign(2*quanta_vect.size(),0.);
			closure_vect[dyn_savecounter] = chk_closure;
			chk_closure = 0.;
//			wp_fftw_cur_tilde.assign(wp_fftw_cur_tilde.size(),0.);
//			save_binary(tun_prob_vect,DYN_CONST::outfile_tun_prob,"append");
//			std::cout	<< dyn_counter	<<	"	"	<<	dyn_savecounter	<< "\n";
//			rtrn_real_imag_norm(psi_real_imag_abs_vect,psi_x);
//			save_binary(psi_real_imag_abs_vect,DYN_CONST::save_outbinfile,"append");
//			psi_real_imag_abs_vect.assign(psi_real_imag_abs_vect.size(),0);
		}
	}
	save_binary(tun_prob_vect,DYN_CONST::outfile_tun_prob,"append");
	save_twodim_vect_binary(Q1_red_density_vect,DYN_CONST::outfile_red_density_Q1,"append");
	save_twodim_vect_binary(Q10_red_density_vect,DYN_CONST::outfile_red_density_Q10,"append");
	save_twodim_vect_binary(Q29_red_density_vect,DYN_CONST::outfile_red_density_Q29,"append");
	save_twodim_vect_binary(Q1_Q10_red_density_vect,DYN_CONST::outfile_red_density_Q1_Q10,"append");
	save_twodim_vect_binary(Q1_Q29_red_density_vect,DYN_CONST::outfile_red_density_Q1_Q29,"append");
	save_twodim_vect_binary(Q10_Q29_red_density_vect,DYN_CONST::outfile_red_density_Q10_Q29,"append");
	save_twodim_vect_binary(expect_all_vect,DYN_CONST::outfile_expect,"append");
	save_twodim_vect_binary(overlap_vect,DYN_CONST::outfile_overlap,"append");
	save_binary(closure_vect,DYN_CONST::outfile_closure,"append");
	return EXIT_SUCCESS;
}
