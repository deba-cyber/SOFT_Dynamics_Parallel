
# include <SOFT/constants.hpp>
# include <SOFT/fileops.hpp>
# include <iostream>
# include <cmath>
# include <string>
# include <vector>
# include <complex>
# include <numeric>


int main()	{
	const int wp_size = std::accumulate(std::begin(
				DYN_CONST::grid_pt_vect),std::end(DYN_CONST::grid_pt_vect),
			1,std::multiplies<int>());
	/***------------------------------------------------***/
	std::vector<int> line_extract_no_vect(wp_size);
	for (int i=0; i<wp_size; i++)	{
		line_extract_no_vect[i] = i+1;
	}
	/*** storing values of all eigenstates in the grid ***/
	/***------------------------------------------------***/
	std::vector<std::vector<double>> eigstates_4_wp;
	std::vector<double> dummyvect;
	for (size_t i=0; i<INIT_STATE_CONST::coeff_4_wp_tun_vect.size(); i++)	{
		std::string curfilename = DYN_CONST::eigstate_4_anal_dyn_file+std::to_string(
				INIT_STATE_CONST::eigstate_srno_wp_tun_vect[i]); 
		std::vector<double> cur_eigstate =  rtrn_vec_from_bin(dummyvect,curfilename,1,wp_size); 
		eigstates_4_wp.push_back(cur_eigstate);
		empty_swap(cur_eigstate);
	}
	/***------------------------------------------------***/
	std::vector<int> desired_eigvals_no;
	for (size_t i=0; i<INIT_STATE_CONST::eigstate_srno_wp_tun_vect.size(); i++)	{
		desired_eigvals_no.emplace_back(INIT_STATE_CONST::
				eigstate_srno_wp_tun_vect[i]);
	}
	std::vector<double> eigvals_vect = get_all_sp_lines(desired_eigvals_no,DYN_CONST::eigvalfile); 
	/***------------------------------------------------***/
	for (size_t i=0; i<DYN_CONST::stepno_vect_4_anal.size(); i++)	{
		double time = DYN_CONST::stepno_vect_4_anal[i]*DYN_CONST::timestep; 
		std::vector<std::complex<double>> psi_x_anal(wp_size);
		std::vector<double> psi_x_anal_abs(wp_size);
		std::vector<double> psi_x_anal_real(wp_size);
		std::vector<double> psi_x_anal_imag(wp_size);
		for (size_t j=0; j<psi_x_anal.size(); j++)	{
			for (size_t k=0; k<INIT_STATE_CONST::coeff_4_wp_tun_vect.size(); k++)	{
				psi_x_anal[j] += INIT_STATE_CONST::coeff_4_wp_tun_vect[k]*
					(exp(-1.*DYN_CONST::EYE*eigvals_vect[k]*DYN_CONST::
						 twopic*time))*eigstates_4_wp[k][j];
			}
			psi_x_anal_abs[j] = abs(psi_x_anal[j]);
			psi_x_anal_real[j] = psi_x_anal[j].real();
			psi_x_anal_imag[j] = psi_x_anal[j].imag();
		}
		save_binary(psi_x_anal_abs,DYN_CONST::anal_outbinfile_abs,"append");
		save_binary(psi_x_anal_real,DYN_CONST::anal_outbinfile_real,"append");
		save_binary(psi_x_anal_imag,DYN_CONST::anal_outbinfile_imag,"append");
		empty_swap(psi_x_anal);
		empty_swap(psi_x_anal_abs);
		empty_swap(psi_x_anal_real);
		empty_swap(psi_x_anal_imag);
	}
	return EXIT_SUCCESS;
}
