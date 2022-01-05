# Wave packet dynamics routine based on Split-Operator Fourier Transform method : Parallel version with OpenMP

# This is a header only library

# Language = C++ (C++17 dependencies)[Parallelized with OpenMP]

# Overview 

**SOFT_dynamics_Parallel** calculates wave packet dynamics based on 
Split-Operator Fourier Transform method.
(Hence can be used for a rectilinear normal coordinate system 
 and an orthogonal curvilinear coordinate system.
 To use for orthogonal curvilinear coordinate system,
 dynamics routines are to be modified for the diagonal scale 
 factors.)

Wave packet dynamics analysis is done on the fly. Some calculations are specific 
to the particular problem I solved during my PhD, which is tunneling dynamics in double well 
potential in "Catecholate Monoanion".

**SOFT_dynamics_Parallel** also calculates initial vibrational wave packet
from linear combination of eigenstates calculated beforehand for a given Hamiltonian or 
direct product state from quantum of excitations in different 
vibrational modes in a multidimensional vibrational problem. 

# The dynamics code can be used separately from initial state preparation code (separate targets are there for initial state preparation and dynamics in the Makefile)

******************************************************************************************************
******************************************************************************************************

# To run dynamics code with initial state prepared differently, file containing initial state values is to be there in "input_data" folder in binary/text format. 

# In the inputs given in **include/SOFT/constants.hpp** and **src/constants.cpp**[are specific to the problem I have worked on which uses initial state preparation and dynamics in dimension-less normal coordinates. The routines can be used in problems with the kinds of coordinate systems mentioned earlier.]
		
# **include/SOFT** contains the headers. 
# **src** contains the corresponding source files and files with the main function
# **test** is the folder for testing the numerical accuracy of the calculation with testing with wave packet prepared from linear combination of eigenstates and hence analytical time-evolution was possible.

## For general use, necessary routines are in **include/SOFT** and **src** folders.

******************************************************************************************************

# Units :
# Dimension-less normal coordinates were used for the calculation
# Potential energy is in wavenumber unit. 

# For initial state preparation, the following scheme is used (useful for rectilinear normal coordinates in vibrational problems)

***primitive DVR (discrete variable representation, Hermite-DVR/sinc-DVR is used)==> (One dimensional basis functions) ==> PODVR(Potential optimised DVR)***

###---------------------------------------###

** whenever a filename is to be given as input in **include/SOFT/constants.hpp** and/or **src/constants.cpp**
	it should be given without the file extension (data will be processed using routines in 
	**include/SOFT/fileops.hpp** and/or **src/fileops.cpp** which has the necessary extensions in the routines.)

## Input instruction for preparing initial state :

## Nothing to be changed in namespace **FUNDAMENTAL_CONST**

## namespace **INIT_STATE_CONST**

	## All extern const types are to be given in corresponding source file **src/constants.cpp**

	coeff_size_max ==> maximum number of eigenstates taken for initial state preparation when
						linear combination of eigenstates is used for initial state.

	NDIM_4_WP_INIT ==> Number of coordinates in the initial state.

	NDIM_4_DYNAMICS ==> Number of coordinates for dynamics.

	{Dynamics dimension may be less than dimension of the wave function.
		For example, one may be interested in dynamics in (x,y) grid at constant value of 
		z; then grid dimension is direct product of sizes in x,y; but wave function 
		is 3-dimensional (x,y,z)}

	dyn_type ==> Information if dynamics type is **TUNNELING** or **IVR**
				if initial state is prepared from linear combination of eigenstates,
				write **TUNNELING**
				else if initial state is direct product state from certain qunata of 
				excitations in different coordinates, 
				write **IVR**

	eigstate_srno_wp_tun_vect ==> eigenstate indices for linear combination is to be provided.
									(1-based indexing)

	coeff_4_wp_tun_vect	==> corresponding coefficients of eigenstates making up the initial 
							wave packet is to be provided.

	**(Size of all these vectors is equal to NDIM_4_WP_INIT)**

	basis_size_vect		==> basis size for individual coordinates

	mode_id_vect		==> Normal mode indices (1-based indexing)(give in increasing order)

	(basis_size_vect & mode_id_vect elements has one-to-one correspondence.)

	## Following few constants are specific to the scheme used for 
		initial state generation (DVR => 1-dim basis => PODVR) 

	Q_shift_vect		==> Value of the local minimum for one dimensional 
							basis calculation (non-zero value only for coordinates 
							for which Hermite-DVR is used.)

	dvr_type_vect(data type string)		==> For a particular mode, if primitive DVR is "HERMITE", write **HO**
							if primtive DVR is "sinc", write **sinc**

	##	(In the scheme for getting 1-dimensional PODVR basis functions, few transformation matrices are 
			required,following are information about number of rows and columns for different such transformations).

	n_row_so_2_anh_vect		==>	vector containing number of rows for different transformation matrices from 
								primitive DVR ==> 1-dim basis functions

	n_col_so_2_anh_vect		==>	vector containing number of columns for different transformation matrices from
								primitive DVR ==> 1-dim basis functions

	n_row_anh_2_po_vect		==> vector containing number of rows for different transformation matrices from 
								1-dim basis => PODVR function.

	n_col_anh_2_po_vect		==>	vector containing number of columns for different transformation matrices from
								1-dim basis => PODVR function.

	tot_prim_dvr_pts		==> total number of primitive DVR points after summing up for all coordinates in the 
								initial state preparation.

	sodvr_2_anharm_file		==> filename corresponding to transformation matrix for primitive DVR => 1-dim basis 

	anharm_2_podvr_file		==> filename corresponding to transformation matrix for 1-dim function => PODVR function
	
	multidim_eigvect_filestring		==> filename corresponding to eigenvector matrix for multidimensional eigenstate 
										calculation (relevant for initial state prepared from linear combination of
										eigenstates.
	
	IVR_init_quanta_vect		==> vector containing quanta of excitation in different coordinates (relevant for 
									initial state prepared as direct product from quanta of excitation in different
									coordinates.

		(Following two input constant vectors are relevant for multidimensional grid preparation.
		As hinted earlier,dynamics can be done in a dimension lower than dimension of the wave function.
		 generating direct product grid in the dynamics dimension e.g. (x,y,z) .. 
		 when full dimensional direct product grid is prepared ..
		In "constants.hpp" NDIM_4_WP_INIT and NDIM_4_DYNAMICS should be equl in this case and equal to desired dimension.
		Input vectors for lower range, upper range & no of points for multidimensional grid are to be modified. 
		For direct product grid in the dynamics dimension,length of the lower range,upper range and no of points 
		vectors have to be equal to dimension of dynamics.)
	
		**To be changed in "constants.cpp" in namespace DYN_CONST**

		lower range ==> "grid_low_dim_vect";  
		upper range ==> "grid_up_dim_vect";  
		no of points ==>"grid_pt_vect";
	
	In namespace INIT_STATE_CONST in "constants.cpp"
	"meshgrid_4_dyn_mode" normal mode numbers are to be given relevant 
	to the dynamics.
	(For dynamics in direct product grid in dynamics dimension, length of 
	this vector has to be equal to dynamics dimension.)
	Length of vector other_coord_val_vect has to be equal to dynamics dimension
	and for this case, all elements are to be zero (i.e. for N-dimensional 
	direct product grid where `N' is dynamics dimension)

	other_coord_val_vect		==> vector containing values of coordinates 
									which may have been kept fixed.
									length of this vector is equal to NDIM_4_WP_INIT
									if dynamics dimension is equal to dimension of the 
									initial wave packet, all elements of this vector 
									has to be zero.
									if dynamics dimension is lower than dimension 
									of the initial wave packet, values corresponding
									to dynamics coordinates are to be kept zero.
									values corresponding to other coordinates are to 
									be given which are kept fixed for dynamics.

	meshgrid_4_dyn_mode			==> Normal mode indices in which dynamics is carried out.
									(length to be equal to dynamics dimension.
									only modes in the dynamics are to be given in
									increasing order.)
										

		**Generating multidimensional grid when dynamics is carried out in 
		a lower dimensional direct product grid than NDIM_4_WP_INIT
		(e.g. in (x,y,z) generating (x,y) variations with constant z for all)**

		In this case => NDIM_4_WP_INIT should be equal to dimension of the wave function.
					NDIM_4_DYNAMICS should be changed to what desired
					no of dimension for direct product grid generation.
					(i.e. should be changed to 2 when (x,y) grid is prepared
					at constant z)

		Input vectors "grid_low_dim_vect";"grid_up_dim_vect","grid_pt_vect"
		are to be changed .. length of each of which equal to NDIM_4_DYNAMICS 
		in **src/constants.cpp**
	
		Input vector "meshgrid_4_dyn_mode" will contain relevant coordinates
		e.g. (x,y) if those are taken for direct product grid generation.
		Length of vector "other_coord_val_vect" again has to be equal to NDIM_4_WP_INIT
		values of coordinate involved in dynamics has to be kept zero.
		other coordinate values will be changed at the desired constant 
		value.(location of different coordinates in the vector has one to 
		one correspondence to vector "mode_id_vect".)

###---------------------------------------###

## Input instruction for doing wave packet dynamics :

	timestep ==> step-size for dynamics (to be given in femtosecond based on current structure.)

	stepno_vect_4_anal	==> vector containing step numbers for which analytical dynamics is to 
							be carried out (for testing) 

	max_step	==> maximum number of steps upto which dynamics will be carried out.

	twopic		==> 2*pi*speed_of_light(in cm/femtosecond) 
				(this is done for necessary unit consistency as energy unit 
				was in wavenumber.
				For general use, unit of choice can be used and routines in 
				**src/TROTTER_prop.cpp** is to be changed for Split-operator 
				dynamics.)

	EYE		==> imaginary 'i'

	step_4_file		==> every n-th step when evolved wave packet is to be saved to file.

	step_4_restart	==> every n-th step when evolved wave packet is to be saved to file for restart (in case dynamics is terminated midway.)

	eigvalfile		==> file containing multidimensional eigenvalues relevant to dynamics
						(useful for testing with analytical evolution.)

	eigstate_4_anal_dyn_file	==> filename for eigenvectors in the relevant multidimension.
									(filenames are to be in the above style and relevant 
									eigenvectors for **TUNNELING** dynamics will be taken 
									using information in **eigstate_srno_wp_tun_vect**)

	grid_pt_vect		==> vector containing number of points along each dimension.
							(length has to be equal to dynamics dimension.)


	grid_low_lim_vect	==> vector containing values of lower limit along each of the 
							coordinates relevant to the dynamics.

	grid_up_lim_vect	==> vector containing values of upper limit along each of the 
							coordinates relevant to the dynamics.

	freq_vect		==> vector containing frequencies of the coordinates relevant to the
						dynamics in wavenumber unit.

	initfiletype	==> whether potential file and initial wave function file are in binary/dat
						form. (write "bin"/"dat" based on file type)

			** following are some filenames for saving output. Give filenames without extensions.
				By default, saving will be done in binary format**

	potfile			==> name of the file containing potential values in the grid 
						where dynamics will be carried out.

	initfile		==> name of the file containing initial state values.

	restartoutbinfile	==> filename for saving wave packet for restart.

	save_outbinfile		==>		filename for saving wave packet during the dynamics.

	anal_outbinfile_abs	==>		filename for saving absolute value of the wave packet during analytical dynamics.

	anal_outbinfile_real	==> filename for saving real part of the wave packet during analytical dynamics.

	anal_outbinfile_imag	==> filename for saving imaginary part of the wave packet during analytical dynamics.


*************************************************************

## How to run

## Run all executables from **bin** folder

## necessary inputs are to be changed in either **src/constants.cpp** and/or **include/SOFT/constants.hpp**

## Necessary input files are in folder **input_data**

## output files will be saved in **output_data** folder.

## compilation

## Before any run :

	make clean

## For initial state preparation :
	make init.exe (Prepared initial state will be saved in **input_data**)	
	make move (moving executables to **bin** folder)


# For running the dynamics :
	make dyn.exe
	make move (moving executables to "bin" folder)


# Author 

Debabrata Bhattacharyya (debabratab@iisc.ac.in & deba.bhat.90@gmail.com)




