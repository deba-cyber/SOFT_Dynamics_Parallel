# if ! defined (FILEOPS)
# define FILEOPS

# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <algorithm>



/*****************************************************************
 * Function for accumulating data from multiple files with many
 * lines and single/multiple columns .. returning a single vector	
 * all values from all files in a single vector	.................
 * This function is specific to this calculation ................
 * working with the files containing transformation matrices 
 * SO => ANH & ANH_2_PO and prim_dvr_points files ...............
 * ***************************************************************/

std::vector<double> vect_accumulate_all_col(int* mode_id_arr, const std::string& s1, int NDIM);

/////////////////////////////////////////////////////

/****************************************************
 * Function template to free up memory allocated 
 * for vector of any data type	********************
 * **************************************************/


template<typename T>
inline auto empty_swap(std::vector<T> & vec)	{
	std::vector<T>().swap(vec);
	return vec;
}

/////////////////////////////////////////////////////

/****************************************************
 * Function template to return specific queried row *
 * or column elements from passed one dimensional 
 * vector (which is obtained from two dimensional 
 * file data ........................***************
 * **************************************************/

template<typename T>
inline auto get_specific_row_or_col_elements(std::vector<T> & vec_onedim,
		int N_row_or_col, int N_tot_rows, int N_tot_cols, 
		const std::string s1)	
{
	/** Indexing for N_row_or_col is 0-based i.e. if 
	 *	N_row_or_col is 100 , actually 101-th row is desired **/
	int start_cnt;
	int end_cnt;
	int cur_row = 0;
	int cur_id_4_col = N_row_or_col;		// index for iterating over for extracting column 
	if (s1 == "ROW")	{
		std::vector<T> rtrn_sp_row_elems;
		start_cnt = N_row_or_col*N_tot_cols;
		end_cnt = (N_row_or_col+1)*N_tot_cols;
		for (int i=start_cnt; i<end_cnt; i++)	{
			rtrn_sp_row_elems.emplace_back(vec_onedim[i]);
		}
		return rtrn_sp_row_elems;
	}
	else if (s1 == "COLUMN")	{
		std::vector<T> rtrn_sp_col_elems;
		while (cur_row != N_tot_rows)	{
			rtrn_sp_col_elems.emplace_back(vec_onedim[cur_id_4_col]);
			cur_id_4_col += N_tot_cols;
			cur_row +=1;
		}
		return rtrn_sp_col_elems;
	}
}

/*******************************************************
 * Function template to return matrix-vector product ***
 * where the matrix is passed in argument as 1-dim 
 * vector	by reference .............................
 * No of rows and columns are also passed by values 
 * *****************************************************/

template<typename T>
inline std::vector<T> get_mat_vec_pdt(std::vector<T> & mat_as_vect, std::vector<T> & onedim_vect,
		int N_row, int N_col)
{
	std::vector<T> rtrn_matvec_pdt_vect;
	std::vector<T> tmp_process_vect;
	int onedim_vect_size = onedim_vect.size();
	T sum;
	if (N_col != onedim_vect_size)	{
		throw "Matrix-Vector multiplication not possible!";
	}
	else
	{
		for (int i=0; i<N_row; i++)	{
			tmp_process_vect = get_specific_row_or_col_elements(mat_as_vect,i,N_row,N_col,"ROW");
			sum = 0.;
			for (size_t j=0; j<tmp_process_vect.size(); j++)	{
				sum += tmp_process_vect[j]*onedim_vect[j];
			}
			rtrn_matvec_pdt_vect.emplace_back(sum);
		}
	}
	return rtrn_matvec_pdt_vect;
}

/************************************************
 * Function for extracting a specific line from 
 * a file of arbitrary length	.. .............
 * this will be used to get the roots and weights
 * for GAUSS-HERMITE quadrature saved in file 
 * calculated before using numpy.polynomial routine
 * specific to data type .. here double
 * **********************************************/

std::vector<double> get_specific_line(int lineno, const std::string& s1);	

/**************************************************
 * Function for extracting specific lines
 * from a file	.................................
 * string s1 is passed by value filename wo ext.
 * vector lineno_vect contains numbers in increasing
 * order corrsponsing to lines to be extracted	
 * lineno_vect has 1-based indexing	...............
 * returning as one dimensional vector ..........
 * specific to data type .. here double
 * ************************************************/

std::vector<double> get_all_sp_lines(const std::vector<int> & lineno_vect, const std::string& s1);	

/***************************************************************************
 * Function for returning stride array as a vector for any given vector as 
 * input argument which gives number of elements for each dimension ....
 * *************************************************************************/

std::vector<int> GET_STRIDE_ARR_4_ANY(const std::vector<int> & size_vect);	

/********************************************************************
 * function template for appending elements to an existing binary 
 * file ... 
 * data is passed by reference in an one dimensional vector
 * alongwith filename & mode which has to be "append" for this 
 * to work
 * ******************************************************************/

template<typename T>
inline void save_binary(std::vector<T> & vect,
		const std::string  filename, const std::string& mode)	{
	std::string ext = ".bin";
	if (mode == "append")	{
		std::ofstream outfile(filename+ext, std::ios::binary|std::ios::app);
		if (outfile.is_open())	{
			outfile.write(reinterpret_cast<char*>(&vect[0]),vect.size()*sizeof(T));
		}
		else	{
			throw "Error opening file for writing!\n";
			exit(1);
		}	
	}
}

/***********************************************************************
 * function template for extracting specific elements from a binary 
 * file ....
 * start_elem_no gives the element from which extraction starts
 * (1-based indexing)
 * vectsize gives number of elements to be extracted including the 
 * starting element ....................................................
 * *********************************************************************/


template<typename T>
inline std::vector<T> rtrn_vec_from_bin(std::vector<T> & dummyvec,
		const std::string& filename,int start_elem_no,
		int vectsize)	{
	/*****************************************
	 * passing a dummy vector of same data type 
	 * of zero size for keeping template structure
	 * ***************************************/
	std::ifstream infile;
	std::string ext = ".bin";
	infile.open(filename+ext, std::ios::binary);
	std::vector<T> vect(vectsize);
	infile.seekg((start_elem_no-1)*sizeof(T),std::ios::beg);
	if (infile.is_open())	{
		infile.read(reinterpret_cast<char*>(&vect[0]),vectsize*sizeof(T));
	}
	else	{
		throw "Error opening file for writing!\n";
		exit(1);
	}
	return vect;
}

/***************************************************
 * function template to get number of elements in 
 * a binary file ..................................
 * *************************************************/


template<typename T>
inline int get_num_elems(std::vector<T> & dummyvec,
		const std::string& filename)	{
	/*****************************************
	 * passing a dummy vector of same data type 
	 * of size zero for keeping template structure
	 * ***************************************/
	std::ifstream infile;
	std::string ext = ".bin";
	infile.open(filename+ext,std::ios::binary);
	infile.seekg(0,std::ios::end);
	int filesize = infile.tellg();
	int num_elems = filesize/sizeof(T);
	return num_elems;
}

/***----------------------------------------------------------***/
/** Function template to print out 1D std::vector **/

template<typename T>
void printonedimvect(const std::vector<T>& vect)	{
	for (unsigned int i=0; i<vect.size(); i++)	
		std::cout	<< std::scientific	<< vect[i]	<< "\n";
}

/***----------------------------------------------------------***/
/** Function template to print out 2D std::vector<std::vector<T>> **/

template<typename T>
void printtwodimvect(const std::vector<std::vector<T>>& vect)	{
	for (unsigned int i=0; i<vect.size(); i++)	{
		for (unsigned int j=0; j<vect[i].size(); j++)	{
			std::cout	<< std::scientific	<<	vect[i][j]	<<	"	";
		}
		std::cout	<< "\n";
	}
}

/***************************************
 * Function template to save 2D vector	
 * i.e. std::vector<std::vector<T>> to 
 * binary file ...
 * *************************************/

template<typename T>
void save_twodim_vect_binary(const std::vector<std::vector<T>>& twodim_invect,
		const std::string filename, const std::string mode)	{
	std::string ext = ".bin";
	if (mode == "append")	{
		std::ofstream outfile(filename+ext, std::ios::binary|std::ios::app);
		if (outfile.is_open())	{
			for (size_t i=0; i<twodim_invect.size(); i++)	{
				outfile.write(reinterpret_cast<const char*>(&twodim_invect[i][0]),twodim_invect[i].size()*sizeof(T));
			}
		}
		else	{
			throw "Error opening file for writing!\n";
			exit(1);
		}
	}
}

/** Function template to convert 1D vector (std::vector<T>) to 2D 
 * vector std::vector<std::vector<T>> with separation at equal intervals 
 * No of rows is equal to size_t partition **/

template<typename T>
inline std::vector<std::vector<T>> get_twodim_vect_from_vect(const std::vector<T>& invect,
		size_t partition)	{
	size_t N_col = invect.size()/partition;
	std::vector<std::vector<T>> rtrn_multvect;
	size_t start_id;
	size_t end_id;
	for (size_t i=0; i<partition; i++)	{
		start_id = i*N_col;
		end_id = (i+1)*N_col;
		std::vector<T> tmpvect;
		for (size_t j=start_id; j<end_id; j++)	{
			tmpvect.emplace_back(invect[j]);
		}
		rtrn_multvect.push_back(tmpvect);
	}
	return rtrn_multvect;
}

/***----------------------------------------------------------***/
/** Function template to generate 1D grid from lower,upper limits and no of grid points **/

template<typename T>
inline std::vector<T> generate_onedim_grid(T low_range_val,T up_range_val,int N_pts)	{
	double grid_spacing = (up_range_val-low_range_val)/(N_pts-1);
	std::vector<T> onedim_grid_vect(N_pts);
	for (int i=0; i<N_pts; i++)	{
		onedim_grid_vect[i] = low_range_val+i*grid_spacing;
	}
	return onedim_grid_vect;
}

/***************************************************
 * Function for generating line number vector from 
 * number of lines ..
 * *************************************************/

std::vector<int> Get_lineno_vect(int NUMLINES);


/*******************************************************
 * Function template for extracting all lines from
 * file and storing in a multidimensional vector
 * Useful when each line in the file has different 
 * number of elements ..
 * lineno_vect has line numbers to be extracted (1-based
 * indexing)
 * each element of the multidimensional vector has length
 * equal to number of elements in the line in the file
 * passing dummy vector to keep template structure
*******************************************************/

template<typename T>
inline std::vector<std::vector<T>> Get_all_lines_multvect(const std::vector<int>& lineno_vect,
		const std::vector<T>& dummyvec,const std::string& fstring)	{
	std::string fileext = ".dat";
	std::string line;
	std::string filepath = fstring + fileext; 
	std::vector<std::vector<T>> rtrn_multline_multvect;
	std::ifstream myfile(filepath);
	int l_cnt = 0;
	for (size_t i=0; i<lineno_vect.size(); i++)	{
		int cur_lineno = lineno_vect[i];
		std::vector<T> tmp_vect;
		if (myfile.is_open())	{
			while (l_cnt != cur_lineno && std::getline(myfile,line))	{
				++ l_cnt;
			}
			if (l_cnt == cur_lineno)	{
				std::stringstream lstream(line);
				T val;
				while (lstream >> val)	{
					tmp_vect.emplace_back(val);
				}
				rtrn_multline_multvect.push_back(tmp_vect);
			}
		}
	}
	return rtrn_multline_multvect;
}

/** Function template to generate 2D meshgrid from 1D grids **/


template<typename T>
inline std::vector<std::vector<T>> generate_2d_meshgrid(const std::vector<T>& first_coord_vect,
		const std::vector<T>& second_coord_vect)	{
	unsigned int meshgrid_size = first_coord_vect.size()*second_coord_vect.size();
	std::vector<std::vector<T>> twodim_meshgrid(meshgrid_size,std::vector<T>(2,first_coord_vect[0]));
	unsigned int counter = 0;
	for (unsigned int i=0; i<first_coord_vect.size(); i++)	{
		for (unsigned int j=0; j<second_coord_vect.size(); j++)	{
			twodim_meshgrid[counter][0] = first_coord_vect[i];
			twodim_meshgrid[counter][1] = second_coord_vect[j];
			counter++;
		}
	}
	return twodim_meshgrid;
}

# endif
