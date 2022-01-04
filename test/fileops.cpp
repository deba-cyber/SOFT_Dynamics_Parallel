# include <SOFT/fileops.hpp>
# include <iostream>
# include <fstream>
# include <sstream>
# include <string>
# include <vector>
# include <numeric>



/*****************************************************************
 * Function for accumulating data from multiple files with many
 * lines and single/multiple columns .. returning a single vector	
 * all values from all files in a single vector	.................
 * This function is specific to this calculation ................
 * working with the files containing transformation matrices 
 * SO => ANH & ANH_2_PO and prim_dvr_points files ...............
 * ***************************************************************/

std::vector<double> vect_accumulate_all_col(int* mode_id_arr, const std::string& s1, int NDIM)	{
	std::string line1;
	std::string ext = ".dat";
	std::vector<double> vect;
	for (int i=0; i<NDIM; i++)	{
		std::string path = s1 + std::to_string(mode_id_arr[i]) + ext; 
		std::ifstream myfile(path);
		if (myfile.is_open())	{
			while (std::getline(myfile,line1))	{
				std::stringstream lstream(line1);
				double val;
				while (lstream >> val)	{
					vect.emplace_back(val);
				}
			}
		}
	}
	return vect;
}

/////////////////////////////////////////////////////


/************************************************
 * Function for extracting a specific line from 
 * a file of arbitrary length	.. .............
 * this will be used to get the roots and weights
 * for GAUSS-HERMITE quadrature saved in file 
 * calculated before using numpy.polynomial routine
 * specific to data type .. here double
 * **********************************************/

std::vector<double> get_specific_line(int lineno, const std::string& s1)	{
	// lineno is given in 1-based indexing //
	// string s1 provides the filename without the extension //
	std::string line1;
	std::string ext = ".dat";
	std::vector<double> rtrn_sp_row;
	int l_cnt = 0;
	std::string path = s1 + ext;
	std::ifstream myfile(path);
	if (myfile.is_open())	{
		while (l_cnt< lineno && std::getline(myfile,line1))	{
			++l_cnt;
		}
		if (l_cnt == lineno)	{
			std::stringstream lstream(line1);
			double val;
			while (lstream >> val)	{
				rtrn_sp_row.emplace_back(val);
			}
		}
	}
	return rtrn_sp_row;
}


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

std::vector<double> get_all_sp_lines(const std::vector<int> & lineno_vect, const std::string& s1)		{
	std::string line1;
	std::string ext = ".dat";
	std::string path = s1 + ext;
	std::vector<double> rtrn_multiple_lines_vect;
	std::ifstream myfile(path);
	int cur_lineno;
	int l_cnt = 0;
	for (size_t i=0; i<lineno_vect.size(); i++)	{
		cur_lineno = lineno_vect[i];
		if (myfile.is_open())	{
			while (l_cnt!= cur_lineno	&& std::getline(myfile,line1))	{
				++l_cnt;
			}
			if (l_cnt == cur_lineno)	{
				std::stringstream lstream(line1);
				double val; 
				while (lstream >> val)	{
					rtrn_multiple_lines_vect.emplace_back(val);
				}
			}
		}
	}
	return rtrn_multiple_lines_vect;
}

/***************************************************************************
 * Function for returning stride array as a vector for any given vector as 
 * input argument which gives number of elements for each dimension ....
 * *************************************************************************/



std::vector<int> GET_STRIDE_ARR_4_ANY(const std::vector<int> & size_vect)	{
	int size = size_vect.size();
	int cur_product;
	int cur_index_pdt = 1;
	for (int j=1; j<size; j++)	{
		cur_index_pdt *= size_vect[j];
	}
	/** Initialising stride array with first element **/
	std::vector<int>stride_arr_init;
	stride_arr_init.emplace_back(cur_index_pdt);
	/** other elements of stride array will be prepared from the first element of the stride array **/
	int cur_index_init = 1;		// initialising current index for generating other elements of stride array 
	while (true)	{
		if (cur_index_init == size-1)	{
			break;
	}
		else 
		{
			cur_product = int(cur_index_pdt/size_vect[cur_index_init]);
			cur_index_init +=1;
			stride_arr_init.emplace_back(cur_product);
			cur_index_pdt = cur_product;
		}
	}
	/////////////////////////////
	int multidim_index = std::accumulate(std::begin(size_vect),std::end(size_vect),1,std::multiplies<double>());  
	std::vector<int> rtrn_all_stride_vect;	// returning as one dimensional vector
	std::vector<int>onedim_index_vect;		// placeholder vector for each integer up tp multidim_index
    /** multidim_index will change for finding each of the 1-d indices in the loop .. 
	 * here it is first initialized **/
	for (int i=1; i<multidim_index+1; i++)	{
		std::vector<int> onedim_index_vect;
		int multidim_index_4_caln = i-1;
		int cur_onedim_index;
		for (size_t j=0; j<stride_arr_init.size(); j++)	{
			cur_onedim_index = int(multidim_index_4_caln/stride_arr_init[j]);
			onedim_index_vect.emplace_back(cur_onedim_index);
			multidim_index_4_caln -= cur_onedim_index*stride_arr_init[j];
		}
		onedim_index_vect.emplace_back(multidim_index_4_caln);
		rtrn_all_stride_vect.insert(end(rtrn_all_stride_vect),begin(onedim_index_vect),end(onedim_index_vect));
	}
	return rtrn_all_stride_vect;
}



