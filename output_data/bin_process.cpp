# include <iostream>
# include <fstream>
# include <vector>
# include <string>
# include <algorithm>


template<typename T>
inline void save_binary(const std::vector<T> & vect,
		const std::string  filename, const std::string& mode)   {
	std::string ext = ".bin";
	if (mode == "append")   {
		std::ofstream outfile(filename+ext, std::ios::binary|std::ios::app);
		if (outfile.is_open())  {
			outfile.write(reinterpret_cast<const char*>(&vect[0]),vect.size()*sizeof(T));
		}
		else    {
			throw "Error opening file for writing!\n";
			exit(1);
		}
	}
}

template<typename T>
inline std::vector<T> rtrn_vec_from_bin(std::vector<T> & dummyvec,
		const std::string& filename,int start_elem_no,int vectsize)   {
	std::ifstream infile;
	std::string ext = ".bin";
	infile.open(filename+ext, std::ios::binary);
	std::vector<T> vect(vectsize);
	infile.seekg((start_elem_no-1)*sizeof(T),std::ios::beg);
	if (infile.is_open())   {
		infile.read(reinterpret_cast<char*>(&vect[0]),vectsize*sizeof(T));
	}
	else    {
		throw "Error opening file for writing!\n";
		exit(1);
	}
	return vect;
}

template<typename T>
inline int get_num_elems(std::vector<T> & dummyvec,
		const std::string& filename)	{
	std::ifstream infile;
	std::string ext = ".bin";
	infile.open(filename+ext,std::ios::binary);
	infile.seekg(0,std::ios::end);
	int filesize = infile.tellg();
	int num_elems = filesize/sizeof(T);
	return num_elems;
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


int main()	{
	std::vector<int> dummyvectint;
	std::vector<double> dummyvectdouble;
	int N_elems = get_num_elems(dummyvectdouble,"Q1_Q10_Q29_closure_case14");
	std::vector<double> vect = rtrn_vec_from_bin(dummyvectdouble,"Q1_Q10_Q29_closure_case14",1,N_elems);
	for_each(vect.begin(),vect.end(),[&](double& elem){std::cout	<<
			std::scientific	<< elem	<< "\n";});
	return 0;
}
