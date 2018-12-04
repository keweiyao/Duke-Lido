#include <vector>
#include <string>
#include <iostream>
#include <boost/multi_array.hpp>
#include <H5Cpp.h>

// helper function for read hdf5 scalar attributes
template <typename T> inline const H5::PredType& type();
template <> inline const H5::PredType& type<size_t>() { return H5::PredType::NATIVE_HSIZE; }
template <> inline const H5::PredType& type<double>() { return H5::PredType::NATIVE_DOUBLE; }
template <typename T>
void hdf5_read_scalar_attr(
  const H5::Group& gp, const std::string& name, T& value) {
  const auto& datatype = type<T>();
  auto attr = gp.openAttribute(name.c_str());
  attr.read(datatype, &value);
}

// function that load a table
template<int rank>
bool load(std::string fname, std::string groupname, int component_index, boost::multi_array<double, rank>& Table){ // filename, groupname, and which component (0 for scalar, 0-3 for vector and 0-15 for tensor)
	// Openfile
	H5::H5File file(fname, H5F_ACC_RDONLY);
	// OpenDataSett
	H5::Group group = H5::Group( file.openGroup(groupname));
	
	size_t temp_rank;
	hdf5_read_scalar_attr(group, "rank", temp_rank);
	if (temp_rank != rank) {	// Make sure the rank is the same as expected
		std::cout<< "Table rank does not match"<< std::endl;
		file.close();
		return false;
	}
	else{	// If rank is as expected, load the table
		hsize_t shape[rank];
		double low[rank], high[rank], step[rank];
		std::cout<< "Rank compitable, loading table"<< std::endl;
		for (auto i=0; i<rank; ++i){ // For each dimension, read its low, high, and number of grids, and calculate the step
			size_t tempN;
			hdf5_read_scalar_attr(group, "shape-"+std::to_string(i), tempN);
			shape[i] = tempN;
			hdf5_read_scalar_attr(group, "low-"+std::to_string(i), low[i]);
			hdf5_read_scalar_attr(group, "high-"+std::to_string(i), high[i]);
			step[i] = (high[i] - low[i])/(shape[i]-1.);
		}
		
		// creat data table
		std::vector<size_t> dims;
		dims.resize(3);
		for (auto i=0; i<rank; ++i) {
			dims[i] = shape[i];
		}
		Table.resize(dims);
		// Let H5 know how large the data we want to load
		H5::DataSpace dataspace(rank, shape);
		// Let H5 know what datatype we want to load (native double type)
		auto datatype(H5::PredType::NATIVE_DOUBLE);
		
		
		H5::DataSet dataset = file.openDataSet("/"+groupname+"/" 
												+std::to_string(component_index));
		dataset.read(Table.data(), H5::PredType::NATIVE_DOUBLE,
					 dataspace, dataset.getSpace());

		file.close();
	}
	return true;
}

int main(){
	const size_t rank=2; //rank=2 because it is a table of variable E and T
	boost::multi_array<double, rank> Table;
	bool status = load<rank>("./table.h5", "/Boltzmann/cg2cg/rate/tensor", 5, Table);
	if (status){
		std::cout << "load successful" << std::endl;
		std::cout << Table[0][0] << " " << Table[0][1] << std::endl;
		return 0;
	}
	return -1;
}
