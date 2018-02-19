#include "TableBase.h"
#include "H5Cpp.h"
#include "predefine.h"
#include <iostream>

template <typename T, size_t N>
TableBase<T, N>::TableBase(std::string Name, T (*f)(Dvec), Svec shape, Dvec low, Dvec high):
_Name(Name), _rank(N), _power_rank(std::pow(2, _rank)), 
_shape(shape), _low(low), _high(high),
_table(_shape), _approximating_function(f) {
	std::cout<<_Name << " dim=" << _rank <<std::endl;
	for(auto i=0; i<_rank; ++i){
		_step.push_back((high[i]-low[i])/(shape[i]-1));
	}
}

template <typename T, size_t N>
T TableBase<T, N>::InterpolateTable(Dvec values){
   Svec start_index;
   Dvec w;
   for(auto i=0; i<_rank; ++i) {
       auto x = (values[i]-_low[i])/_step[i];
       x = std::min(std::max(x, 0.), _shape[i]-2.); // cut at lower and higher bounds bounds
       size_t nx = size_t(std::floor(x));
       double rx = x-nx;
       start_index.push_back(nx);
       w.push_back(rx);
   }
   Svec index(_rank);
   T result{0.};
   for(auto i=0; i<_power_rank; ++i) {
        auto W = 1.0;
        for (auto j=0; j<_rank; ++j) {
            index[j] = start_index[j] + ((i & ( 1 << j )) >> j);
            W *= (index[j]==start_index[j])?(1.-w[j]):w[j];
        }
        result = result + _table(index)*W;
   }
   return result;
}

template <typename T, size_t N>
void TableBase<T, N>::SetTableValue(Svec index, T v){
    _table(index) = v;
}

template <typename T, size_t N>
bool TableBase<T, N>::Save(std::string fname){
	// the is a dumb implementation
	H5::H5File file(fname, H5F_ACC_TRUNC);
	H5::Group group = H5::Group( file.createGroup( "/"+_Name ));
	
	boost::multi_array<double, N> buffer(_shape);
	hsize_t dims[_rank];
	for (auto i=0; i<_rank; ++i) dims[i]=_shape[i];
	H5::DSetCreatPropList proplist{};
	proplist.setChunk(_rank, dims);

	H5::DataSpace dataspace(_rank, dims);
	auto datatype(H5::PredType::NATIVE_DOUBLE);
	
	for(auto comp=0; comp<T::size(); ++comp) {
		for(auto i=0; i<_table.num_elements(); ++i) {
			T item = _table.data()[i];
			buffer.data()[i] = item.get(comp);
		}
		H5::DataSet dataset = file.createDataSet("/"+_Name+"/"+std::to_string(comp), 
													datatype, dataspace, proplist);
		dataset.write(buffer.data(), datatype);
	}
	hdf5_add_scalar_attr(group, "rank", _rank);
	for (auto i=0; i<_rank; ++i){
		hdf5_add_scalar_attr(group, "shape-"+std::to_string(i), _shape[i]);
		hdf5_add_scalar_attr(group, "low-"+std::to_string(i), _low[i]);
		hdf5_add_scalar_attr(group, "high-"+std::to_string(i), _high[i]);
	}
	file.close();
	return true;
}

template <typename T, size_t N>
bool TableBase<T, N>::Load(std::string fname){
	H5::H5File file(fname, H5F_ACC_RDONLY);
	H5::Group group = H5::Group( file.openGroup( "/"+_Name ));
	size_t temp_rank;
	hdf5_read_scalar_attr(group, "rank", temp_rank);
	if (temp_rank != _rank) {
		std::cout << "Table rank does not match" << std::endl;
		file.close();
		return false;
	}
	else{
		std::cout << "Rank compitable, loading table" << std::endl;
		for (auto i=0; i<_rank; ++i){
			hdf5_read_scalar_attr(group, "shape-"+std::to_string(i), _shape[i]);
			hdf5_read_scalar_attr(group, "low-"+std::to_string(i), _low[i]);
			hdf5_read_scalar_attr(group, "high-"+std::to_string(i), _high[i]);
			_step[i] = (_high[i] - _low[i])/(_shape[i]-1.);
		}
		_table.resize(_shape);
		hsize_t dims[_rank];
		for (auto i=0; i<_rank; ++i) dims[i]=_shape[i];
		boost::multi_array<double, N> buffer(_shape);
		H5::DataSpace dataspace(_rank, dims);
		auto datatype(H5::PredType::NATIVE_DOUBLE);
		for(auto comp=0; comp<T::size(); ++comp) {
			H5::DataSet dataset = file.openDataSet("/"+_Name+"/"+std::to_string(comp));
			dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE,
						 dataspace, dataset.getSpace());
			for(auto i=0; i<_table.num_elements(); ++i) {
				_table.data()[i].set(comp, buffer.data()[i]);
			}
		}
		file.close();
	}
	return true;
}

template class TableBase<scalar, 2>;
template class TableBase<scalar, 3>;
template class TableBase<scalar, 4>;
template class TableBase<fourvec, 2>;
template class TableBase<fourvec, 3>;


