#include "StochasticBase.h"
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <thread>

template<size_t N>
StochasticBase<N>::StochasticBase(std::string Name, 
								boost::property_tree::ptree config):
_Name(Name)
{
	std::vector<std::string> strs, slots;
	boost::split(strs, _Name, boost::is_any_of("/"));
	std::string process_name = strs[0];
	std::string quantity_name = strs[1];
	std::cout << __func__ << " " << _Name << std::endl;
	auto tree = config.get_child(process_name+"."+quantity_name);
	std::string allslots = tree.get<std::string>("<xmlattr>.slots");
	boost::split(slots, allslots, boost::is_any_of(",") );

	std::vector<size_t> shape;
	std::vector<double> low, high;
	for(auto & v : slots){
		shape.push_back(tree.get<size_t>("N"+v));
		low.push_back(tree.get<double>("L"+v));
		high.push_back(tree.get<double>("H"+v));
	}
	
	StochasticBase<N>::_ZeroMoment = 
		std::make_shared<TableBase<scalar, N>>(Name+"/scalar", shape, low, high);
	StochasticBase<N>::_FirstMoment = 
		std::make_shared<TableBase<fourvec, N>>(Name+"/vector", shape, low, high);
	StochasticBase<N>::_SecondMoment = 
		std::make_shared<TableBase<tensor, N>>(Name+"/tensor", shape, low, high);
_ZeroMoment->SetTableValue({0, 0}, {10});
}

template<size_t N>
void StochasticBase<N>::save(std::string fname){
	_ZeroMoment->Save(fname);
	_FirstMoment->Save(fname);
	_SecondMoment->Save(fname);
}

template<size_t N>
void StochasticBase<N>::init(void){
	std::cout << "init table" << std::endl;
	auto code = [this](int start, int end) { this->compute(start, end); };
	std::vector<std::thread> threads;
	// ZeroMoment table
	size_t nthreads = 9;
	size_t padding = size_t(std::ceil(_ZeroMoment->length()*1./nthreads));
	for(auto i=0; i<nthreads; ++i) {
		int start = i*padding;
		int end = std::min(padding*(i+1), _ZeroMoment->length());
		threads.push_back( std::thread(code, start, end) );
	}
	for(auto& t : threads) t.join();
}

template<size_t N>
void StochasticBase<N>::compute(int start, int end){
	std::cout << "start= " << start << ", end = " << end << ", computing table" << std::endl;
	std::vector<size_t> index;
	index.resize(N);
	for(auto i=start; i<end; ++i){
		size_t q = i;
		for(int d=N-1; d>=0; d--){
			size_t dim = _ZeroMoment->shape(d);	
			size_t n = q%dim;
			q = q/dim;
			index[d] = n;
		}
		scalar res = calculate_scalar(_ZeroMoment->parameters(index));
		_ZeroMoment->SetTableValue(index, res);
	}
}

template class StochasticBase<2>;
template class StochasticBase<3>;
template class StochasticBase<4>;
